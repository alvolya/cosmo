"""Stored sparse H and J+ built directly from bitstrings: the stage-1 reference path.

Never used at scale. Its role is to be the oracle for the factorized on-the-fly kernels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt

import numba
import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray

from cosmo_gpu.basis import MBasis, build_mbasis
from cosmo_gpu.interaction import MSchemeInteraction
from cosmo_gpu.reference_kernels import build_csc

_EMPTY4 = np.zeros((0, 4), dtype=np.int32)
_EMPTY2 = np.zeros((0, 2), dtype=np.int32)
_EMPTYV = np.zeros(0, dtype=np.float64)


@dataclass(frozen=True)
class OperatorTerms:
    """Second-quantized terms of an operator in pn m-scheme form (see ``reference_kernels``)."""

    diag_p: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    diag_n: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    ob_p: NDArray[np.int32] = field(default_factory=lambda: _EMPTY2)
    ob_p_v: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    ob_n: NDArray[np.int32] = field(default_factory=lambda: _EMPTY2)
    ob_n_v: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    tb_p: NDArray[np.int32] = field(default_factory=lambda: _EMPTY4)
    tb_p_v: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    tb_n: NDArray[np.int32] = field(default_factory=lambda: _EMPTY4)
    tb_n_v: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
    tb_pn: NDArray[np.int32] = field(default_factory=lambda: _EMPTY4)
    tb_pn_v: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)


IndexedTerms = tuple[NDArray[np.int64], NDArray[np.int32], NDArray[np.float64]]


def _index_terms(
    keys: NDArray[np.int64], created: NDArray[np.int32], vals: NDArray[np.float64], n_keys: int
) -> IndexedTerms:
    """Group terms into CSR over their annihilator key: ``(ptr, created, values)``."""
    order = np.argsort(keys, kind="stable")
    ptr = np.zeros(n_keys + 1, dtype=np.int64)
    np.add.at(ptr, keys + 1, 1)
    return (
        np.cumsum(ptr),
        np.ascontiguousarray(created[order], dtype=np.int32),
        np.ascontiguousarray(vals[order], dtype=np.float64),
    )


def _one_body(rows: NDArray[np.int32], vals: NDArray[np.float64], nm: int) -> IndexedTerms:
    rows = rows.reshape(-1, 2)
    return _index_terms(rows[:, 1].astype(np.int64), rows[:, 0], vals, nm)


def _two_body(rows: NDArray[np.int32], vals: NDArray[np.float64], nm: int) -> IndexedTerms:
    """Like-particle rows ``(c1, c2, a1, a2)`` or pn rows ``(cp, ap, cn, an)`` -> indexed."""
    rows = rows.reshape(-1, 4)
    return _index_terms(
        rows[:, 2].astype(np.int64) * nm + rows[:, 3], rows[:, [0, 1]], vals, nm * nm
    )


def _pn(rows: NDArray[np.int32], vals: NDArray[np.float64], nm: int) -> IndexedTerms:
    rows = rows.reshape(-1, 4)
    return _index_terms(
        rows[:, 1].astype(np.int64) * nm + rows[:, 3], rows[:, [0, 2]], vals, nm * nm
    )


def _max_bucket(t: IndexedTerms) -> int:
    return int(np.diff(t[0]).max(initial=0))


def build_operator(dst: MBasis, src: MBasis, terms: OperatorTerms) -> sp.csr_matrix:
    """Assemble the stored matrix of ``terms`` from basis ``src`` to basis ``dst``."""
    if terms.diag_p.shape != terms.diag_n.shape:
        msg = "diag_p and diag_n must both be given or both be empty"
        raise ValueError(msg)
    if dst.protons is not src.protons or dst.neutrons is not src.neutrons:
        msg = "source and destination bases must share their species bases"
        raise ValueError(msg)
    P, N = src.protons, src.neutrons
    nm = src.space.n_morbits
    src_t = (
        src.block_offset,
        src.block_psub,
        src.block_nsub,
        P.sub_start,
        N.sub_start,
        P.dets,
        N.dets,
    )
    dst_t = (dst.block_of, dst.block_offset)
    lookup = (
        P.sorted_dets,
        P.det_sub,
        P.det_local,
        N.sorted_dets,
        N.det_sub,
        N.det_local,
        N.sub_start,
    )
    ob_p = _one_body(terms.ob_p, terms.ob_p_v, nm)
    ob_n = _one_body(terms.ob_n, terms.ob_n_v, nm)
    tb_p = _two_body(terms.tb_p, terms.tb_p_v, nm)
    tb_n = _two_body(terms.tb_n, terms.tb_n_v, nm)
    tb_pn = _pn(terms.tb_pn, terms.tb_pn_v, nm)
    z, n = P.n_particles, N.n_particles
    buf = (
        1
        + z * _max_bucket(ob_p)
        + n * _max_bucket(ob_n)
        + z * (z - 1) // 2 * _max_bucket(tb_p)
        + n * (n - 1) // 2 * _max_bucket(tb_n)
        + z * n * _max_bucket(tb_pn)
    )
    n_chunks = max(1, min(src.dim, 64 * numba.get_num_threads()))
    indptr, indices, data = build_csc(
        src.dim,
        nm,
        buf,
        n_chunks,
        src_t,
        dst_t,
        lookup,
        np.ascontiguousarray(terms.diag_p, dtype=np.float64),
        np.ascontiguousarray(terms.diag_n, dtype=np.float64),
        ob_p,
        ob_n,
        tb_p,
        tb_n,
        tb_pn,
    )
    mat = sp.csc_matrix((data, indices, indptr), shape=(dst.dim, src.dim))
    return sp.csr_matrix(mat)


def hamiltonian_csr(basis: MBasis, inter: MSchemeInteraction) -> sp.csr_matrix:
    """Build the stored Hamiltonian on ``basis``."""
    terms = OperatorTerms(
        diag_p=inter.spe_p,
        diag_n=inter.spe_n,
        tb_p=inter.pp,
        tb_p_v=inter.pp_v,
        tb_n=inter.nn,
        tb_n_v=inter.nn_v,
        tb_pn=inter.pn,
        tb_pn_v=inter.pn_v,
    )
    return build_operator(basis, basis, terms)


def jplus_terms(basis: MBasis) -> OperatorTerms:
    """One-body J+ = sum sqrt(j(j+1) - m(m+1)) a+_{m+1} a_m, for both species."""
    space = basis.space
    rows: list[tuple[int, int]] = []
    vals: list[float] = []
    for k in range(space.n_morbits - 1):
        same_orbit = space.morb_orbit[k] == space.morb_orbit[k + 1]
        if same_orbit:
            j2, m2 = int(space.morb_j2[k]), int(space.morb_m2[k])
            rows.append((k + 1, k))
            vals.append(sqrt((j2 * (j2 + 2) - m2 * (m2 + 2)) / 4.0))
    ob = np.array(rows, dtype=np.int32).reshape(-1, 2)
    v = np.array(vals, dtype=np.float64)
    return OperatorTerms(ob_p=ob, ob_p_v=v, ob_n=ob, ob_n_v=v)


def jplus_csr(basis: MBasis) -> tuple[sp.csr_matrix, MBasis]:
    """J+ from ``basis`` (total 2M) to the basis with total 2M + 2, which is also returned."""
    up = build_mbasis(
        basis.space, basis.protons, basis.neutrons, basis.m2 + 2, basis.parity, basis.truncation
    )
    return build_operator(up, basis, jplus_terms(basis)), up


def j_from_jsq(jsq: NDArray[np.float64]) -> NDArray[np.float64]:
    """Invert ``<J^2> = J (J + 1)``."""
    return 0.5 * (np.sqrt(1.0 + 4.0 * np.maximum(jsq, 0.0)) - 1.0)


def j_expectation(
    vectors: NDArray[np.float64], jplus: sp.csr_matrix, m2: int
) -> NDArray[np.float64]:
    """``J`` of each row of ``vectors`` from ``<J^2> = |J+ psi|^2 + M (M + 1)``."""
    M = m2 / 2.0
    jp = jplus @ vectors.T
    jsq = np.einsum("ij,ij->j", jp, jp) / np.einsum("ij,ij->i", vectors, vectors) + M * (M + 1)
    return j_from_jsq(np.asarray(jsq, dtype=np.float64))
