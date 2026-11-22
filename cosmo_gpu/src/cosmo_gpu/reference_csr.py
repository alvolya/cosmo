"""Stored sparse H and J+ built directly from bitstrings: the stage-1 reference path.

Never used at scale. Its role is to be the oracle for the factorized on-the-fly kernels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import sqrt

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

    diag: NDArray[np.float64] = field(default_factory=lambda: _EMPTYV)
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


def build_operator(dst: MBasis, src: MBasis, terms: OperatorTerms) -> sp.csr_matrix:
    """Assemble the stored matrix of ``terms`` from basis ``src`` to basis ``dst``."""
    if dst.protons is not src.protons or dst.neutrons is not src.neutrons:
        msg = "source and destination bases must share their species bases"
        raise ValueError(msg)
    P, N = src.protons, src.neutrons
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
    lookup = (P.sorted_dets, P.det_sub, P.det_local, N.sorted_dets, N.det_sub, N.det_local)
    indptr, indices, data = build_csc(
        src.dim,
        src_t,
        dst_t,
        lookup,
        np.ascontiguousarray(terms.diag, dtype=np.float64),
        np.ascontiguousarray(terms.ob_p, dtype=np.int32),
        np.ascontiguousarray(terms.ob_p_v, dtype=np.float64),
        np.ascontiguousarray(terms.ob_n, dtype=np.int32),
        np.ascontiguousarray(terms.ob_n_v, dtype=np.float64),
        np.ascontiguousarray(terms.tb_p, dtype=np.int32),
        np.ascontiguousarray(terms.tb_p_v, dtype=np.float64),
        np.ascontiguousarray(terms.tb_n, dtype=np.int32),
        np.ascontiguousarray(terms.tb_n_v, dtype=np.float64),
        np.ascontiguousarray(terms.tb_pn, dtype=np.int32),
        np.ascontiguousarray(terms.tb_pn_v, dtype=np.float64),
    )
    mat = sp.csc_matrix((data, indices, indptr), shape=(dst.dim, src.dim))
    return sp.csr_matrix(mat)


def hamiltonian_csr(basis: MBasis, inter: MSchemeInteraction) -> sp.csr_matrix:
    """Build the stored Hamiltonian on ``basis``."""
    terms = OperatorTerms(
        diag=inter.spe,
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
    up = build_mbasis(basis.space, basis.protons, basis.neutrons, basis.m2 + 2, basis.parity)
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
