from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import load_interaction
from cosmo_gpu.operator_tables import (
    expand_to_csr,
    hamiltonian_tables,
    jplus_tables,
    make_moves,
)
from cosmo_gpu.reference_csr import hamiltonian_csr, jplus_csr

from .systems import build_system

SNT = Path(__file__).parent / "data" / "kshell" / "snt"
TOL = 1e-13  # relative to the largest matrix element (summation order differs)


def _max_diff(a: sp.spmatrix, b: sp.spmatrix) -> float:
    """Largest elementwise difference relative to the largest |element| of ``b``."""
    b = sp.csr_matrix(b)
    d = abs(sp.csr_matrix(a) - b)
    scale = max(1.0, float(abs(b).max()))
    return float(d.max()) / scale if d.nnz else 0.0


def _systems():
    yield "20Ne", build_system("20Ne")
    yield "25Mg", build_system("25Mg")
    yield "20Ne_1hw-", build_system("20Ne", "spsdpf", "fsu9", parity=-1, hw=(1, 1))
    yield "33Al_1hw-", build_system("33Al", "spsdpf", "fsu9", parity=-1, hw=(1, 1))


CASES = ["20Ne", "25Mg", "20Ne_1hw-", "33Al_1hw-"]


def _system(name: str):
    return dict(_systems())[name]


@pytest.mark.parametrize("name", CASES)
def test_hamiltonian_tables_reproduce_stored_h(name: str):
    s = _system(name)
    tables = hamiltonian_tables(s.basis, s.inter)
    assert _max_diff(expand_to_csr(tables), s.H) < TOL


def test_hamiltonian_tables_from_snt_44ti():
    space, zv, nv, inter = load_interaction("44Ti", None, None, SNT / "kb3g_from_int.snt")
    basis = build_basis(space, zv, nv, parity=1)
    H = hamiltonian_csr(basis, inter)
    assert _max_diff(expand_to_csr(hamiltonian_tables(basis, inter)), H) < TOL


@pytest.mark.parametrize("name", ["20Ne", "25Mg", "33Al_1hw-"])
def test_jplus_tables_reproduce_stored_jplus(name: str):
    s = _system(name)
    tables, up = jplus_tables(s.basis)
    ref, up_ref = jplus_csr(s.basis)
    assert up.dim == up_ref.dim
    assert _max_diff(expand_to_csr(tables), ref) < TOL


def test_task_list_is_well_formed():
    s = build_system("24Mg")
    t = hamiltonian_tables(s.basis, s.inter)
    assert t.task_ptr[0] == 0
    assert t.task_ptr[-1] == t.n_tasks
    assert np.all(np.diff(t.task_ptr) >= 0)
    # every destination block has its diagonal pp and nn tasks (single-particle energies)
    for B in range(s.basis.n_blocks):
        lo, hi = t.task_ptr[B], t.task_ptr[B + 1]
        srcs = t.task_src_block[lo:hi]
        kinds = t.task_kind[lo:hi]
        assert np.any((srcs == B) & (kinds == 1))
        assert np.any((srcs == B) & (kinds == 2))
    # labels fit in 15 bits and table ids are in range
    assert np.all((t.jump_op.astype(np.int32) & 0x7FFF) >= 0)
    n_jump_tables = len(t.jump_table_start)
    assert np.all(t.task_pjump < n_jump_tables)
    assert np.all(t.task_njump < n_jump_tables)


def test_moves_classes_and_labels(sd_space):
    pairs = np.array([(1, 0), (2, 1), (0, 1), (3, 3)])
    mv = make_moves(sd_space, pairs)
    keys = mv.class_keys
    for c, a, k, q in zip(mv.c, mv.a, mv.cls, mv.label, strict=True):
        d2m = int(sd_space.morb_m2[c] - sd_space.morb_m2[a])
        assert keys[k][0] == d2m
        assert 0 <= q < mv.class_size(k)


def test_tables_scale_with_species_dimensions():
    """Tables of 24Mg sd are much smaller than its stored H."""
    s = build_system("24Mg")
    t = hamiltonian_tables(s.basis, s.inter)
    stored = s.H.nnz * 12  # 8-byte value + 4-byte column
    assert sum(t.nbytes().values()) < stored / 5
