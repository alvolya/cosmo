from __future__ import annotations

import numpy as np
import pytest

from cosmo_gpu.onthefly import hamiltonian_operator, j_expectation_otf, jplus_operator
from cosmo_gpu.reference_csr import j_expectation, jplus_csr

from .systems import build_system

CASES = {
    "25Mg": ("25Mg", "sd", "usdb", None, None),
    "20Ne_1hw-": ("20Ne", "spsdpf", "fsu9", -1, (1, 1)),
    "33Al_1hw-": ("33Al", "spsdpf", "fsu9", -1, (1, 1)),
}


def _system(name: str):
    nucleus, space, ham, parity, hw = CASES[name]
    return build_system(nucleus, space, ham, parity=parity, hw=hw)


def _rel(Y: np.ndarray, Yref: np.ndarray) -> float:
    return float(np.max(np.abs(Y - Yref)) / np.max(np.abs(Yref)))


@pytest.mark.parametrize("name", list(CASES))
def test_otf_matches_stored_h(name: str):
    s = _system(name)
    X = np.random.default_rng(5).standard_normal((4, s.basis.dim))
    Y = hamiltonian_operator(s.basis, s.inter).matvec(X)
    assert _rel(Y, np.asarray((s.H @ X.T).T)) < 1e-13


@pytest.mark.parametrize("name", ["25Mg", "33Al_1hw-"])
def test_otf_jplus_matches_stored(name: str):
    s = _system(name)
    X = np.random.default_rng(6).standard_normal((2, s.basis.dim))
    ref, _ = jplus_csr(s.basis)
    op = jplus_operator(s.basis)
    assert op.shape == ref.shape
    assert _rel(op.matvec(X), np.asarray((ref @ X.T).T)) < 1e-13
    np.testing.assert_allclose(
        j_expectation_otf(X, op, s.basis.m2), j_expectation(X, ref, s.basis.m2), rtol=1e-12
    )


def test_otf_is_deterministic():
    s = _system("33Al_1hw-")
    op = hamiltonian_operator(s.basis, s.inter)
    X = np.random.default_rng(7).standard_normal((3, s.basis.dim))
    assert np.array_equal(op.matvec(X), op.matvec(X))


def test_otf_shape_errors():
    s = _system("25Mg")
    op = hamiltonian_operator(s.basis, s.inter)
    with pytest.raises(ValueError, match="must be"):
        op.apply(np.zeros((2, 5)), np.zeros((2, s.basis.dim)))
