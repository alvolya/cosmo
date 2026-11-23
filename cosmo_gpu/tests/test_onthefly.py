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


@pytest.mark.parametrize(
    ("nucleus", "space", "ham", "parity", "hw"),
    [("24Mg", "sd", "usdb", 1, None), ("20Ne", "spsdpf", "fsu9", -1, (1, 1))],
)
def test_compute_levels_otf_equals_stored(nucleus, space, ham, parity, hw):
    from cosmo_gpu.cli import compute_levels
    from cosmo_gpu.onthefly import OnTheFlyOperator

    kw = {"nev": 6, "parity": parity, "hw": hw, "tol": 1e-10}
    a = compute_levels(nucleus, space, ham, matvec="stored", **kw)
    b = compute_levels(nucleus, space, ham, matvec="otf", **kw)
    assert isinstance(b.H, OnTheFlyOperator)
    assert "build_tables" in b.timings
    np.testing.assert_allclose(b.energies, a.energies, atol=1e-10)
    np.testing.assert_allclose(b.spins, a.spins, atol=1e-8)
    X = np.random.default_rng(8).standard_normal((2, a.basis.dim))
    assert _rel(b.matvec(X), a.matvec(X)) < 1e-13


def test_lanczos_keep_and_history():
    from cosmo_gpu.backend import NumpyBackend
    from cosmo_gpu.lanczos import block_lanczos

    s = _system("25Mg")
    calls: list[tuple[int, np.ndarray, np.ndarray]] = []
    res = block_lanczos(
        lambda X: np.ascontiguousarray((s.H @ X.T).T),
        s.basis.dim,
        5,
        NumpyBackend(),
        block=4,
        ncv=40,
        n_keep=15,
        tol=1e-10,
        history=lambda n, th, r: calls.append((n, th, r)),
    )
    assert res.converged
    assert calls
    assert calls[-1][0] == res.n_matvec
    assert np.all(np.diff([c[0] for c in calls]) > 0)
    np.testing.assert_allclose(calls[-1][1], res.eigenvalues)
