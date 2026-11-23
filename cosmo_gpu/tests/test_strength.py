from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as sla

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.operators import e2, m1, tensor_matrix
from cosmo_gpu.strength import lanczos_strength

from .systems import build_system

TOL = 1e-10


def _ground_state_pivot(name: str, op_name: str):
    s = build_system(name)
    w, V = sla.eigh(s.H.toarray()) if s.basis.dim <= 2000 else (None, None)
    if V is None:
        from cosmo_gpu.lanczos import block_lanczos

        res = block_lanczos(
            lambda X: np.ascontiguousarray((s.H @ X.T).T), s.basis.dim, 1, NumpyBackend(), tol=1e-12
        )
        psi = res.eigenvectors[0]
    else:
        psi = V[:, 0]
    op = e2(s.space, 1.5, 0.5) if op_name == "E2" else m1(s.space, 1.0, 5.5857, 0.0, -3.82608)
    pivot = tensor_matrix(op, 0, s.basis, s.basis) @ psi
    return s, np.ascontiguousarray(pivot[None, :])


def _matvec(H):
    return lambda X: np.ascontiguousarray((H @ X.T).T)


@pytest.mark.parametrize(("name", "op"), [("20Ne", "E2"), ("24Mg", "M1")])
@pytest.mark.parametrize("reorth", [False, True])
def test_total_strength_and_moments(name: str, op: str, reorth: bool):
    s, p = _ground_state_pivot(name, op)
    n = 12
    st = lanczos_strength(_matvec(s.H), p, n, NumpyBackend(), reorthogonalize=reorth)
    assert st.weights.sum() == pytest.approx(float((p @ p.T)[0, 0]), rel=TOL)
    # moments about the pivot's mean energy keep the powers well conditioned
    e0 = st.moment(1) / st.norm2
    x = p[0].copy()
    hx = [x]
    for _ in range(2 * st.n_steps - 1):
        hx.append(s.H @ hx[-1] - e0 * hx[-1])
    for k in range(2 * st.n_steps):
        direct = float(hx[k // 2] @ hx[k - k // 2])
        lanczos = float(np.sum(st.weights * (st.energies - e0) ** k))
        assert lanczos == pytest.approx(direct, rel=1e-9, abs=1e-9 * st.norm2), k


def test_full_lanczos_reproduces_dense_distribution():
    s, p = _ground_state_pivot("20Ne", "E2")
    H = s.H.toarray()
    w, V = sla.eigh(H)
    amp2 = (V.T @ p[0]) ** 2
    st = lanczos_strength(_matvec(s.H), p, s.basis.dim, NumpyBackend(), reorthogonalize=True)
    keep = amp2 > 1e-12 * amp2.sum()
    dense_e, dense_w = w[keep], amp2[keep]
    lanc = st.weights > 1e-12 * st.weights.sum()
    np.testing.assert_allclose(st.energies[lanc], dense_e, atol=1e-8)
    np.testing.assert_allclose(st.weights[lanc], dense_w, rtol=1e-7, atol=1e-12)


def test_continued_fraction_matches_resolvent():
    s, p = _ground_state_pivot("20Ne", "E2")
    st = lanczos_strength(_matvec(s.H), p, s.basis.dim, NumpyBackend(), reorthogonalize=True)
    z = -30.0 + 0.5j
    H = s.H.toarray()
    direct = p[0] @ np.linalg.solve(z * np.eye(H.shape[0]) - H, p[0])
    assert st.continued_fraction(z) == pytest.approx(direct, rel=1e-10)


def test_smoothed_integrates_to_total():
    s, p = _ground_state_pivot("20Ne", "E2")
    st = lanczos_strength(_matvec(s.H), p, 30, NumpyBackend())
    grid = np.linspace(st.energies.min() - 20, st.energies.max() + 20, 20001)
    area = np.trapezoid(st.smoothed(grid, 0.5, "gaussian"), grid)
    assert area == pytest.approx(st.norm2, rel=1e-6)
