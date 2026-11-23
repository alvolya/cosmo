import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.lanczos import block_lanczos


def _random_symmetric(n, density, seed):
    A = sp.random(n, n, density=density, random_state=seed, format="csr")
    return (A + A.T) * 0.5 + sp.diags(np.linspace(-5, 5, n))


@pytest.mark.parametrize(("block", "nev"), [(1, 4), (4, 6), (3, 10)])
def test_matches_eigsh(block, nev):
    H = _random_symmetric(2000, 0.002, 7)
    res = block_lanczos(
        lambda X: (H @ X.T).T, H.shape[0], nev, NumpyBackend(), block=block, tol=1e-11
    )
    assert res.converged
    ref = np.sort(spla.eigsh(H, k=nev, which="SA", tol=1e-13)[0])
    np.testing.assert_allclose(res.eigenvalues, ref, atol=1e-10)
    Y = res.eigenvectors
    true_res = np.linalg.norm((H @ Y.T).T - res.eigenvalues[:, None] * Y, axis=1)
    assert np.all(true_res < 1e-9)
    np.testing.assert_allclose(Y @ Y.T, np.eye(nev), atol=1e-10)


def test_degenerate_spectrum():
    d = np.repeat(np.arange(300.0), 3)  # every eigenvalue triply degenerate
    H = sp.diags(d) + sp.random(900, 900, density=0.0, random_state=1)
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.standard_normal((900, 900)))
    Hd = Q @ H.toarray() @ Q.T
    res = block_lanczos(lambda X: (Hd @ X.T).T, 900, 6, NumpyBackend(), block=4, tol=1e-10)
    assert res.converged
    np.testing.assert_allclose(res.eigenvalues, [0, 0, 0, 1, 1, 1], atol=1e-9)


def test_too_small_raises():
    # nev = 4, b = 4 needs ncv >= 12 and dim > ncv + b
    with pytest.raises(ValueError, match="dense"):
        block_lanczos(lambda X: X, 14, 4, NumpyBackend())


def test_default_ncv_is_capped_by_dimension():
    rng = np.random.default_rng(3)
    A = rng.standard_normal((20, 20))
    A = A + A.T
    res = block_lanczos(lambda X: X @ A, 20, 4, NumpyBackend())
    np.testing.assert_allclose(res.eigenvalues, np.linalg.eigvalsh(A)[:4], atol=1e-8)
