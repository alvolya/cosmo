"""Block thick-restart Lanczos for the lowest eigenpairs of a real symmetric operator.

The operator is any callable mapping a block ``(b, D)`` to ``H`` applied to each row, so the
stored reference matrix and the on-the-fly kernels are interchangeable. Vectors stay in the
backend; only small projected matrices come to the host.

Each step projects ``H V_j`` on the whole basis (full reorthogonalization, applied twice). The
projections fill the projected matrix ``T`` directly, so after a thick restart (which keeps
Ritz vectors and appends the residual block) ``T`` needs no special arrowhead bookkeeping.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla
from numpy.typing import NDArray

from cosmo_gpu.backend import Array, VectorBackend

MatVec = Callable[[Array], Array]

_RANK_TOL = 1e-24  # relative eigenvalue of the Gram matrix below which a direction is lost


@dataclass
class LanczosResult:
    """Eigenvalues, eigenvectors (rows, in the backend) and convergence data."""

    eigenvalues: NDArray[np.float64]
    eigenvectors: Array
    residuals: NDArray[np.float64]
    converged: bool
    n_matvec: int
    n_restart: int


def _orthonormalize(
    w: Array, basis: Array, k: int, backend: VectorBackend, seed: int
) -> tuple[Array, NDArray[np.float64]]:
    """Return ``(Q, R)`` with ``w = R @ Q``, ``Q`` orthonormal rows, orthogonal to ``basis[:k]``.

    Directions of ``w`` that are numerically zero are replaced by random vectors (with zero
    coefficient in ``R``), so the block size never shrinks.
    """
    b = w.shape[0]
    R_total = np.eye(b)
    q = w
    for sweep in range(2):
        G = backend.inner(q, q)
        lam, U = sla.eigh(G)
        lam_max = max(float(lam[-1]), 0.0)
        bad = lam <= _RANK_TOL * max(lam_max, 1e-300)
        if sweep == 0 and np.any(bad):
            # replace lost directions by fresh random vectors orthogonal to everything so far
            nbad = int(bad.sum())
            fresh = backend.random(nbad, q.shape[1], seed)
            good = backend.combine((U[:, ~bad] / np.sqrt(lam[~bad])).T, q)
            for _ in range(2):
                if k > 0:
                    backend.subtract_combine(fresh, backend.inner(fresh, basis[:k]), basis[:k])
                backend.subtract_combine(fresh, backend.inner(fresh, good), good)
            Gf = backend.inner(fresh, fresh)
            lf, Uf = sla.eigh(Gf)
            fresh = backend.combine((Uf / np.sqrt(lf)).T, fresh)
            R = np.zeros((b, b))
            R[:, : b - nbad] = U[:, ~bad] * np.sqrt(lam[~bad])
            q = backend.from_host(np.concatenate([backend.to_host(good), backend.to_host(fresh)]))
        else:
            q = backend.combine((U / np.sqrt(lam)).T, q)
            R = U * np.sqrt(lam)
        R_total = R_total @ R
    return q, R_total


def _thick_restart(
    V: Array,
    T: NDArray[np.float64],
    S: NDArray[np.float64],
    theta: NDArray[np.float64],
    coupling: NDArray[np.float64],
    q: Array,
    k: int,
    backend: VectorBackend,
) -> None:
    """Keep the Ritz vectors ``S.T @ V[:k]`` and append the residual block ``q``, in place."""
    p, b = theta.shape[0], q.shape[0]
    V[:p] = backend.combine(np.ascontiguousarray(S.T), V[:k])
    V[p : p + b] = q
    T[:] = 0.0
    T[:p, :p] = np.diag(theta)
    T[p : p + b, :p] = coupling
    T[:p, p : p + b] = coupling.T


def block_lanczos(
    matvec: MatVec,
    dim: int,
    nev: int,
    backend: VectorBackend,
    block: int = 4,
    ncv: int | None = None,
    tol: float = 1e-8,
    max_restarts: int = 1000,
    seed: int = 12345,
    start: Array | None = None,
) -> LanczosResult:
    """Lowest ``nev`` eigenpairs of ``matvec`` by block thick-restart Lanczos.

    Parameters
    ----------
    matvec
        Maps a ``(b, dim)`` block to ``H`` applied to each row.
    dim
        Dimension of the space.
    nev
        Number of lowest eigenpairs wanted.
    backend
        Vector backend holding all ``dim``-sized arrays.
    block
        Block size ``b``.
    ncv
        Maximum basis size before a restart; default ``max(6 nev, nev + 6 b)``, capped below
        ``dim - b``. A larger subspace means fewer restarts: on 33Al FSU9 1ħω (10 states,
        b = 4) it cut matvecs from 3356 at ``max(2 nev + 2 b, nev + 6 b)`` to 1260.
    tol
        Convergence when ``|H y - theta y| <= tol * max(1, |theta|)`` for each wanted pair.
    max_restarts
        Give up (``converged=False``) after this many restarts.
    seed
        Seed of the random start block and of replacement directions.
    start
        Optional ``(b, dim)`` start block.

    Returns
    -------
    LanczosResult
        Ritz values ascending, Ritz vectors as rows, residual norms and counters.
    """
    b = block
    if ncv is None:
        ncv = max(nev + 2 * b, min(max(6 * nev, nev + 6 * b), dim - b - 1))
    if nev + 2 * b > ncv:
        msg = f"ncv={ncv} too small for nev={nev}, block={b}"
        raise ValueError(msg)
    if dim <= ncv + b:
        msg = f"dimension {dim} too small for Lanczos with ncv={ncv}; use a dense solver"
        raise ValueError(msg)

    V = backend.zeros(ncv, dim)
    T = np.zeros((ncv, ncv))
    x0 = backend.random(b, dim, seed) if start is None else start
    q, _ = _orthonormalize(x0, V, 0, backend, seed + 1)
    V[:b] = q
    k, cur, n_matvec = b, 0, 0
    theta = np.zeros(0)
    resid = np.zeros(0)

    for restart in range(max_restarts + 1):
        while True:
            W = matvec(V[cur : cur + b])
            n_matvec += b
            C = backend.inner(V[:k], W)
            T[:k, cur : cur + b] = C
            T[cur : cur + b, :k] = C.T
            backend.subtract_combine(W, C.T, V[:k])
            backend.subtract_combine(W, backend.inner(V[:k], W).T, V[:k])
            q, R = _orthonormalize(W, V, k, backend, seed + 2 + n_matvec)
            if k + b > ncv:
                break
            V[k : k + b] = q
            T[k : k + b, cur : cur + b] = R.T
            T[cur : cur + b, k : k + b] = R
            cur = k
            k += b

        theta, S = sla.eigh(T[:k, :k])
        coupling = R.T @ S[cur : cur + b, :]
        resid = np.linalg.norm(coupling, axis=0)
        wanted = resid[:nev] <= tol * np.maximum(1.0, np.abs(theta[:nev]))
        if bool(np.all(wanted)) or restart == max_restarts:
            vecs = backend.combine(np.ascontiguousarray(S[:, :nev].T), V[:k])
            return LanczosResult(
                eigenvalues=theta[:nev].copy(),
                eigenvectors=vecs,
                residuals=resid[:nev].copy(),
                converged=bool(np.all(wanted)),
                n_matvec=n_matvec,
                n_restart=restart,
            )

        p = min(ncv - 2 * b, max(nev + b, (k + nev) // 2))
        _thick_restart(V, T, S[:, :p], theta[:p], coupling[:, :p], q, k, backend)
        k, cur = p + b, p

    msg = "unreachable"
    raise AssertionError(msg)
