"""Lanczos strength functions: ``S(E)`` of an operator acting on a state, from a pivot vector.

For a pivot ``p = O psi_0``, N plain Lanczos steps give a tridiagonal ``T`` with eigenvalues
``theta_k`` and eigenvectors ``y_k``. The weights ``|p|^2 y_(1k)^2`` reproduce the first 2N
moments ``<p|H^m|p>`` of the strength distribution exactly in exact arithmetic, and ``T`` gives
the resolvent ``<p|(z - H)^-1|p>`` as a continued fraction. Only three vectors are needed
without reorthogonalization; with it, all Lanczos vectors are kept and Ritz values do not repeat.

The pivot has good J only if ``psi_0`` has J = 0 or the pivot is projected; otherwise ``S(E)``
sums over final J and must be labelled as such.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla
from numpy.typing import NDArray

from cosmo_gpu.backend import Array, VectorBackend

BREAKDOWN = 1e-12


@dataclass(frozen=True)
class Strength:
    """Lanczos representation of a strength distribution."""

    alpha: NDArray[np.float64]
    beta: NDArray[np.float64]  # beta[k] couples step k and k + 1; len(beta) == len(alpha) - 1
    norm2: float  # |p|^2, the total strength
    energies: NDArray[np.float64]
    weights: NDArray[np.float64]

    @property
    def n_steps(self) -> int:
        """Number of Lanczos steps taken (fewer than asked at an invariant subspace)."""
        return int(self.alpha.shape[0])

    def moment(self, k: int) -> float:
        """``sum_k w_k theta_k^m`` (equals ``<p|H^m|p>`` for ``m <= 2 N - 1``)."""
        return float(np.sum(self.weights * self.energies**k))

    def continued_fraction(self, z: complex) -> complex:
        """``<p|(z - H)^-1|p>`` from the tridiagonal coefficients."""
        g: complex = 0.0
        for k in range(self.n_steps - 1, -1, -1):
            b2 = self.beta[k] ** 2 if k < self.n_steps - 1 else 0.0
            g = 1.0 / (z - self.alpha[k] - b2 * g)
        return self.norm2 * g

    def smoothed(
        self, grid: NDArray[np.float64], width: float, shape: str = "lorentzian"
    ) -> NDArray[np.float64]:
        """``S(E)`` broadened by a Lorentzian (half width) or Gaussian (sigma) of ``width``."""
        d = grid[:, None] - self.energies[None, :]
        if shape == "lorentzian":
            kern = width / np.pi / (d**2 + width**2)
        elif shape == "gaussian":
            kern = np.exp(-0.5 * (d / width) ** 2) / (width * np.sqrt(2 * np.pi))
        else:
            msg = f"unknown shape {shape!r}"
            raise ValueError(msg)
        return np.asarray(kern @ self.weights, dtype=np.float64)


def lanczos_strength(
    matvec: Callable[[Array], Array],
    pivot: Array,
    n_steps: int,
    backend: VectorBackend,
    reorthogonalize: bool = False,
) -> Strength:
    """Run up to ``n_steps`` plain Lanczos steps from ``pivot`` (shape ``(1, D)``).

    Stops early when the Krylov space becomes invariant (``beta < BREAKDOWN * |alpha|``).
    """
    norm2 = float(backend.inner(pivot, pivot)[0, 0])
    if norm2 <= 0.0:
        return Strength(np.zeros(0), np.zeros(0), 0.0, np.zeros(0), np.zeros(0))
    v = backend.combine(np.array([[1.0 / np.sqrt(norm2)]]), pivot)
    dim = v.shape[1]
    keep = backend.zeros(min(n_steps, dim), dim) if reorthogonalize else None
    v_prev = None
    alpha: list[float] = []
    beta: list[float] = []
    for k in range(min(n_steps, dim)):
        if keep is not None:
            keep[k : k + 1] = v
        w = matvec(v)
        a = float(backend.inner(v, w)[0, 0])
        alpha.append(a)
        backend.subtract_combine(w, np.array([[a]]), v)
        if v_prev is not None:
            backend.subtract_combine(w, np.array([[beta[-1]]]), v_prev)
        if keep is not None:
            for _ in range(2):
                c = backend.inner(w, keep[: k + 1])
                backend.subtract_combine(w, c, keep[: k + 1])
        b = float(np.sqrt(max(backend.inner(w, w)[0, 0], 0.0)))
        if k == min(n_steps, dim) - 1 or b <= BREAKDOWN * max(1.0, abs(a)):
            break
        beta.append(b)
        v_prev = v
        v = backend.combine(np.array([[1.0 / b]]), w)
    al = np.array(alpha)
    be = np.array(beta)
    theta, y = sla.eigh_tridiagonal(al, be) if len(al) > 1 else (al.copy(), np.ones((1, 1)))
    return Strength(al, be, norm2, theta, norm2 * y[0, :] ** 2)
