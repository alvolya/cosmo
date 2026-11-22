"""Vector backends: the only place the solver touches an array library.

A block of vectors is a C-contiguous 2D array of shape ``(nvec, D)``. Small dense matrices
(projections, Ritz coefficients) always live on the host as NumPy arrays. A CuPy backend with
the same methods comes with the GPU kernel in stage 4.
"""

from __future__ import annotations

from typing import Any, Protocol

import numpy as np
from numpy.typing import NDArray

#: A backend array (``numpy.ndarray`` here, ``cupy.ndarray`` on GPU).
Array = Any


class VectorBackend(Protocol):
    """Operations the solver needs on blocks of vectors."""

    def zeros(self, nvec: int, dim: int) -> Array:
        """Allocate a zero block."""
        ...

    def random(self, nvec: int, dim: int, seed: int) -> Array:
        """Allocate a block of standard-normal entries from a seeded generator."""
        ...

    def inner(self, a: Array, b: Array) -> NDArray[np.float64]:
        """Return ``a @ b.T`` on the host."""
        ...

    def combine(self, coeffs: NDArray[np.float64], a: Array) -> Array:
        """Return the new block ``coeffs @ a``."""
        ...

    def subtract_combine(self, w: Array, coeffs: NDArray[np.float64], a: Array) -> None:
        """In place ``w -= coeffs @ a``."""
        ...

    def to_host(self, a: Array) -> NDArray[np.float64]:
        """Copy to a host NumPy array."""
        ...

    def from_host(self, a: NDArray[np.float64]) -> Array:
        """Copy a host array to the backend."""
        ...


class NumpyBackend:
    """Host backend on NumPy."""

    def zeros(self, nvec: int, dim: int) -> NDArray[np.float64]:
        """Allocate a zero block."""
        return np.zeros((nvec, dim), dtype=np.float64)

    def random(self, nvec: int, dim: int, seed: int) -> NDArray[np.float64]:
        """Allocate a block of standard-normal entries from a seeded generator."""
        return np.random.default_rng(seed).standard_normal((nvec, dim))

    def inner(self, a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return ``a @ b.T``."""
        return np.asarray(a @ b.T, dtype=np.float64)

    def combine(self, coeffs: NDArray[np.float64], a: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return the new block ``coeffs @ a``."""
        return np.ascontiguousarray(coeffs @ a)

    def subtract_combine(
        self, w: NDArray[np.float64], coeffs: NDArray[np.float64], a: NDArray[np.float64]
    ) -> None:
        """In place ``w -= coeffs @ a``."""
        w -= coeffs @ a

    def to_host(self, a: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return ``a`` itself (already on the host)."""
        return a

    def from_host(self, a: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return a contiguous float64 copy of ``a``."""
        return np.ascontiguousarray(a, dtype=np.float64)
