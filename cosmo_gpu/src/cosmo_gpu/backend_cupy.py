"""CuPy vector backend: the same operations as ``NumpyBackend`` on device arrays.

Imported only when a GPU is used, so the package works without CuPy installed. Small dense
matrices (coefficients, projections) stay on the host as NumPy arrays and are copied to the
device per call; they are at most ``ncv x ncv``.
"""

from __future__ import annotations

from typing import Any

import cupy as cp
import numpy as np
from numpy.typing import NDArray

#: A device array (``cupy.ndarray``).
DeviceArray = Any


class CupyBackend:
    """Device backend on CuPy, FP64."""

    def zeros(self, nvec: int, dim: int) -> DeviceArray:
        """Allocate a zero block on the device."""
        return cp.zeros((nvec, dim), dtype=cp.float64)

    def random(self, nvec: int, dim: int, seed: int) -> DeviceArray:
        """Standard-normal block from the same seeded host generator as ``NumpyBackend``."""
        return cp.asarray(np.random.default_rng(seed).standard_normal((nvec, dim)))

    def inner(self, a: DeviceArray, b: DeviceArray) -> NDArray[np.float64]:
        """Return ``a @ b.T`` on the host."""
        return np.asarray(cp.asnumpy(a @ b.T), dtype=np.float64)

    def combine(self, coeffs: NDArray[np.float64], a: DeviceArray) -> DeviceArray:
        """Return the new block ``coeffs @ a``."""
        return cp.ascontiguousarray(cp.asarray(coeffs) @ a)

    def subtract_combine(self, w: DeviceArray, coeffs: NDArray[np.float64], a: DeviceArray) -> None:
        """In place ``w -= coeffs @ a``."""
        w -= cp.asarray(coeffs) @ a

    def to_host(self, a: DeviceArray) -> NDArray[np.float64]:
        """Copy to a host NumPy array."""
        return np.asarray(cp.asnumpy(a), dtype=np.float64)

    def from_host(self, a: NDArray[np.float64]) -> DeviceArray:
        """Copy a host array to the device."""
        return cp.ascontiguousarray(cp.asarray(a, dtype=cp.float64))
