"""Numba bit-manipulation primitives on uint64 Slater determinants.

All arithmetic stays in ``np.uint64``: mixing signed and unsigned integers in Numba promotes to
float64, which would silently corrupt determinants.
"""

from __future__ import annotations

import numpy as np
from numba import njit
from numpy.typing import NDArray

ONE = np.uint64(1)
M1 = np.uint64(0x5555555555555555)
M2 = np.uint64(0x3333333333333333)
M4 = np.uint64(0x0F0F0F0F0F0F0F0F)
H01 = np.uint64(0x0101010101010101)


@njit(cache=True, inline="always")
def popcount(x: np.uint64) -> np.int64:
    """Count the set bits of a uint64."""
    x = x - ((x >> ONE) & M1)
    x = (x & M2) + ((x >> np.uint64(2)) & M2)
    x = (x + (x >> np.uint64(4))) & M4
    return np.int64((x * H01) >> np.uint64(56))


@njit(cache=True, inline="always")
def bit(k: int) -> np.uint64:
    """Return the uint64 with only bit ``k`` set."""
    return ONE << np.uint64(k)


@njit(cache=True, inline="always")
def annihilate(s: np.uint64, k: int) -> tuple[np.uint64, int]:
    """Apply a_k to determinant ``s``; return (new determinant, sign), sign 0 if it vanishes."""
    b = ONE << np.uint64(k)
    if (s & b) == np.uint64(0):
        return s, 0
    sign = -1 if popcount(s & (b - ONE)) & 1 else 1
    return s ^ b, sign


@njit(cache=True, inline="always")
def create(s: np.uint64, k: int) -> tuple[np.uint64, int]:
    """Apply a+_k to determinant ``s``; return (new determinant, sign), sign 0 if it vanishes."""
    b = ONE << np.uint64(k)
    if (s & b) != np.uint64(0):
        return s, 0
    sign = -1 if popcount(s & (b - ONE)) & 1 else 1
    return s | b, sign


@njit(cache=True)
def enumerate_determinants(n_bits: int, n_particles: int, count: int) -> NDArray[np.uint64]:
    """All ``count`` = C(n_bits, n_particles) bit patterns with ``n_particles`` bits set, ascending.

    Uses Gosper's hack.
    """
    out = np.empty(count, dtype=np.uint64)
    if n_particles == 0:
        out[0] = np.uint64(0)
        return out
    v = (ONE << np.uint64(n_particles)) - ONE
    for i in range(count):
        out[i] = v
        if i + 1 < count:
            c = v & (~v + ONE)
            r = v + c
            v = (((r ^ v) >> np.uint64(2)) // c) | r
    return out


@njit(cache=True)
def occupation_stats(
    dets: NDArray[np.uint64],
    morb_orbit: NDArray[np.int32],
    morb_m2: NDArray[np.int32],
    n_orbits: int,
) -> tuple[NDArray[np.int32], NDArray[np.int16]]:
    """Per determinant: total 2M and occupation of each orbit."""
    n = dets.shape[0]
    nm = morb_orbit.shape[0]
    m2 = np.zeros(n, dtype=np.int32)
    occ = np.zeros((n, n_orbits), dtype=np.int16)
    for i in range(n):
        s = dets[i]
        for k in range(nm):
            if (s >> np.uint64(k)) & ONE:
                m2[i] += morb_m2[k]
                occ[i, morb_orbit[k]] += 1
    return m2, occ


@njit(cache=True, inline="always")
def find_det(
    s: np.uint64,
    sorted_dets: NDArray[np.uint64],
    det_sub: NDArray[np.int32],
    det_local: NDArray[np.int32],
) -> tuple[int, int]:
    """Locate determinant ``s``; return (sub-block id, local index) or (-1, -1)."""
    i = np.searchsorted(sorted_dets, s)
    if i < sorted_dets.shape[0] and sorted_dets[i] == s:
        return det_sub[i], det_local[i]
    return -1, -1
