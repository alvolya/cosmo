"""M-scheme basis: per-species sub-blocks and proton x neutron blocks.

A *sub-block* of one species is the set of its Slater determinants with a given partition
(occupation of each spherical orbit) and total 2M. Its determinants are stored contiguously and
sorted, so a determinant is addressed as (sub-block id, local index).

A *block* pairs a proton sub-block with a neutron sub-block whose 2M values add to the target and
whose parity product equals the target parity. Its amplitudes form a dense row-major
``(d_p, d_n)`` matrix, and the global index of amplitude ``(ip, in)`` of block ``B`` is
``block_offset[B] + ip * d_n + in``. Truncation removes whole blocks only.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import comb

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.bits_kernels import enumerate_determinants, occupation_stats
from cosmo_gpu.modelspace import ModelSpace

#: Predicate on (proton partition, neutron partition) deciding whether a block is kept.
Truncation = Callable[[NDArray[np.int16], NDArray[np.int16]], bool]


@dataclass(frozen=True)
class SpeciesBasis:
    """All determinants of one species, grouped into (partition, 2M) sub-blocks."""

    n_particles: int
    dets: NDArray[np.uint64]
    sub_start: NDArray[np.int64]
    sub_m2: NDArray[np.int32]
    sub_parity: NDArray[np.int8]
    sub_partition: NDArray[np.int16]
    sorted_dets: NDArray[np.uint64]
    det_sub: NDArray[np.int32]
    det_local: NDArray[np.int32]

    @property
    def n_sub(self) -> int:
        """Number of sub-blocks."""
        return int(self.sub_m2.shape[0])

    @property
    def sub_dim(self) -> NDArray[np.int64]:
        """Dimension of each sub-block."""
        return np.diff(self.sub_start)


def build_species_basis(space: ModelSpace, n_particles: int) -> SpeciesBasis:
    """Enumerate every determinant of ``n_particles`` in ``space`` and sort into sub-blocks."""
    nm = space.n_morbits
    if not 0 <= n_particles <= nm:
        msg = f"{n_particles} particles do not fit in {nm} m-orbits"
        raise ValueError(msg)
    dets = enumerate_determinants(nm, n_particles, comb(nm, n_particles))
    m2, occ = occupation_stats(dets, space.morb_orbit, space.morb_m2, space.n_orbits)
    partitions, part_id = np.unique(occ, axis=0, return_inverse=True)
    part_id = part_id.reshape(-1)
    # sub-block key: (2M, partition id); determinants sorted inside each sub-block
    order = np.lexsort((dets, part_id, m2))
    dets_o = dets[order]
    keys = np.stack([m2[order], part_id[order]], axis=1)
    change = np.ones(len(order), dtype=bool)
    change[1:] = np.any(keys[1:] != keys[:-1], axis=1)
    starts = np.flatnonzero(change)
    sub_start = np.append(starts, len(order)).astype(np.int64)
    sub_m2 = keys[starts, 0].astype(np.int32)
    sub_partition = partitions[keys[starts, 1]].astype(np.int16)
    odd = sub_partition[:, space.orbit_parity < 0].sum(axis=1) % 2
    sub_parity = np.where(odd == 1, -1, 1).astype(np.int8)
    det_sub = np.repeat(np.arange(len(starts), dtype=np.int32), np.diff(sub_start))
    det_local = (np.arange(len(order)) - sub_start[det_sub]).astype(np.int32)
    lookup = np.argsort(dets_o)
    return SpeciesBasis(
        n_particles=n_particles,
        dets=dets_o,
        sub_start=sub_start,
        sub_m2=sub_m2,
        sub_parity=sub_parity,
        sub_partition=sub_partition,
        sorted_dets=dets_o[lookup],
        det_sub=det_sub[lookup],
        det_local=det_local[lookup],
    )


@dataclass(frozen=True)
class MBasis:
    """The many-body M-scheme basis: a list of (proton sub-block, neutron sub-block) blocks."""

    space: ModelSpace
    protons: SpeciesBasis
    neutrons: SpeciesBasis
    m2: int
    parity: int
    block_psub: NDArray[np.int32]
    block_nsub: NDArray[np.int32]
    block_offset: NDArray[np.int64]
    block_of: NDArray[np.int32]

    @property
    def dim(self) -> int:
        """Total dimension D."""
        return int(self.block_offset[-1])

    @property
    def n_blocks(self) -> int:
        """Number of blocks."""
        return int(self.block_psub.shape[0])


def possible_parities(protons: SpeciesBasis, neutrons: SpeciesBasis) -> set[int]:
    """Parities reachable by some proton x neutron sub-block pair."""
    return {int(a) * int(b) for a in set(protons.sub_parity) for b in set(neutrons.sub_parity)}


def build_mbasis(
    space: ModelSpace,
    protons: SpeciesBasis,
    neutrons: SpeciesBasis,
    m2: int,
    parity: int | None = None,
    truncation: Truncation | None = None,
) -> MBasis:
    """Pair sub-blocks into blocks with total 2M ``m2`` and the requested parity."""
    if parity is None:
        options = possible_parities(protons, neutrons)
        if len(options) != 1:
            msg = "this model space mixes parities: pass parity=+1 or -1"
            raise ValueError(msg)
        parity = options.pop()
    psub: list[int] = []
    nsub: list[int] = []
    sizes: list[int] = []
    pdim, ndim = protons.sub_dim, neutrons.sub_dim
    by_m2: dict[int, list[int]] = {}
    for j in range(neutrons.n_sub):
        by_m2.setdefault(int(neutrons.sub_m2[j]), []).append(j)
    for i in range(protons.n_sub):
        for j in by_m2.get(m2 - int(protons.sub_m2[i]), []):
            if int(protons.sub_parity[i]) * int(neutrons.sub_parity[j]) != parity:
                continue
            if truncation is not None and not truncation(
                protons.sub_partition[i], neutrons.sub_partition[j]
            ):
                continue
            psub.append(i)
            nsub.append(j)
            sizes.append(int(pdim[i] * ndim[j]))
    block_of = np.full((protons.n_sub, neutrons.n_sub), -1, dtype=np.int32)
    block_of[psub, nsub] = np.arange(len(psub), dtype=np.int32)
    return MBasis(
        space=space,
        protons=protons,
        neutrons=neutrons,
        m2=m2,
        parity=parity,
        block_psub=np.array(psub, dtype=np.int32),
        block_nsub=np.array(nsub, dtype=np.int32),
        block_offset=np.concatenate([[0], np.cumsum(sizes, dtype=np.int64)]).astype(np.int64),
        block_of=block_of,
    )
