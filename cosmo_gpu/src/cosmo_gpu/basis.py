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

from dataclasses import dataclass
from functools import cache
from itertools import combinations
from math import comb

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.bits_kernels import enumerate_determinants, occupation_stats
from cosmo_gpu.modelspace import ModelSpace
from cosmo_gpu.truncation import Truncation, orbit_quanta


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


def _assemble(
    space: ModelSpace,
    n_particles: int,
    dets: NDArray[np.uint64],
    m2: NDArray[np.int32],
    part_id: NDArray[np.int64],
    partitions: NDArray[np.int16],
) -> SpeciesBasis:
    """Sort determinants into (2M, partition) sub-blocks and build the lookup arrays."""
    order = np.lexsort((dets, part_id, m2))
    dets_o = dets[order]
    keys = np.stack([m2[order].astype(np.int64), part_id[order]], axis=1)
    change = np.ones(len(order), dtype=bool)
    change[1:] = np.any(keys[1:] != keys[:-1], axis=1)
    starts = np.flatnonzero(change)
    sub_start = np.append(starts, len(order)).astype(np.int64)
    sub_m2 = keys[starts, 0].astype(np.int32)
    sub_partition = partitions[keys[starts, 1]].astype(np.int16).reshape(-1, space.n_orbits)
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


def build_species_basis_full(space: ModelSpace, n_particles: int) -> SpeciesBasis:
    """Enumerate all C(n_morbits, n) determinants (small spaces and tests only)."""
    nm = space.n_morbits
    if not 0 <= n_particles <= nm:
        msg = f"{n_particles} particles do not fit in {nm} m-orbits"
        raise ValueError(msg)
    dets = enumerate_determinants(nm, n_particles, comb(nm, n_particles))
    m2, occ = occupation_stats(dets, space.morb_orbit, space.morb_m2, space.n_orbits)
    partitions, part_id = np.unique(occ, axis=0, return_inverse=True)
    return _assemble(space, n_particles, dets, m2, part_id.reshape(-1), partitions)


def species_partitions(
    space: ModelSpace, n_particles: int, max_quanta: int | None = None
) -> list[tuple[int, ...]]:
    """All occupations of the orbits with ``n_particles`` in total and quanta <= ``max_quanta``."""
    q = [int(x) for x in orbit_quanta(space)]
    deg = [o.degeneracy for o in space.orbits]
    norb = space.n_orbits
    # room[i]: particles that orbits i.. can still hold; qmin_rest[i][k]: least quanta for k there
    room = [sum(deg[i:]) for i in range(norb + 1)]
    out: list[tuple[int, ...]] = []
    cap = max_quanta if max_quanta is not None else 10**9

    def least_quanta(i: int, k: int) -> int:
        levels = sorted(qq for j in range(i, norb) for qq in [q[j]] * deg[j])
        return sum(levels[:k])

    least = [[least_quanta(i, k) for k in range(n_particles + 1)] for i in range(norb + 1)]

    def rec(i: int, left: int, used_q: int, occ: list[int]) -> None:
        if i == norb:
            if left == 0:
                out.append(tuple(occ))
            return
        for k in range(min(deg[i], left), -1, -1):
            rest = left - k
            if rest > room[i + 1]:
                break
            qk = used_q + k * q[i]
            if qk + least[i + 1][rest] > cap:
                continue
            occ.append(k)
            rec(i + 1, rest, qk, occ)
            occ.pop()

    rec(0, n_particles, 0, [])
    out.sort()  # lexicographic, the same order as np.unique in the full enumeration
    return out


@cache
def _orbit_patterns(
    offset: int, degeneracy: int, j2: int, k: int
) -> tuple[NDArray[np.uint64], NDArray[np.int32]]:
    """Bit patterns and 2M of all ways to put ``k`` particles in one orbit's m-substates."""
    m2 = list(range(-j2, j2 + 1, 2))
    bits, ms = [], []
    for c in combinations(range(degeneracy), k):
        bits.append(sum(1 << (offset + i) for i in c))
        ms.append(sum(m2[i] for i in c))
    return np.array(bits, dtype=np.uint64), np.array(ms, dtype=np.int32)


def build_species_basis(
    space: ModelSpace, n_particles: int, max_quanta: int | None = None
) -> SpeciesBasis:
    """Build one species' basis partition by partition, never enumerating the full space.

    Only partitions with total quanta <= ``max_quanta`` are generated, so an Nħω-truncated
    psdpf basis costs time proportional to what is kept.
    """
    nm = space.n_morbits
    if not 0 <= n_particles <= nm:
        msg = f"{n_particles} particles do not fit in {nm} m-orbits"
        raise ValueError(msg)
    parts = species_partitions(space, n_particles, max_quanta)
    if not parts:
        msg = f"no partition of {n_particles} particles within {max_quanta} quanta"
        raise ValueError(msg)
    offsets = np.concatenate([[0], np.cumsum([o.degeneracy for o in space.orbits])])
    det_chunks, m2_chunks, id_chunks = [], [], []
    for pid, occ in enumerate(parts):
        bits = np.zeros(1, dtype=np.uint64)
        m2 = np.zeros(1, dtype=np.int32)
        for i, k in enumerate(occ):
            if k == 0:
                continue
            o = space.orbits[i]
            b, m = _orbit_patterns(int(offsets[i]), o.degeneracy, o.j2, k)
            bits = (bits[:, None] | b[None, :]).ravel()
            m2 = (m2[:, None] + m[None, :]).ravel()
        det_chunks.append(bits)
        m2_chunks.append(m2)
        id_chunks.append(np.full(bits.shape[0], pid, dtype=np.int64))
    partitions = np.array(parts, dtype=np.int16).reshape(-1, space.n_orbits)
    return _assemble(
        space,
        n_particles,
        np.concatenate(det_chunks),
        np.concatenate(m2_chunks),
        np.concatenate(id_chunks),
        partitions,
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
    truncation: Truncation | None = None

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
    """Pair sub-blocks into blocks with total 2M ``m2``, the requested parity and truncation."""
    if parity is None:
        options = possible_parities(protons, neutrons)
        if len(options) != 1:
            msg = "this model space mixes parities: pass parity=+1 or -1"
            raise ValueError(msg)
        parity = options.pop()
    keep_p = _truncation_features(protons, truncation, proton=True)
    keep_n = _truncation_features(neutrons, truncation, proton=False)
    psub_chunks: list[NDArray[np.int64]] = []
    nsub_chunks: list[NDArray[np.int64]] = []
    by_m2: dict[int, NDArray[np.int64]] = {}
    for mval in np.unique(neutrons.sub_m2):
        by_m2[int(mval)] = np.flatnonzero(neutrons.sub_m2 == mval)
    for i in range(protons.n_sub):
        cand = by_m2.get(m2 - int(protons.sub_m2[i]))
        if cand is None:
            continue
        ok = neutrons.sub_parity[cand] * int(protons.sub_parity[i]) == parity
        if truncation is not None:
            ok &= _pair_mask(truncation, keep_p[i], keep_n[cand])
        js = cand[ok]
        if js.size:
            psub_chunks.append(np.full(js.size, i, dtype=np.int64))
            nsub_chunks.append(js)
    psub = np.concatenate(psub_chunks) if psub_chunks else np.zeros(0, dtype=np.int64)
    nsub = np.concatenate(nsub_chunks) if nsub_chunks else np.zeros(0, dtype=np.int64)
    sizes = protons.sub_dim[psub] * neutrons.sub_dim[nsub]
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
        truncation=truncation,
    )


def _truncation_features(
    sb: SpeciesBasis, trunc: Truncation | None, proton: bool
) -> NDArray[np.int64]:
    """Per sub-block: total quanta, then the occupation counted by each rule (columns)."""
    if trunc is None:
        return np.zeros((sb.n_sub, 0), dtype=np.int64)
    part = sb.sub_partition.astype(np.int64)
    cols = [part @ orbit_quanta(trunc.space)]
    for rule in trunc.occupation:
        mask = np.array(rule.proton_orbits if proton else rule.neutron_orbits, dtype=np.int64)
        cols.append(part @ mask)
    return np.stack(cols, axis=1)


def _pair_mask(
    trunc: Truncation, fp: NDArray[np.int64], fn: NDArray[np.int64]
) -> NDArray[np.bool_]:
    """Vectorized ``Truncation.keeps`` for one proton sub-block against many neutron ones."""
    mask = np.ones(fn.shape[0], dtype=bool)
    if trunc.hw is not None:
        excess = fp[0] + fn[:, 0] - trunc.min_quanta_p - trunc.min_quanta_n
        mask &= (excess >= trunc.hw[0]) & (excess <= trunc.hw[1])
    for r, rule in enumerate(trunc.occupation):
        occ = fp[1 + r] + fn[:, 1 + r]
        mask &= (occ >= rule.nmin) & (occ <= rule.nmax)
    return mask


def build_basis(
    space: ModelSpace,
    z_val: int,
    n_val: int,
    m2: int | None = None,
    parity: int | None = None,
    truncation: Truncation | None = None,
) -> MBasis:
    """Build both species (pruned by the truncation's quanta cap) and the block list.

    ``m2`` defaults to ``(z_val + n_val) % 2`` (M = 0 or 1/2), which contains every J.
    """
    cap_p = truncation.species_quanta_cap(proton=True) if truncation else None
    cap_n = truncation.species_quanta_cap(proton=False) if truncation else None
    protons = build_species_basis(space, z_val, cap_p)
    neutrons = build_species_basis(space, n_val, cap_n)
    m2 = (z_val + n_val) % 2 if m2 is None else m2
    return build_mbasis(space, protons, neutrons, m2, parity, truncation)
