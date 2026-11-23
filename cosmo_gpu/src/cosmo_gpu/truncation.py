"""Basis truncations: oscillator quanta (Nħω) and orbit-occupation windows.

Semantics follow CoSMo's ``XsysmbsOMP``:

* Each orbit carries ``2n + l`` oscillator quanta. ``min_quanta(space, n)`` is the smallest total
  for ``n`` particles of one species. A block (proton partition, neutron partition) is kept when
  ``Q_p + Q_n - min_p - min_n`` lies in ``[hw_min, hw_max]``.
* An occupation rule sums the occupations of chosen orbits (protons, neutrons or both) and keeps
  the block when the sum lies in ``[nmin, nmax]``.

CoSMo's ``CoreQuanta`` (``include/sm/SingleParticleStates.cxx``) sorts an empty range, so it sums
the first ``n`` m-orbits in file order instead of the lowest ones. The two agree whenever the
valence space lists orbits in ascending quanta, as all spaces used here do; ``min_quanta`` is the
true minimum.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.modelspace import ModelSpace


@dataclass(frozen=True)
class OccupationRule:
    """Keep blocks whose occupation of the masked orbits lies in ``[nmin, nmax]``."""

    proton_orbits: tuple[bool, ...]
    neutron_orbits: tuple[bool, ...]
    nmin: int
    nmax: int


@dataclass(frozen=True)
class Truncation:
    """A conjunction of an optional Nħω window and occupation rules."""

    space: ModelSpace
    hw: tuple[int, int] | None = None
    occupation: tuple[OccupationRule, ...] = ()
    min_quanta_p: int = 0
    min_quanta_n: int = 0

    def keeps(self, part_p: NDArray[np.int16], part_n: NDArray[np.int16]) -> bool:
        """Return True if the block (proton partition, neutron partition) is in the space."""
        if self.hw is not None:
            q = orbit_quanta(self.space)
            excess = int(part_p @ q) + int(part_n @ q) - self.min_quanta_p - self.min_quanta_n
            if not self.hw[0] <= excess <= self.hw[1]:
                return False
        for rule in self.occupation:
            occ = int(part_p[list(rule.proton_orbits)].sum()) + int(
                part_n[list(rule.neutron_orbits)].sum()
            )
            if not rule.nmin <= occ <= rule.nmax:
                return False
        return True

    def species_quanta_cap(self, proton: bool) -> int | None:
        """Largest quanta one species can carry in any kept block, or None if unbounded."""
        if self.hw is None:
            return None
        return (self.min_quanta_p if proton else self.min_quanta_n) + self.hw[1]

    def describe(self) -> str:
        """Short human-readable form, e.g. ``hw=1:1 occ[...]``."""
        parts = []
        if self.hw is not None:
            parts.append(f"hw={self.hw[0]}:{self.hw[1]}")
        for r in self.occupation:
            names = [
                o.name + t
                for o, p, n in zip(
                    self.space.orbits, r.proton_orbits, r.neutron_orbits, strict=True
                )
                for t, on in (("p", p), ("n", n))
                if on
            ]
            parts.append(f"occ({','.join(names)})={r.nmin}:{r.nmax}")
        return " ".join(parts) if parts else "none"


def orbit_quanta(space: ModelSpace) -> NDArray[np.int64]:
    """Oscillator quanta ``2n + l`` of each orbit."""
    return np.array([2 * o.n + o.l for o in space.orbits], dtype=np.int64)


def min_quanta(space: ModelSpace, n_particles: int) -> int:
    """Smallest total quanta of ``n_particles`` of one species (lowest m-orbits filled)."""
    q = np.sort(orbit_quanta(space)[space.morb_orbit])
    if n_particles > q.shape[0]:
        msg = f"{n_particles} particles do not fit in {q.shape[0]} m-orbits"
        raise ValueError(msg)
    return int(q[:n_particles].sum())


def make_truncation(
    space: ModelSpace,
    z_val: int,
    n_val: int,
    hw: tuple[int, int] | None = None,
    occupation: tuple[OccupationRule, ...] = (),
) -> Truncation:
    """Build a truncation for ``z_val`` protons and ``n_val`` neutrons."""
    return Truncation(
        space=space,
        hw=hw,
        occupation=occupation,
        min_quanta_p=min_quanta(space, z_val),
        min_quanta_n=min_quanta(space, n_val),
    )


def parse_hw(text: str) -> tuple[int, int]:
    """Parse ``MIN:MAX`` (or a single ``N`` meaning ``0:N``)."""
    lo, sep, hi = text.partition(":")
    if not sep:
        return (0, int(lo))
    if not lo.strip() or not hi.strip():
        msg = f"bad Nhw window {text!r}: give MIN:MAX or a single N"
        raise ValueError(msg)
    pair = (int(lo), int(hi))
    if not 0 <= pair[0] <= pair[1]:
        msg = f"bad Nhw window {text!r}: need 0 <= MIN <= MAX"
        raise ValueError(msg)
    return pair


def parse_rejection(text: str, space: ModelSpace) -> OccupationRule:
    """Parse ``ORBS:TYPE:NMIN:NMAX``, e.g. ``0f5,1p3,1p1:pn:0:4`` (as CoSMo's make_system)."""
    try:
        orbs, kind, nmin, nmax = text.split(":")
    except ValueError:
        msg = f"bad rejection {text!r}, expected ORBS:TYPE:NMIN:NMAX"
        raise ValueError(msg) from None
    if kind not in ("p", "n", "pn"):
        msg = f"bad rejection type {kind!r}, expected p, n or pn"
        raise ValueError(msg)
    names = {o.name: i for i, o in enumerate(space.orbits)}
    mask = [False] * space.n_orbits
    for name in orbs.split(","):
        if name not in names:
            msg = f"orbit {name!r} not in {space.name}: {sorted(names)}"
            raise ValueError(msg)
        mask[names[name]] = True
    none = (False,) * space.n_orbits
    return OccupationRule(
        proton_orbits=tuple(mask) if "p" in kind else none,
        neutron_orbits=tuple(mask) if "n" in kind else none,
        nmin=int(nmin),
        nmax=int(nmax),
    )
