"""Model spaces: spherical orbits, m-orbits and the CoSMo interaction catalogue.

Angular momenta are stored doubled (``j2 = 2j``, ``m2 = 2m``) so that every quantity is an
integer. Orbits keep CoSMo's 1-based ``index``, which is the label used in OXBASH ``.int``
files. Both species share the same orbit list (``type="pn"`` valence spaces).

M-orbit ``k`` of a species is bit ``k`` of a determinant. M-orbits are ordered by orbit (in
file order) and, within an orbit, by ascending ``m``.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cached_property
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

MAX_MORBITS = 64

# fmt: off
ELEMENTS = (  # noqa: SIM905 (index = Z; "n" is the neutron)
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As "
    "Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd "
    "Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U"
).split()
# fmt: on


@dataclass(frozen=True)
class Orbit:
    """A spherical single-particle orbit (n, l, j)."""

    name: str
    n: int
    l: int
    j2: int
    index: int

    @property
    def parity(self) -> int:
        """Parity of the orbit, +1 or -1."""
        return 1 if self.l % 2 == 0 else -1

    @property
    def degeneracy(self) -> int:
        """Number of m-substates, 2j + 1."""
        return self.j2 + 1


@dataclass(frozen=True)
class ModelSpace:
    """A valence space: orbits shared by protons and neutrons, plus the inert core."""

    name: str
    orbits: tuple[Orbit, ...]
    core_z: int
    core_n: int
    interactions_dir: Path | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        if self.n_morbits > MAX_MORBITS:
            msg = f"{self.name}: {self.n_morbits} m-orbits per species exceeds {MAX_MORBITS}"
            raise ValueError(msg)

    @property
    def n_orbits(self) -> int:
        """Number of spherical orbits."""
        return len(self.orbits)

    @property
    def n_morbits(self) -> int:
        """Number of m-orbits per species (bits per determinant)."""
        return sum(o.degeneracy for o in self.orbits)

    @cached_property
    def morb_orbit(self) -> NDArray[np.int32]:
        """Orbit position (0-based, in ``orbits`` order) of each m-orbit."""
        return np.array(
            [i for i, o in enumerate(self.orbits) for _ in range(o.degeneracy)], dtype=np.int32
        )

    @cached_property
    def morb_m2(self) -> NDArray[np.int32]:
        """Doubled projection 2m of each m-orbit."""
        return np.array(
            [m2 for o in self.orbits for m2 in range(-o.j2, o.j2 + 1, 2)], dtype=np.int32
        )

    @cached_property
    def morb_j2(self) -> NDArray[np.int32]:
        """Doubled angular momentum 2j of each m-orbit."""
        return np.array([self.orbits[i].j2 for i in self.morb_orbit], dtype=np.int32)

    @cached_property
    def orbit_parity(self) -> NDArray[np.int8]:
        """Parity (+1 or -1) of each orbit."""
        return np.array([o.parity for o in self.orbits], dtype=np.int8)

    def orbit_position(self, cosmo_index: int) -> int:
        """Map a CoSMo / OXBASH 1-based orbit label to its position in ``orbits``."""
        for i, o in enumerate(self.orbits):
            if o.index == cosmo_index:
                return i
        msg = f"orbit index {cosmo_index} not in model space {self.name}"
        raise KeyError(msg)


@dataclass(frozen=True)
class InteractionFile:
    """One OXBASH file of a Hamiltonian, with CoSMo's optional rank scalings."""

    path: Path
    spe_scale: float = 1.0
    tbme_scale: float = 1.0


@dataclass(frozen=True)
class HamiltonianSpec:
    """A named Hamiltonian from the CoSMo catalogue: its valence space and files."""

    name: str
    valence_space: str
    files: tuple[InteractionFile, ...]


def default_interactions_dir() -> Path:
    """Return the CoSMo ``interactions`` directory.

    Uses ``$COSMO/interactions`` if set, otherwise the ``interactions`` directory of the repository
    this package lives in.
    """
    cosmo = os.environ.get("COSMO")
    if cosmo and (Path(cosmo) / "interactions").is_dir():
        return Path(cosmo) / "interactions"
    return Path(__file__).resolve().parents[3] / "interactions"


def find_interaction_file(name: str, interactions_dir: Path | None = None) -> Path:
    """Locate an interaction file the way CoSMo does.

    Searches ``$SMINTPATH`` (colon-separated) first, then ``interactions_dir`` and each of its
    subdirectories.
    """
    candidates: list[Path] = [Path(p) for p in os.environ.get("SMINTPATH", "").split(":") if p]
    base = interactions_dir or default_interactions_dir()
    candidates.append(base)
    if base.is_dir():
        candidates.extend(sorted(p for p in base.iterdir() if p.is_dir()))
    for d in candidates:
        f = d / name
        if f.is_file():
            return f
    msg = f"interaction file {name!r} not found in SMINTPATH or {base}"
    raise FileNotFoundError(msg)


def _catalogue(interactions_dir: Path | None) -> ET.Element:
    base = interactions_dir or default_interactions_dir()
    return ET.parse(base / "int.xml").getroot()


def spin2(text: str) -> int:
    """Doubled angular momentum of a spin string such as ``7/2`` or ``0``."""
    return _j2(text)


def _j2(text: str) -> int:
    return int(2 * Fraction(text))


def parse_valence_space(vs: ET.Element, interactions_dir: Path | None = None) -> ModelSpace:
    """Build a ``ModelSpace`` from a ``<valencespace>`` element (``int.xml`` or a system XML)."""
    name = vs.get("name", "")
    core = vs.find("core")
    if core is None:
        msg = f"valence space {name} has no core"
        raise ValueError(msg)
    orbits = []
    for o in vs.findall("orbital"):
        if o.get("type", "pn") != "pn":
            msg = f"valence space {name}: only type='pn' orbitals are supported"
            raise NotImplementedError(msg)
        orbits.append(
            Orbit(
                name=o.get("name", ""),
                n=int(o.get("n", "0")),
                l=int(o.get("l", "0")),
                j2=_j2(o.get("j", "0")),
                index=int(o.get("index", "0")),
            )
        )
    return ModelSpace(
        name=name,
        orbits=tuple(orbits),
        core_z=int(core.get("Z", "0")),
        core_n=int(core.get("N", "0")),
        interactions_dir=interactions_dir,
    )


def load_model_space(name: str, interactions_dir: Path | None = None) -> ModelSpace:
    """Read valence space ``name`` from the CoSMo catalogue ``int.xml``."""
    root = _catalogue(interactions_dir)
    for vs in root.iter("valencespace"):
        if vs.get("name") == name:
            return parse_valence_space(vs, interactions_dir)
    msg = f"valence space {name!r} not found in int.xml"
    raise KeyError(msg)


def parse_hamiltonian(h: ET.Element, interactions_dir: Path | None = None) -> HamiltonianSpec:
    """Build a ``HamiltonianSpec`` from a ``<hamiltonian>`` element.

    ``<scale rank="1|2">`` factors are kept; ``<formulascale>`` is ignored, as CoSMo does.
    """
    name = h.get("name", "")
    files = []
    for f in h.findall("file"):
        fname = (f.findtext("name") or "").strip()
        scales = {1: 1.0, 2: 1.0}
        for s in f.findall("scale"):
            rank = s.get("rank")
            if rank not in ("1", "2"):
                msg = f"hamiltonian {name}, file {fname}: <scale> needs rank 1 or 2, got {rank!r}"
                raise ValueError(msg)
            scales[int(rank)] *= float((s.text or "1").strip())
        files.append(
            InteractionFile(
                path=find_interaction_file(fname, interactions_dir),
                spe_scale=scales[1],
                tbme_scale=scales[2],
            )
        )
    if not files:
        msg = f"hamiltonian {name}: only OXBASH <file> entries are supported"
        raise NotImplementedError(msg)
    return HamiltonianSpec(name=name, valence_space=h.get("valencespace", ""), files=tuple(files))


def load_hamiltonian_spec(
    name: str, valence_space: str | None = None, interactions_dir: Path | None = None
) -> HamiltonianSpec:
    """Read Hamiltonian ``name`` (and its file list and scalings) from ``int.xml``."""
    root = _catalogue(interactions_dir)
    for h in root.iter("hamiltonian"):
        if h.get("name") != name:
            continue
        if valence_space is not None and h.get("valencespace", "") != valence_space:
            continue
        return parse_hamiltonian(h, interactions_dir)
    msg = f"hamiltonian {name!r} (valence space {valence_space}) not found in int.xml"
    raise KeyError(msg)


def parse_nucleus(label: str) -> tuple[int, int]:
    """Parse a label such as ``24Mg`` into (Z, N)."""
    digits = "".join(c for c in label if c.isdigit())
    symbol = "".join(c for c in label if c.isalpha())
    if not digits or symbol not in ELEMENTS:
        msg = f"cannot parse nucleus {label!r}"
        raise ValueError(msg)
    a = int(digits)
    z = ELEMENTS.index(symbol)
    return z, a - z
