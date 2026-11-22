"""OXBASH isospin-format interactions and their conversion to pn m-scheme matrix elements.

Conventions
-----------
* Coupled matrix elements ``V_JT(ab, cd) = <ab; JT | V | cd; JT>`` are normalized and
  antisymmetrized, with the pair state ``|ab; JM> = N_ab sum CG(ja ma jb mb | J M) a+_a a+_b |0>``
  and ``N_ab = 1 / sqrt(1 + delta_ab)``.
* Exchanging the orbits of a pair gives ``|ba; JT> = (-1)^(ja + jb - J - T) |ab; JT>``.
* The m-scheme matrix element is ``<alpha beta | V | gamma delta> = <0| a_beta a_alpha V
  a+_gamma a+_delta |0>``, and the two-body operator is
  ``sum_{alpha<beta, gamma<delta} <alpha beta|V|gamma delta> a+_alpha a+_beta a_delta a_gamma``
  for like particles and ``sum <pi nu|V|pi' nu'> a+_pi a_pi' a+_nu a_nu'`` for protons and
  neutrons (proton operators to the left of neutron operators).
* Isospin-format files carry no Coulomb term, so pp and nn matrix elements are equal.

Mass scaling follows CoSMo (``sm/sphericalint.cxx``): a negative line count in the header means
the last three header numbers are ``Acore A0 p`` and every TBME of that file is multiplied by
``(A0 / (n_valence + Acore))^p``. CoSMo's ``<scale rank="1|2">`` elements multiply the s.p.
energies and TBMEs of a file.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial, sqrt
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.modelspace import HamiltonianSpec, ModelSpace

ZERO_TOL = 1e-13


def clebsch_gordan(j1: int, m1: int, j2: int, m2: int, j: int, m: int) -> float:
    """Clebsch-Gordan coefficient ``<j1 m1 j2 m2 | j m>`` with all arguments doubled."""
    if m1 + m2 != m or not abs(j1 - j2) <= j <= j1 + j2 or (j1 + j2 + j) % 2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return 0.0
    if (j1 + m1) % 2 or (j2 + m2) % 2 or (j + m) % 2:
        return 0.0
    a = (j1 + j2 - j) // 2
    b = (j1 - j2 + j) // 2
    c = (-j1 + j2 + j) // 2
    pref = (j + 1) * factorial(a) * factorial(b) * factorial(c) / factorial((j1 + j2 + j) // 2 + 1)
    pref *= (
        factorial((j1 + m1) // 2)
        * factorial((j1 - m1) // 2)
        * factorial((j2 + m2) // 2)
        * factorial((j2 - m2) // 2)
        * factorial((j + m) // 2)
        * factorial((j - m) // 2)
    )
    total = 0.0
    for k in range(0, a + 1):
        d = [
            k,
            a - k,
            (j1 - m1) // 2 - k,
            (j2 + m2) // 2 - k,
            (j - j2 + m1) // 2 + k,
            (j - j1 - m2) // 2 + k,
        ]
        if min(d) < 0:
            continue
        den = 1
        for x in d:
            den *= factorial(x)
        total += (-1) ** k / den
    return sqrt(pref) * total


# ------------------------------------------------------------------ coupled (JT) interaction


@dataclass
class CoupledInteraction:
    """Single-particle energies and canonical coupled TBMEs in isospin form.

    Orbits are referred to by their position in ``ModelSpace.orbits``. Keys of ``tbme`` are
    canonical: ``a <= b``, ``c <= d`` and ``(a, b) <= (c, d)``.
    """

    space: ModelSpace
    spe: NDArray[np.float64]
    tbme: dict[tuple[int, int, int, int, int, int], float]

    def get(self, a: int, b: int, c: int, d: int, J: int, T: int) -> float:
        """Return ``V_JT(ab, cd)`` for any ordering of the orbits (``J`` not doubled)."""
        phase = 1
        orb = self.space.orbits
        if a > b:
            a, b = b, a
            phase *= _exchange_phase(orb[a].j2, orb[b].j2, J, T)
        if c > d:
            c, d = d, c
            phase *= _exchange_phase(orb[c].j2, orb[d].j2, J, T)
        if (a, b) > (c, d):
            a, b, c, d = c, d, a, b
        return phase * self.tbme.get((a, b, c, d, J, T), 0.0)


def _exchange_phase(ja2: int, jb2: int, J: int, T: int) -> int:
    return -1 if ((ja2 + jb2) // 2 - J - T) % 2 else 1


def _read_oxbash(
    path: Path,
    space: ModelSpace,
    n_valence: int,
    spe_scale: float,
    tbme_scale: float,
    cosmo_float32: bool = False,
) -> tuple[NDArray[np.float64], dict[tuple[int, int, int, int, int, int], float]]:
    num = _parse_float32 if cosmo_float32 else float
    lines = [ln.strip() for ln in path.read_text().splitlines()]
    lines = [ln for ln in lines if ln and ln[0] not in "!#$%"]
    header = lines[0].split()
    count = int(float(header[0]))
    values = [num(x) for x in header[1:]]
    mass_factor = 1.0
    if count < 0:
        a_core, a_ref, power = values[-3:]
        values = values[:-3]
        mass_factor = (a_ref / (n_valence + a_core)) ** power
    spe = np.zeros(space.n_orbits)
    for i, e in enumerate(values):
        spe[space.orbit_position(i + 1)] = e * spe_scale
    tbme: dict[tuple[int, int, int, int, int, int], float] = {}
    scale = tbme_scale * mass_factor
    orb = space.orbits
    for ln in lines[1:]:
        tok = ln.split()
        if len(tok) < 7:
            continue
        a, b, c, d = (space.orbit_position(int(x)) for x in tok[:4])
        J, T = int(tok[4]), int(tok[5])
        v = num(tok[6]) * scale
        if a > b:
            a, b = b, a
            v *= _exchange_phase(orb[a].j2, orb[b].j2, J, T)
        if c > d:
            c, d = d, c
            v *= _exchange_phase(orb[c].j2, orb[d].j2, J, T)
        if (a, b) > (c, d):
            a, b, c, d = c, d, a, b
        key = (a, b, c, d, J, T)
        if key in tbme:
            if abs(tbme[key] - v) > 1e-6 * max(1.0, abs(v)):
                msg = f"{path.name}: inconsistent duplicate matrix element {key}"
                raise ValueError(msg)
            continue
        tbme[key] = v
    return spe, tbme


def _parse_float32(text: str) -> float:
    return float(np.float32(text))


def read_interaction(
    spec: HamiltonianSpec, space: ModelSpace, n_valence: int, cosmo_float32: bool = False
) -> CoupledInteraction:
    """Read and sum all OXBASH files of a Hamiltonian, with CoSMo's scalings applied.

    With ``cosmo_float32`` every number in the files is rounded to single precision first, as
    CoSMo's ``std::stof`` parser does. This reproduces CoSMo energies to about 1e-7 MeV instead
    of about 1e-6 MeV and is meant for code-code comparisons only.
    """
    spe = np.zeros(space.n_orbits)
    tbme: dict[tuple[int, int, int, int, int, int], float] = {}
    for f in spec.files:
        s, t = _read_oxbash(f.path, space, n_valence, f.spe_scale, f.tbme_scale, cosmo_float32)
        spe += s
        for k, v in t.items():
            tbme[k] = tbme.get(k, 0.0) + v
    return CoupledInteraction(space=space, spe=spe, tbme=tbme)


# ------------------------------------------------------------------ m-scheme interaction


@dataclass(frozen=True)
class MSchemeInteraction:
    """Hamiltonian in pn m-scheme form, indexed by m-orbit.

    Like-particle terms are ``v a+_c1 a+_c2 a_a2 a_a1`` with ``c1 < c2`` and ``a1 < a2``.
    Proton-neutron terms are ``v a+_cp a_ap a+_cn a_an``. Here ``c`` labels created and ``a``
    annihilated m-orbits.
    """

    spe: NDArray[np.float64]
    pp: NDArray[np.int32]
    pp_v: NDArray[np.float64]
    nn: NDArray[np.int32]
    nn_v: NDArray[np.float64]
    pn: NDArray[np.int32]
    pn_v: NDArray[np.float64]


def _coupled_to_m(inter: CoupledInteraction, T: int) -> NDArray[np.float64]:
    """Dense ``W[alpha, beta, gamma, delta]`` from the ``T`` channel, any species labels."""
    space = inter.space
    nm = space.n_morbits
    starts = np.concatenate([[0], np.cumsum([o.degeneracy for o in space.orbits])])
    m2s = [np.arange(-o.j2, o.j2 + 1, 2) for o in space.orbits]
    W = np.zeros((nm, nm, nm, nm))
    norb = space.n_orbits
    jmax2 = 2 * max(o.j2 for o in space.orbits)
    cg_cache: dict[tuple[int, int, int], NDArray[np.float64]] = {}

    def cg_pair(a: int, b: int, J: int) -> NDArray[np.float64]:
        key = (a, b, J)
        if key not in cg_cache:
            ja, jb = space.orbits[a].j2, space.orbits[b].j2
            cg_cache[key] = np.array(
                [
                    [clebsch_gordan(ja, int(ma), jb, int(mb), 2 * J, int(ma + mb)) for mb in m2s[b]]
                    for ma in m2s[a]
                ]
            )
        return cg_cache[key]

    def m_total(a: int, b: int) -> NDArray[np.int64]:
        return np.asarray(np.add.outer(m2s[a], m2s[b]), dtype=np.int64)

    for a in range(norb):
        for b in range(norb):
            for c in range(norb):
                for d in range(norb):
                    block = np.zeros((len(m2s[a]), len(m2s[b]), len(m2s[c]), len(m2s[d])))
                    mask = None
                    for J in range(0, jmax2 // 2 + 1):
                        v = inter.get(a, b, c, d, J, T)
                        if v == 0.0:
                            continue
                        if mask is None:
                            mask = np.equal.outer(m_total(a, b), m_total(c, d))
                        norm = sqrt((1 + (a == b)) * (1 + (c == d)))
                        block += norm * v * np.multiply.outer(cg_pair(a, b, J), cg_pair(c, d, J))
                    if mask is not None:
                        W[
                            starts[a] : starts[a + 1],
                            starts[b] : starts[b + 1],
                            starts[c] : starts[c + 1],
                            starts[d] : starts[d + 1],
                        ] += block * mask
    return W


def to_mscheme(inter: CoupledInteraction) -> MSchemeInteraction:
    """Expand an isospin-coupled interaction into pp, nn and pn m-scheme term lists."""
    space = inter.space
    nm = space.n_morbits
    spe = inter.spe[space.morb_orbit].astype(np.float64)
    W1 = _coupled_to_m(inter, 1)
    W0 = _coupled_to_m(inter, 0)
    Wpn = 0.5 * (W1 + W0)

    iu = [(i, j) for i in range(nm) for j in range(i + 1, nm)]
    pp_rows: list[tuple[int, int, int, int]] = []
    pp_vals: list[float] = []
    for c1, c2 in iu:
        for a1, a2 in iu:
            v = W1[c1, c2, a1, a2]
            if abs(v) > ZERO_TOL:
                pp_rows.append((c1, c2, a1, a2))
                pp_vals.append(v)
    idx = np.argwhere(np.abs(Wpn) > ZERO_TOL)  # (cp, cn, ap, an)
    pn = np.ascontiguousarray(idx[:, [0, 2, 1, 3]], dtype=np.int32)  # -> (cp, ap, cn, an)
    pn_v = Wpn[idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]].astype(np.float64)
    pp = np.array(pp_rows, dtype=np.int32).reshape(-1, 4)
    pp_v = np.array(pp_vals, dtype=np.float64)
    return MSchemeInteraction(
        spe=spe, pp=pp, pp_v=pp_v, nn=pp.copy(), nn_v=pp_v.copy(), pn=pn, pn_v=pn_v
    )
