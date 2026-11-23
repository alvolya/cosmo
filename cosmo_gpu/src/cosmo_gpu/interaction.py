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

from collections.abc import Callable
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

    def to_pn(self) -> PNCoupledInteraction:
        """Return the same interaction in proton-neutron form.

        Like-particle matrix elements are the ``T = 1`` ones. The pn product-state matrix
        elements are ``V_pn(ab, cd; J) = 1/2 sqrt((1 + d_ab)(1 + d_cd)) sum_T V_JT(ab, cd)``.
        """
        orb = self.space.orbits
        norb = len(orb)
        like = {(a, b, c, d, J): v for (a, b, c, d, J, T), v in self.tbme.items() if T == 1}
        pn: dict[tuple[int, int, int, int, int], float] = {}
        for a in range(norb):
            for b in range(norb):
                for c in range(norb):
                    for d in range(norb):
                        if (a, b) > (c, d):
                            continue
                        jlo = max(abs(orb[a].j2 - orb[b].j2), abs(orb[c].j2 - orb[d].j2)) // 2
                        jhi = min(orb[a].j2 + orb[b].j2, orb[c].j2 + orb[d].j2) // 2
                        norm = 0.5 * sqrt((1 + (a == b)) * (1 + (c == d)))
                        for J in range(jlo, jhi + 1):
                            v = self.get(a, b, c, d, J, 1) + self.get(a, b, c, d, J, 0)
                            if v != 0.0:
                                pn[(a, b, c, d, J)] = norm * v
        return PNCoupledInteraction(
            space=self.space,
            spe_p=self.spe.copy(),
            spe_n=self.spe.copy(),
            pp=like,
            nn=dict(like),
            pn=pn,
        )


def _exchange_phase(ja2: int, jb2: int, J: int, T: int) -> int:
    return -1 if ((ja2 + jb2) // 2 - J - T) % 2 else 1


TBMEKey = tuple[int, int, int, int, int]


@dataclass
class PNCoupledInteraction:
    """Single-particle energies and coupled TBMEs in proton-neutron form (as KSHELL ``.snt``).

    Orbits are positions in ``ModelSpace.orbits``, shared by both species.

    * ``pp`` and ``nn`` hold normalized antisymmetrized ``<ab; J|V|cd; J>`` with canonical keys
      ``a <= b``, ``c <= d`` and ``(a, b) <= (c, d)``.
    * ``pn`` holds product-state ``<a_p b_n; J|V|c_p d_n; J>``, where
      ``|a_p b_n; JM> = sum CG(ja ma jb mb | J M) a+_(a,p) a+_(b,n) |0>`` (proton first), with
      keys ``(a, b) <= (c, d)``.
    """

    space: ModelSpace
    spe_p: NDArray[np.float64]
    spe_n: NDArray[np.float64]
    pp: dict[TBMEKey, float]
    nn: dict[TBMEKey, float]
    pn: dict[TBMEKey, float]

    def get_like(self, proton: bool, a: int, b: int, c: int, d: int, J: int) -> float:
        """Return ``<ab; J|V|cd; J>`` for protons or neutrons, any orbit order."""
        orb = self.space.orbits
        phase = 1
        if a > b:
            a, b = b, a
            phase *= _exchange_phase(orb[a].j2, orb[b].j2, J, 1)
        if c > d:
            c, d = d, c
            phase *= _exchange_phase(orb[c].j2, orb[d].j2, J, 1)
        if (a, b) > (c, d):
            a, b, c, d = c, d, a, b
        table = self.pp if proton else self.nn
        return phase * table.get((a, b, c, d, J), 0.0)

    def get_pn(self, a: int, b: int, c: int, d: int, J: int) -> float:
        """Return ``<a_p b_n; J|V|c_p d_n; J>``."""
        if (a, b) > (c, d):
            a, b, c, d = c, d, a, b
        return self.pn.get((a, b, c, d, J), 0.0)


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
    # the header may continue on further lines, each ending with "-" (OXBASH)
    header = lines[0].split()
    first_tbme = 1
    while header and header[-1] == "-":
        header = header[:-1] + lines[first_tbme].split()
        first_tbme += 1
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
    for ln in lines[first_tbme:]:
        tok = ln.split()
        if len(tok) < 7:
            msg = f"{path.name}: line {ln!r} has {len(tok)} fields; expected i j k l J T V"
            raise ValueError(msg)
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

    ``spe_p`` and ``spe_n`` are single-particle energies per m-orbit. Like-particle terms are
    ``v a+_c1 a+_c2 a_a2 a_a1`` with ``c1 < c2`` and ``a1 < a2``. Proton-neutron terms are
    ``v a+_cp a_ap a+_cn a_an``. Here ``c`` labels created and ``a`` annihilated m-orbits.
    """

    spe_p: NDArray[np.float64]
    spe_n: NDArray[np.float64]
    pp: NDArray[np.int32]
    pp_v: NDArray[np.float64]
    nn: NDArray[np.int32]
    nn_v: NDArray[np.float64]
    pn: NDArray[np.int32]
    pn_v: NDArray[np.float64]


PairGetter = Callable[[int, int, int, int, int], float]


def _pairs_to_m(space: ModelSpace, get: PairGetter, like: bool) -> NDArray[np.float64]:
    """Dense ``W[alpha, beta, gamma, delta] = <alpha beta|V|gamma delta>`` from coupled TBMEs.

    ``get(a, b, c, d, J)`` returns the coupled matrix element. For like particles the pair
    states are normalized, which brings the factor ``sqrt((1 + d_ab)(1 + d_cd))``; pn product
    states need none.
    """
    nm = space.n_morbits
    starts = np.concatenate([[0], np.cumsum([o.degeneracy for o in space.orbits])])
    m2s = [np.arange(-o.j2, o.j2 + 1, 2) for o in space.orbits]
    W = np.zeros((nm, nm, nm, nm))
    norb = space.n_orbits
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

    orb = space.orbits
    for a in range(norb):
        for b in range(norb):
            for c in range(norb):
                for d in range(norb):
                    jlo = max(abs(orb[a].j2 - orb[b].j2), abs(orb[c].j2 - orb[d].j2)) // 2
                    jhi = min(orb[a].j2 + orb[b].j2, orb[c].j2 + orb[d].j2) // 2
                    block = None
                    for J in range(jlo, jhi + 1):
                        v = get(a, b, c, d, J)
                        if v == 0.0:
                            continue
                        if block is None:
                            block = np.zeros((len(m2s[a]), len(m2s[b]), len(m2s[c]), len(m2s[d])))
                        norm = sqrt((1 + (a == b)) * (1 + (c == d))) if like else 1.0
                        block += norm * v * np.multiply.outer(cg_pair(a, b, J), cg_pair(c, d, J))
                    if block is not None:
                        mask = np.equal.outer(m_total(a, b), m_total(c, d))
                        W[
                            starts[a] : starts[a + 1],
                            starts[b] : starts[b + 1],
                            starts[c] : starts[c + 1],
                            starts[d] : starts[d + 1],
                        ] += block * mask
    return W


def like_to_m(inter: PNCoupledInteraction, proton: bool) -> NDArray[np.float64]:
    """Dense antisymmetric m-scheme matrix elements of the pp (or nn) interaction."""

    def get(a: int, b: int, c: int, d: int, J: int) -> float:
        return inter.get_like(proton, a, b, c, d, J)

    return _pairs_to_m(inter.space, get, like=True)


def pn_to_m(inter: PNCoupledInteraction) -> NDArray[np.float64]:
    """Dense m-scheme ``W[pi, nu, pi', nu'] = <pi nu|V|pi' nu'>`` of the pn interaction."""
    return _pairs_to_m(inter.space, inter.get_pn, like=False)


def _like_terms(W: NDArray[np.float64]) -> tuple[NDArray[np.int32], NDArray[np.float64]]:
    nm = W.shape[0]
    upper = np.triu(np.ones((nm, nm), dtype=bool), k=1)
    keep = (np.abs(W) > ZERO_TOL) & upper[:, :, None, None] & upper[None, None, :, :]
    idx = np.argwhere(keep)  # (c1, c2, a1, a2), lexicographic
    rows = np.ascontiguousarray(idx, dtype=np.int32).reshape(-1, 4)
    return rows, W[keep].astype(np.float64)


def to_mscheme(inter: CoupledInteraction | PNCoupledInteraction) -> MSchemeInteraction:
    """Expand a coupled interaction into pp, nn and pn m-scheme term lists."""
    pn_inter = inter.to_pn() if isinstance(inter, CoupledInteraction) else inter
    space = pn_inter.space
    pp, pp_v = _like_terms(like_to_m(pn_inter, proton=True))
    if pn_inter.nn == pn_inter.pp:
        nn, nn_v = pp.copy(), pp_v.copy()
    else:
        nn, nn_v = _like_terms(like_to_m(pn_inter, proton=False))
    Wpn = pn_to_m(pn_inter)
    idx = np.argwhere(np.abs(Wpn) > ZERO_TOL)  # (cp, cn, ap, an)
    pn = np.ascontiguousarray(idx[:, [0, 2, 1, 3]], dtype=np.int32)  # -> (cp, ap, cn, an)
    pn_v = Wpn[idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]].astype(np.float64)
    return MSchemeInteraction(
        spe_p=pn_inter.spe_p[space.morb_orbit].astype(np.float64),
        spe_n=pn_inter.spe_n[space.morb_orbit].astype(np.float64),
        pp=pp,
        pp_v=pp_v,
        nn=nn,
        nn_v=nn_v,
        pn=pn,
        pn_v=pn_v,
    )
