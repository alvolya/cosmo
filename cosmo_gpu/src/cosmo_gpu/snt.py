"""KSHELL ``.snt`` interaction files: proton-neutron form, read and write.

Layout (comments start with ``!`` or ``#``)::

    n_p_orbits n_n_orbits core_Z core_N
    k  n  l  2j  tz            (one line per orbit, protons tz = -1 first, then neutrons tz = +1)
    n_1b method1 [hw]
    i  j  <i|H|j>              (one-body terms)
    n_2b method2 [A0 p]
    i  j  k  l  J  V           (TBMEs)

TBMEs are ``<ij; J|V|kl; J>``. For two protons or two neutrons they are normalized and
antisymmetrized; for a proton and a neutron (``i`` proton, ``j`` neutron) they are product-state
matrix elements, as in ``PNCoupledInteraction``. ``method2 = 1`` scales every TBME by
``(A / A0)^p`` with ``A`` the mass number (core + valence), as KSHELL's
``operator_jscheme.f90`` does. One-body methods 0 and 1 carry no mass dependence.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.interaction import PNCoupledInteraction, TBMEKey, _exchange_phase
from cosmo_gpu.modelspace import ModelSpace, Orbit

L_NAMES = "spdfghijklmno"
SYMMETRY_TOL = 1e-6


def _data_lines(path: Path) -> list[list[str]]:
    rows = []
    for raw in path.read_text().splitlines():
        line = raw
        for mark in "!#":
            line = line.split(mark, 1)[0]
        tok = line.split()
        if tok:
            rows.append(tok)
    return rows


def _orbit_name(n: int, l: int, j2: int) -> str:
    return f"{n}{L_NAMES[l]}{j2}"


class _Rows:
    """Sequential access to the data lines of a file."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.rows = _data_lines(path)
        self.pos = 0

    def take(self) -> list[str]:
        if self.pos >= len(self.rows):
            msg = f"{self.path.name}: unexpected end of file"
            raise ValueError(msg)
        self.pos += 1
        return self.rows[self.pos - 1]


Where = dict[int, tuple[bool, int]]


def read_snt(path: Path, n_valence: int) -> tuple[ModelSpace, PNCoupledInteraction]:
    """Read a KSHELL ``.snt`` file for a nucleus with ``n_valence`` valence nucleons.

    Protons and neutrons must have the same set of orbits. The model space takes the proton
    orbit order and is named after the file.

    Raises
    ------
    NotImplementedError
        For off-diagonal one-body terms, or a mass-dependence method other than 0 or 1.
    ValueError
        For malformed files or inconsistent duplicate matrix elements.
    """
    rows = _Rows(path)
    space, where = _read_orbits(rows)
    spe_p, spe_n = _read_one_body(rows, space, where)
    pp, nn, pn = _read_two_body(rows, space, where, n_valence)
    return space, PNCoupledInteraction(space, spe_p, spe_n, pp, nn, pn)


def _read_orbits(rows: _Rows) -> tuple[ModelSpace, Where]:
    """Read the orbit header; map file orbit k (1-based) to (is proton, position in space)."""
    n_p, n_n, core_z, core_n = (int(x) for x in rows.take()[:4])
    proton_nlj: list[tuple[int, int, int]] = []
    neutron_nlj: list[tuple[int, int, int]] = []
    for _ in range(n_p + n_n):
        _k, n, l, j2, tz = (int(x) for x in rows.take()[:5])
        (proton_nlj if tz == -1 else neutron_nlj).append((n, l, j2))
    if len(proton_nlj) != n_p or sorted(proton_nlj) != sorted(neutron_nlj):
        msg = f"{rows.path.name}: protons and neutrons must share one set of orbits"
        raise NotImplementedError(msg)
    orbits = tuple(
        Orbit(name=_orbit_name(n, l, j2), n=n, l=l, j2=j2, index=i + 1)
        for i, (n, l, j2) in enumerate(proton_nlj)
    )
    space = ModelSpace(name=rows.path.stem, orbits=orbits, core_z=core_z, core_n=core_n)
    where: Where = {k + 1: (True, k) for k in range(n_p)}
    for k, nlj in enumerate(neutron_nlj):
        where[n_p + k + 1] = (False, proton_nlj.index(nlj))
    return space, where


def _read_one_body(
    rows: _Rows, space: ModelSpace, where: Where
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    head = rows.take()
    n_1b, method1 = int(head[0]), int(head[1])
    if method1 not in (0, 1):
        msg = f"{rows.path.name}: one-body method {method1} is not supported"
        raise NotImplementedError(msg)
    spe_p = np.zeros(space.n_orbits)
    spe_n = np.zeros(space.n_orbits)
    for _ in range(n_1b):
        tok = rows.take()
        i, j, v = int(tok[0]), int(tok[1]), float(tok[2])
        if i != j:
            if v != 0.0:
                msg = f"{rows.path.name}: off-diagonal one-body term ({i}, {j})"
                raise NotImplementedError(msg)
            continue
        proton, a = where[i]
        (spe_p if proton else spe_n)[a] += v
    return spe_p, spe_n


def _read_two_body(
    rows: _Rows, space: ModelSpace, where: Where, n_valence: int
) -> tuple[dict[TBMEKey, float], dict[TBMEKey, float], dict[TBMEKey, float]]:
    head = rows.take()
    n_2b, method2 = int(head[0]), int(head[1])
    if method2 == 0:
        factor = 1.0
    elif method2 == 1:
        a0, power = float(head[2]), float(head[3])
        factor = ((space.core_z + space.core_n + n_valence) / a0) ** power
    else:
        msg = f"{rows.path.name}: two-body method {method2} is not supported"
        raise NotImplementedError(msg)
    pp: dict[TBMEKey, float] = {}
    nn: dict[TBMEKey, float] = {}
    pn: dict[TBMEKey, float] = {}
    for _ in range(n_2b):
        tok = rows.take()
        idx = [where[int(x)] for x in tok[:4]]
        J, v = int(tok[4]), float(tok[5]) * factor
        species = tuple(p for p, _a in idx)
        a, b, c, d = (o for _p, o in idx)
        if species == (True, True, True, True):
            _store(pp, *_canonical_like(space, a, b, c, d, J, v), rows.path)
        elif species == (False, False, False, False):
            _store(nn, *_canonical_like(space, a, b, c, d, J, v), rows.path)
        elif species[0] != species[1] and species[2] != species[3]:
            _store(pn, *_canonical_pn(space, species, a, b, c, d, J, v), rows.path)
        else:
            msg = f"{rows.path.name}: TBME {tok[:5]} does not conserve charge"
            raise ValueError(msg)
    return pp, nn, pn


def _canonical_pn(
    space: ModelSpace,
    species: tuple[bool, ...],
    a: int,
    b: int,
    c: int,
    d: int,
    J: int,
    v: float,
) -> tuple[TBMEKey, float]:
    """Put the proton first in each pair: ``|b_n a_p; J> = -(-1)^(ja+jb-J) |a_p b_n; J>``."""
    orb = space.orbits
    if not species[0]:
        v *= -_exchange_phase(orb[a].j2, orb[b].j2, J, 0)
        a, b = b, a
    if not species[2]:
        v *= -_exchange_phase(orb[c].j2, orb[d].j2, J, 0)
        c, d = d, c
    if (a, b) > (c, d):
        a, b, c, d = c, d, a, b
    return (a, b, c, d, J), v


def _canonical_like(
    space: ModelSpace, a: int, b: int, c: int, d: int, J: int, v: float
) -> tuple[TBMEKey, float]:
    orb = space.orbits
    if a > b:
        a, b = b, a
        v *= _exchange_phase(orb[a].j2, orb[b].j2, J, 1)
    if c > d:
        c, d = d, c
        v *= _exchange_phase(orb[c].j2, orb[d].j2, J, 1)
    if (a, b) > (c, d):
        a, b, c, d = c, d, a, b
    return (a, b, c, d, J), v


def _store(table: dict[TBMEKey, float], key: TBMEKey, v: float, path: Path) -> None:
    if key in table:
        if abs(table[key] - v) > SYMMETRY_TOL * max(1.0, abs(v)):
            msg = f"{path.name}: inconsistent duplicate matrix element {key}"
            raise ValueError(msg)
        return
    table[key] = v


def write_snt(path: Path, inter: PNCoupledInteraction, comment: str = "") -> None:
    """Write ``inter`` as a ``.snt`` file with no mass dependence (method 0).

    Values are written with 17 significant digits, so reading the file back is exact.
    """
    space = inter.space
    no = space.n_orbits
    out = [f"! {line}" for line in comment.splitlines()]
    out.append(f" {no:3d} {no:3d}   {space.core_z:3d} {space.core_n:3d}")
    for tz, shift in ((-1, 0), (1, no)):
        for i, o in enumerate(space.orbits):
            out.append(f"{i + 1 + shift:5d}   {o.n:3d} {o.l:3d} {o.j2:3d} {tz:3d}")
    out.append(f" {2 * no:3d} {0:3d}")
    for shift, spe in ((0, inter.spe_p), (no, inter.spe_n)):
        out.extend(f"{i + 1 + shift:3d} {i + 1 + shift:3d} {spe[i]:.17g}" for i in range(no))
    lines = []
    for table, sa, sb in ((inter.pp, 0, 0), (inter.nn, no, no), (inter.pn, 0, no)):
        for (a, b, c, d, J), v in sorted(table.items()):
            lines.append(
                f"{a + 1 + sa:3d} {b + 1 + sb:3d} {c + 1 + sa:3d} {d + 1 + sb:3d} {J:3d} {v:.17g}"
            )
    out.append(f" {len(lines):10d} {0:3d}")
    out.extend(lines)
    path.write_text("\n".join(out) + "\n")
