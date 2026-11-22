"""KSHELL as a second oracle: same ``.snt`` input in both codes, double precision throughout.

The ``.snt`` files were converted from the repo's OXBASH files by KSHELL's own
``nushell2snt.py`` (``tests/data/kshell/make_snt.py``) and the reference energies come from
``kshell.exe`` (``tests/data/kshell/run_kshell.py``, read from the ``.wav`` header at full
precision).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import compute_levels
from cosmo_gpu.interaction import PNCoupledInteraction
from cosmo_gpu.modelspace import parse_nucleus
from cosmo_gpu.snt import read_snt
from cosmo_gpu.truncation import make_truncation

from .systems import build_system

DATA = Path(__file__).parent / "data" / "kshell"
SNT = DATA / "snt"
REFERENCE = DATA / "kshell_reference.json"
KSHELL_TOL = 1e-8  # MeV
TBME_TOL = 1e-7  # nushell2snt prints 8 decimals

RECORDS = json.loads(REFERENCE.read_text()) if REFERENCE.exists() else {}


def _cases(slow_ok: bool = True) -> list[Any]:
    out = []
    for name, rec in sorted(RECORDS.items()):
        marks = [pytest.mark.slow] if rec["slow"] and slow_ok else []
        out.append(pytest.param(rec, id=name, marks=marks))
    return out or [pytest.param(None, id="no-reference", marks=pytest.mark.skip("no KSHELL data"))]


def _valence(rec: dict[str, Any]) -> tuple[int, int]:
    return rec["z_val"], rec["n_val"]


@pytest.mark.parametrize("rec", _cases(slow_ok=False))
def test_dimension_matches_kshell(rec):
    zv, nv = _valence(rec)
    space, _ = read_snt(SNT / rec["snt"], zv + nv)
    hw = tuple(rec["hw"]) if rec["hw"] is not None else None
    trunc = make_truncation(space, zv, nv, hw) if hw is not None else None
    basis = build_basis(space, zv, nv, m2=rec["m2"], parity=rec["parity"], truncation=trunc)
    assert basis.dim == rec["dimension"]


@pytest.mark.parametrize("rec", _cases())
def test_levels_match_kshell(rec):
    assert rec["converged"]
    Z, N = parse_nucleus(rec["nucleus"])
    zv, nv = _valence(rec)
    space, _ = read_snt(SNT / rec["snt"], 0)
    assert (Z - space.core_z, N - space.core_n) == (zv, nv)
    hw = tuple(rec["hw"]) if rec["hw"] is not None else None
    lev = compute_levels(
        rec["nucleus"],
        None,
        None,
        nev=len(rec["energies"]),
        parity=rec["parity"],
        m2=rec["m2"],
        tol=1e-10,
        hw=hw,
        snt=SNT / rec["snt"],
    )
    assert lev.basis.dim == rec["dimension"]
    np.testing.assert_allclose(lev.energies, rec["energies"], atol=KSHELL_TOL, rtol=0)
    np.testing.assert_array_equal(np.rint(2 * lev.spins).astype(int), rec["jj"])


def _assert_equivalent(ours: PNCoupledInteraction, theirs: PNCoupledInteraction):
    """Compare every coupled matrix element, matching orbits by (n, l, 2j)."""
    o1, o2 = ours.space.orbits, theirs.space.orbits
    lab2 = {(o.n, o.l, o.j2): i for i, o in enumerate(o2)}
    m = [lab2[(o.n, o.l, o.j2)] for o in o1]
    np.testing.assert_allclose(theirs.spe_p[m], ours.spe_p, atol=TBME_TOL)
    np.testing.assert_allclose(theirs.spe_n[m], ours.spe_n, atol=TBME_TOL)
    n = len(o1)
    worst = {"pp": 0.0, "nn": 0.0, "pn": 0.0}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    jlo = max(abs(o1[a].j2 - o1[b].j2), abs(o1[c].j2 - o1[d].j2)) // 2
                    jhi = min(o1[a].j2 + o1[b].j2, o1[c].j2 + o1[d].j2) // 2
                    for J in range(jlo, jhi + 1):
                        for name, f1, f2 in (
                            ("pp", ours.get_like, theirs.get_like),
                            ("nn", ours.get_like, theirs.get_like),
                        ):
                            proton = name == "pp"
                            x = f1(proton, a, b, c, d, J)
                            y = f2(proton, m[a], m[b], m[c], m[d], J)
                            worst[name] = max(worst[name], abs(x - y))
                        x = ours.get_pn(a, b, c, d, J)
                        y = theirs.get_pn(m[a], m[b], m[c], m[d], J)
                        worst["pn"] = max(worst["pn"], abs(x - y))
    assert max(worst.values()) < TBME_TOL, worst


@pytest.mark.parametrize(
    ("nucleus", "space", "ham", "snt"),
    [
        ("24Mg", "sd", "usdb", "usdb_from_int.snt"),
        ("44Ti", "fp", "kb3g", "kb3g_from_int.snt"),
        ("33Al", "spsdpf", "fsu9", "fsu9_A33.snt"),
        ("20Ne", "spsdpf", "fsu9", "fsu9_A20.snt"),
    ],
)
def test_isospin_to_pn_matches_kshell_conversion(nucleus, space, ham, snt):
    """Our OXBASH reader + to_pn against KSHELL's nushell2snt on the same files."""
    path = SNT / snt
    if not path.exists():
        pytest.skip(f"{snt} not generated yet")
    ours = build_system(
        nucleus,
        space,
        ham,
        parity=1 if space != "spsdpf" else -1,
        hw=(1, 1) if space == "spsdpf" else None,
    ).coupled.to_pn()
    zv_nv = sum(parse_nucleus(nucleus)) - ours.space.core_z - ours.space.core_n
    _, theirs = read_snt(path, zv_nv)
    _assert_equivalent(ours, theirs)
