"""E2, M1, occupations and spectroscopic factors against CoSMo (XSHLEMB, XSHLAO, XSHLSF).

The references (``tests/data/cosmo_ops_reference.json``) come from fresh CoSMo runs in
``tests/data/cosmo_ops`` (``ops.sb``, ``extract_ops.py``). States are matched by J and order,
and only states CoSMo converged (eigenvector error < 1e-6) are compared. Inputs are rounded to
float32 as CoSMo's parser does. CoSMo prints 6 significant digits.

CoSMo's M1 strengths are uniformly 2.5e-5 (relative) below ours, while ours match KSHELL to
1e-7 (``test_kshell.py``), so the M1 tolerance is 5e-5: a constant inside CoSMo's operator.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as sla

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.lanczos import block_lanczos
from cosmo_gpu.operators import (
    MatrixCache,
    State,
    e2,
    expectation,
    m1,
    number,
    spectroscopic_factor,
    transition_strength,
)
from cosmo_gpu.reference_csr import j_expectation, jplus_csr

from .systems import build_system

REF = json.loads((Path(__file__).parent / "data" / "cosmo_ops_reference.json").read_text())
CONVERGED = 1e-6
HEAVY = {"24Mg", "25Mg"}  # Lanczos on 28k/44k states plus ~100 transitions: run with -m slow


def _param(nucleus: str):
    return pytest.param(nucleus, marks=pytest.mark.slow) if nucleus in HEAVY else nucleus


def _name(j2: int, index: int) -> str:
    js = str(j2 // 2) if j2 % 2 == 0 else f"{j2}/2"
    return f"{js}+({index})"


@cache
def _states(nucleus: str) -> dict[str, State]:
    """Our lowest states, named like CoSMo (J, parity, index within J)."""
    s = build_system(nucleus, cosmo_float32=True)
    nev = len(REF["nuclei"][nucleus]["states"]) + 2
    if s.basis.dim <= 2000:
        _, V = sla.eigh(s.H.toarray())
        vecs = np.ascontiguousarray(V[:, :nev].T)
    else:
        res = block_lanczos(
            lambda X: np.ascontiguousarray((s.H @ X.T).T),
            s.basis.dim,
            nev,
            NumpyBackend(),
            tol=1e-11,
        )
        vecs = res.eigenvectors
    jp, _ = jplus_csr(s.basis)
    j2 = np.rint(2 * j_expectation(vecs, jp, s.basis.m2)).astype(int)
    out: dict[str, State] = {}
    for k in range(nev):
        idx = int(np.sum(j2[: k + 1] == j2[k]))
        out[_name(int(j2[k]), idx)] = State(vecs[k], s.basis, int(j2[k]))
    return out


def _good(nucleus: str) -> set[str]:
    return {
        st["name"]
        for st in REF["nuclei"][nucleus]["states"]
        if st.get("vector_error", 1.0) < CONVERGED
    }


@pytest.mark.parametrize("nucleus", [_param(n) for n in REF["nuclei"]])
@pytest.mark.parametrize("kind", ["E2", "M1"])
def test_transitions_match_cosmo(nucleus: str, kind: str):
    states = _states(nucleus)
    space = build_system(nucleus).space
    if kind == "E2":
        op = e2(space, *REF["e2_charges"])
        rel = 2e-5
    else:
        gl = REF["m1"]["gl"]
        mu = REF["m1"]["mu"]
        op = m1(space, gl[0], 2 * mu[0], gl[1], 2 * mu[1])
        rel = 5e-5
    good = _good(nucleus)
    cache_: MatrixCache = {}
    n = 0
    for ini, fin, b in REF["nuclei"][nucleus][kind]:
        if ini not in good or fin not in good:
            continue
        ours = transition_strength(op, states[ini], states[fin], cache_)
        assert ours == pytest.approx(b, rel=rel, abs=1e-9), (ini, fin)
        n += 1
    assert n > 10


@pytest.mark.parametrize("nucleus", [_param(n) for n in ("22Na", "24Mg", "25Mg")])
def test_occupations_match_cosmo(nucleus: str):
    states = _states(nucleus)
    space = build_system(nucleus).space
    names = [o.name for o in space.orbits]
    for st_name, occ in REF["nuclei"][nucleus]["occupations"].items():
        if st_name not in _good(nucleus):
            continue
        for orb, (z_occ, n_occ) in occ.items():
            a = names.index(orb)
            assert expectation(number(space, a, True), states[st_name]) == pytest.approx(
                z_occ,
                abs=6e-6,  # 6 significant digits on values up to 10
            )
            assert expectation(number(space, a, False), states[st_name]) == pytest.approx(
                n_occ, abs=6e-6
            )


@pytest.mark.slow
def test_spectroscopic_factors_24mg_to_25mg_match_cosmo():
    """S = |<25Mg||a+_j||24Mg>|^2 / (2 J_25 + 1) equals CoSMo's C2S for 25Mg -> 24Mg removal."""
    i_states = _states("24Mg")
    f_states = _states("25Mg")
    space = build_system("24Mg").space
    lab = {(o.n, o.l, o.j2): a for a, o in enumerate(space.orbits)}
    good_i, good_f = _good("24Mg"), _good("25Mg")
    cache_: dict[int, MatrixCache] = {}
    n = 0
    for parent, daughter, nlj, c2s in REF["sf_25Mg_24Mg"]:
        if parent not in good_f or daughter not in good_i:
            continue
        a = lab[tuple(nlj)]
        ours = spectroscopic_factor(
            space, a, False, i_states[daughter], f_states[parent], cache_.setdefault(a, {})
        )
        assert ours == pytest.approx(c2s, rel=2e-5, abs=1e-9), (parent, daughter, nlj)
        n += 1
    assert n > 50
