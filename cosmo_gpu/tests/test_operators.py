"""Sum rules and internal consistency of one-body operators (all final states, dense)."""

from __future__ import annotations

from functools import cache
from math import sqrt

import numpy as np
import pytest
import scipy.linalg as sla

from cosmo_gpu.basis import build_mbasis, build_species_basis
from cosmo_gpu.interaction import read_interaction, to_mscheme
from cosmo_gpu.modelspace import load_hamiltonian_spec, load_model_space
from cosmo_gpu.operators import (
    MatrixCache,
    State,
    coupled_density,
    creation_matrix,
    e2,
    expectation,
    gamow_teller,
    ho_radial,
    m1,
    number,
    sp_reduced,
    spectroscopic_factor,
    tensor_rme,
    total_strength,
    transition_strength,
)
from cosmo_gpu.reference_csr import hamiltonian_csr, j_expectation, jplus_csr

SUM_TOL = 1e-10


@cache
def _space():
    return load_model_space("sd")


@cache
def _dense(zv: int, nv: int, m2: int | None = None) -> tuple[State, ...]:
    """All eigenstates of USDB for (zv, nv) valence nucleons at the lowest M, with J."""
    space = _space()
    inter = to_mscheme(read_interaction(load_hamiltonian_spec("usdb", "sd"), space, zv + nv))
    P = build_species_basis(space, zv)
    N = build_species_basis(space, nv)
    m2 = (zv + nv) % 2 if m2 is None else m2
    basis = build_mbasis(space, P, N, m2, parity=1)
    H = hamiltonian_csr(basis, inter).toarray()
    _, V = sla.eigh(H)
    jp, _ = jplus_csr(basis)
    J = j_expectation(np.ascontiguousarray(V.T), jp, basis.m2)
    j2 = np.rint(2 * J).astype(int)
    assert np.allclose(2 * J, j2, atol=1e-6), "degenerate states with mixed J"
    return tuple(State(np.ascontiguousarray(V[:, k]), basis, int(j2[k])) for k in range(V.shape[1]))


def _species(z: int, n: int):
    space = _space()
    return build_species_basis(space, z), build_species_basis(space, n)


def test_radial_integrals_match_closed_form():
    # <0d|r^2|1s> = -sqrt(10), <0d|r^2|0d> = 7/2 (b = 1), as KSHELL's radius_power
    assert ho_radial(0, 2, 1, 0, 2) == pytest.approx(-sqrt(10.0), abs=1e-12)
    assert ho_radial(0, 2, 0, 2, 2) == pytest.approx(3.5, abs=1e-12)
    assert ho_radial(1, 0, 1, 0, 0) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("which", ["E2", "M1"])
@pytest.mark.parametrize("initial", [0, 1])
def test_newsr_e2_m1_20ne(which: str, initial: int):
    space = _space()
    op = e2(space, 1.5, 0.5) if which == "E2" else m1(space, 1.0, 5.586, 0.0, -3.826)
    states = _dense(2, 2)
    i = states[initial]
    cache_ = MatrixCache()
    total = sum(transition_strength(op, i, f, cache_) for f in states)
    assert total == pytest.approx(total_strength(op, i), rel=SUM_TOL, abs=SUM_TOL)


@pytest.mark.parametrize(("zv", "nv"), [(1, 3), (2, 3)])
def test_ikeda_sum_rule(zv: int, nv: int):
    """S(n->p) - S(p->n) = 3 (N - Z) in a space closed under spin-orbit partners."""
    space = _space()
    i = _dense(zv, nv)[0]
    s_minus = total_strength(gamow_teller(space, "n->p"), i, *_species(zv + 1, nv - 1))
    s_plus = total_strength(gamow_teller(space, "p->n"), i, *_species(zv - 1, nv + 1))
    assert s_minus - s_plus == pytest.approx(3 * (nv - zv), abs=SUM_TOL)


def test_gt_strengths_sum_to_total_20ne_to_20na():
    space = _space()
    op = gamow_teller(space, "n->p")
    i = _dense(2, 2)[0]
    finals = _dense(3, 1, m2=0)
    cache_ = MatrixCache()
    total = sum(transition_strength(op, i, f, cache_) for f in finals)
    assert total == pytest.approx(total_strength(op, i, *_species(3, 1)), rel=SUM_TOL)


@pytest.mark.parametrize("orbit", [0, 1, 2])
def test_macfarlane_french_18o_to_19o(orbit: int):
    """sum_f (2J_f+1)/(2J_i+1) S = holes in orbit j for neutron addition."""
    space = _space()
    i = _dense(0, 2)[0]
    finals = _dense(0, 3)
    cache_ = MatrixCache()
    total = sum(
        (f.j2 + 1) / (i.j2 + 1) * spectroscopic_factor(space, orbit, False, i, f, cache_)
        for f in finals
    )
    holes = (space.orbits[orbit].j2 + 1) - expectation(number(space, orbit, False), i)
    assert total == pytest.approx(holes, abs=SUM_TOL)
    # the same from sum_m ||a+_(j m) psi||^2
    direct = 0.0
    for x in np.flatnonzero(space.morb_orbit == orbit):
        m2 = int(space.morb_m2[x])
        dst = build_mbasis(space, i.basis.protons, finals[0].basis.neutrons, i.basis.m2 + m2, 1)
        v = creation_matrix(space, orbit, m2, False, dst, i.basis) @ i.vector
        direct += float(v @ v)
    assert direct == pytest.approx(holes, abs=SUM_TOL)


def test_occupations_sum_to_valence_numbers():
    space = _space()
    s = _dense(2, 3)[0]
    zp = sum(expectation(number(space, a, True), s) for a in range(space.n_orbits))
    zn = sum(expectation(number(space, a, False), s) for a in range(space.n_orbits))
    assert (zp, zn) == pytest.approx((2.0, 3.0), abs=1e-12)


def test_e2_from_one_body_densities():
    """<f||T||i> = sum_ab <a||t||b> <f||[a+_a ~a_b]^k||i> / sqrt(2k+1) (Edmonds)."""
    space = _space()
    op = e2(space, 1.5, 0.5)
    states = _dense(2, 2)
    i = states[0]
    f = next(s for s in states if s.j2 == 4)
    direct = tensor_rme(op, i, f)
    via = 0.0
    for proton in (True, False):
        for a in range(space.n_orbits):
            for b in range(space.n_orbits):
                t_ab = sp_reduced(space, op, a, b, proton)
                if t_ab == 0.0:
                    continue
                rho = tensor_rme(coupled_density(space, a, b, 4, proton), i, f)
                via += t_ab * rho / sqrt(5)
    assert via == pytest.approx(direct, rel=1e-10)


def test_rme_independent_of_m():
    """The same reduced matrix element from states computed at different M."""
    space = _space()
    op = e2(space, 1.5, 0.5)
    states0 = _dense(2, 2, m2=0)
    i0 = states0[0]
    f0 = next(s for s in states0 if s.j2 == 4)
    states2 = _dense(2, 2, m2=2)  # J >= 1 only
    f2 = next(s for s in states2 if s.j2 == 4)
    # same state up to sign: compare B, which is sign-free
    b0 = transition_strength(op, i0, f0)
    b2 = transition_strength(op, i0, f2)
    assert b2 == pytest.approx(b0, rel=1e-10)
