import numpy as np
import scipy.linalg as sla

from cosmo_gpu.basis import build_mbasis
from cosmo_gpu.reference_csr import hamiltonian_csr, j_expectation, jplus_csr
from tests.systems import build_system


def test_h_symmetric():
    H = build_system("24Mg").H
    assert abs(H - H.T).max() < 1e-12


def test_h_commutes_with_jplus():
    """J+ H = H' J+ with H' on the 2M+2 basis: a strict check of signs and phases."""
    s = build_system("20Ne")
    jp, up = jplus_csr(s.basis)
    H_up = hamiltonian_csr(up, s.inter)
    assert abs(jp @ s.H - H_up @ jp).max() < 1e-11


def test_h_commutes_with_jplus_odd():
    s = build_system("21Ne")
    jp, up = jplus_csr(s.basis)
    H_up = hamiltonian_csr(up, s.inter)
    assert abs(jp @ s.H - H_up @ jp).max() < 1e-11


def test_spins_are_exact():
    s = build_system("20Ne")
    w, v = sla.eigh(s.H.toarray())
    jp, _ = jplus_csr(s.basis)
    J = j_expectation(np.ascontiguousarray(v[:, :40].T), jp, s.basis.m2)
    np.testing.assert_allclose(J, np.round(J), atol=1e-6)


def test_jplus_annihilates_top_m():
    s = build_system("18O")  # 2 neutrons in sd: maximum 2M = 8 (J = 4)
    top = build_mbasis(s.space, s.basis.protons, s.basis.neutrons, 8)
    assert top.dim == 2  # d5/2^2 and d5/2 d3/2, both J = 4
    jp, up = jplus_csr(top)
    assert up.dim == 0
    assert jp.shape == (0, 2)
