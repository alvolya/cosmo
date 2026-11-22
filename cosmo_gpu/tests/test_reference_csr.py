import numpy as np
import pytest
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


def _brute_force_h(s):
    """Dense H from the flat term lists with plain Python bit operations (independent check)."""
    b = s.basis
    P, N = b.protons, b.neutrons
    states = []
    for B in range(b.n_blocks):
        pd = P.dets[P.sub_start[b.block_psub[B]] : P.sub_start[b.block_psub[B] + 1]]
        nd = N.dets[N.sub_start[b.block_nsub[B]] : N.sub_start[b.block_nsub[B] + 1]]
        states += [(int(p), int(n)) for p in pd for n in nd]
    index = {st: i for i, st in enumerate(states)}

    def apply(s0, ops):
        sign = 1
        for kind, k in ops:
            bit = 1 << int(k)
            occupied = bool(s0 & bit)
            if occupied != (kind == "a"):
                return None, 0
            if bin(s0 & (bit - 1)).count("1") % 2:
                sign = -sign
            s0 ^= bit
        return s0, sign

    inter = s.inter
    H = np.zeros((len(states), len(states)))
    for j, (p, n) in enumerate(states):
        H[j, j] += sum(inter.spe[k] for k in range(64) if (p >> k) & 1)
        H[j, j] += sum(inter.spe[k] for k in range(64) if (n >> k) & 1)
        for (c1, c2, a1, a2), v in zip(inter.pp, inter.pp_v, strict=True):
            ops = [("a", a1), ("a", a2), ("c", c2), ("c", c1)]
            p2, g = apply(p, ops)
            if g and (p2, n) in index:
                H[index[(p2, n)], j] += g * v
            n2, g = apply(n, ops)
            if g and (p, n2) in index:
                H[index[(p, n2)], j] += g * v
        for (cp, ap, cn, an), v in zip(inter.pn, inter.pn_v, strict=True):
            p2, g1 = apply(p, [("a", ap), ("c", cp)])
            n2, g2 = apply(n, [("a", an), ("c", cn)])
            if g1 and g2 and (p2, n2) in index:
                H[index[(p2, n2)], j] += g1 * g2 * v
    return H


def test_indexed_kernel_matches_brute_force():
    s = build_system("20Ne")
    np.testing.assert_allclose(s.H.toarray(), _brute_force_h(s), atol=1e-12)


@pytest.mark.parametrize(("parity", "hw"), [(1, (0, 0)), (-1, (1, 1))])
def test_h_commutes_with_jplus_cross_shell(parity, hw):
    """Truncated psdpf basis: tests phases of parity-changing pn terms and of truncation."""
    s = build_system("18O", "spsdpf", "fsu9", parity=parity, hw=hw)
    jp, up = jplus_csr(s.basis)
    H_up = hamiltonian_csr(up, s.inter)
    assert abs(jp @ s.H - H_up @ jp).max() < 1e-10
    assert abs(s.H - s.H.T).max() < 1e-12
