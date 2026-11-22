from math import sqrt

import numpy as np
import pytest

from cosmo_gpu.interaction import _coupled_to_m, clebsch_gordan, read_interaction
from cosmo_gpu.modelspace import load_hamiltonian_spec
from tests.systems import build_system


def test_cg_known_values():
    assert clebsch_gordan(1, 1, 1, -1, 0, 0) == pytest.approx(1 / sqrt(2))
    assert clebsch_gordan(1, -1, 1, 1, 0, 0) == pytest.approx(-1 / sqrt(2))
    assert clebsch_gordan(2, 2, 2, 0, 4, 2) == pytest.approx(1 / sqrt(2))
    assert clebsch_gordan(5, 5, 5, 5, 10, 10) == pytest.approx(1.0)
    assert clebsch_gordan(3, 1, 1, 1, 2, 4) == 0.0


@pytest.mark.parametrize(("j1", "j2"), [(1, 1), (3, 5), (5, 5), (7, 3)])
def test_cg_orthogonality(j1, j2):
    for M in range(-(j1 + j2), j1 + j2 + 1, 2):
        js = range(abs(j1 - j2), j1 + j2 + 1, 2)
        ms = [m1 for m1 in range(-j1, j1 + 1, 2) if abs(M - m1) <= j2]
        C = np.array(
            [[clebsch_gordan(j1, m1, j2, M - m1, J, M) for m1 in ms] for J in js if abs(M) <= J]
        )
        np.testing.assert_allclose(C @ C.T, np.eye(C.shape[0]), atol=1e-13)
        np.testing.assert_allclose(C.T @ C, np.eye(C.shape[1]), atol=1e-13)


def test_mass_scaling_matches_cosmo(sd_space):
    # CoSMo prints "Lines in interaction file :63 scale: 0.906149" for 25Mg (9 valence)
    spec = load_hamiltonian_spec("usdb", "sd")
    raw = read_interaction(spec, sd_space, 2)  # (18/18)^0.3 = 1
    scaled = read_interaction(spec, sd_space, 9)
    key = next(iter(raw.tbme))
    assert scaled.tbme[key] / raw.tbme[key] == pytest.approx(0.906149, abs=1e-6)
    np.testing.assert_array_equal(raw.spe, scaled.spe)


def test_exchange_symmetry(sd_space):
    inter = build_system("20Ne").coupled
    for J in range(6):
        for T in (0, 1):
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        for d in range(3):
                            v = inter.get(a, b, c, d, J, T)
                            assert inter.get(c, d, a, b, J, T) == v
                            ph = (-1) ** (
                                (sd_space.orbits[a].j2 + sd_space.orbits[b].j2) // 2 - J - T
                            )
                            assert inter.get(b, a, c, d, J, T) == pytest.approx(ph * v)


def test_mscheme_antisymmetry_and_hermiticity():
    inter = build_system("20Ne").coupled
    W = _coupled_to_m(inter, 1)
    np.testing.assert_allclose(W, -W.transpose(1, 0, 2, 3), atol=1e-12)
    np.testing.assert_allclose(W, -W.transpose(0, 1, 3, 2), atol=1e-12)
    np.testing.assert_allclose(W, W.transpose(2, 3, 0, 1), atol=1e-12)
