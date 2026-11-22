import numpy as np
import pytest

from cosmo_gpu.modelspace import find_interaction_file, load_hamiltonian_spec, parse_nucleus


def test_sd_orbits(sd_space):
    assert [o.name for o in sd_space.orbits] == ["0d3", "0d5", "1s1"]
    assert [o.j2 for o in sd_space.orbits] == [3, 5, 1]
    assert (sd_space.core_z, sd_space.core_n) == (8, 8)
    assert sd_space.n_morbits == 12
    np.testing.assert_array_equal(sd_space.morb_m2[:4], [-3, -1, 1, 3])
    np.testing.assert_array_equal(sd_space.morb_orbit, [0] * 4 + [1] * 6 + [2] * 2)


def test_parse_nucleus():
    assert parse_nucleus("24Mg") == (12, 12)
    assert parse_nucleus("18O") == (8, 10)
    with pytest.raises(ValueError, match="cannot parse"):
        parse_nucleus("24Xx")


def test_hamiltonian_lookup():
    spec = load_hamiltonian_spec("usdb", "sd")
    assert spec.valence_space == "sd"
    assert spec.files[0].path.name == "usdb.int"
    assert find_interaction_file("usdb.int") == spec.files[0].path
