import numpy as np
import pytest

from cosmo_gpu.basis import (
    build_basis,
    build_species_basis,
    build_species_basis_full,
    species_partitions,
)
from cosmo_gpu.modelspace import load_model_space, parse_nucleus
from cosmo_gpu.truncation import (
    make_truncation,
    min_quanta,
    orbit_quanta,
    parse_hw,
    parse_rejection,
)


@pytest.fixture(scope="module")
def psdpf():
    return load_model_space("spsdpf")


def test_quanta(psdpf):
    np.testing.assert_array_equal(orbit_quanta(psdpf), [0, 1, 1, 2, 2, 2, 3, 3, 3, 3])
    assert min_quanta(psdpf, 2) == 0
    assert min_quanta(psdpf, 8) == 6  # 0s^2 0p^6
    assert min_quanta(psdpf, 10) == 10  # + 2 in sd


def test_parse():
    assert parse_hw("1:1") == (1, 1)
    assert parse_hw("2") == (0, 2)
    sd = load_model_space("sd")
    r = parse_rejection("1s1,0d3:pn:0:2", sd)
    assert r.proton_orbits == (True, False, True)
    assert r.neutron_orbits == (True, False, True)
    assert parse_rejection("0d5:p:1:3", sd).neutron_orbits == (False, False, False)
    with pytest.raises(ValueError, match="not in sd"):
        parse_rejection("0f7:pn:0:1", sd)


@pytest.mark.parametrize("n", range(13))
def test_partition_first_equals_full_enumeration(sd_space, n):
    a = build_species_basis(sd_space, n)
    b = build_species_basis_full(sd_space, n)
    for f in (
        "dets",
        "sub_start",
        "sub_m2",
        "sub_partition",
        "sorted_dets",
        "det_sub",
        "det_local",
    ):
        np.testing.assert_array_equal(getattr(a, f), getattr(b, f), err_msg=f)


def test_partition_quanta_cap(psdpf):
    parts = species_partitions(psdpf, 10, max_quanta=11)
    q = orbit_quanta(psdpf)
    assert all(int(np.dot(p, q)) <= 11 for p in parts)
    assert all(sum(p) == 10 for p in parts)


# Dimensions from CoSMo's XsysmbsOMP (bench-psdpf/*_1hw-.xml and tests/data/cosmo_fsu9)
@pytest.mark.parametrize(
    ("nucleus", "parity", "hw", "dim"),
    [
        ("20Ne", -1, (1, 1), 26_600),
        ("22Na", -1, (1, 1), 218_924),
        ("24Mg", -1, (1, 1), 978_948),
        ("28Si", -1, (1, 1), 3_834_452),
        ("33Al", 1, (0, 0), 119),
        ("33Al", -1, (1, 1), 23_842),
    ],
)
def test_psdpf_dimensions_match_cosmo(psdpf, nucleus, parity, hw, dim):
    Z, N = parse_nucleus(nucleus)
    b = build_basis(psdpf, Z, N, parity=parity, truncation=make_truncation(psdpf, Z, N, hw))
    assert b.dim == dim


def test_block_predicate_matches_vectorized(psdpf):
    Z, N = parse_nucleus("20Ne")
    tr = make_truncation(psdpf, Z, N, (1, 1), (parse_rejection("0f7,0f5,1p3,1p1:pn:0:1", psdpf),))
    b = build_basis(psdpf, Z, N, parity=-1, truncation=tr)
    for B in range(b.n_blocks):
        assert tr.keeps(
            b.protons.sub_partition[b.block_psub[B]], b.neutrons.sub_partition[b.block_nsub[B]]
        )


def test_wrong_parity_is_empty(psdpf):
    Z, N = parse_nucleus("20Ne")
    b = build_basis(psdpf, Z, N, parity=-1, truncation=make_truncation(psdpf, Z, N, (0, 0)))
    assert b.dim == 0


def test_occupation_rejection_matches_cosmo(sd_space):
    """24Mg USDB with at most 2 nucleons in 1s1 + 0d3: CoSMo XsysmbsOMP gives 2829 states."""
    Z, N = parse_nucleus("24Mg")
    zv, nv = Z - 8, N - 8
    rule = parse_rejection("1s1,0d3:pn:0:2", sd_space)
    b = build_basis(sd_space, zv, nv, truncation=make_truncation(sd_space, zv, nv, None, (rule,)))
    assert b.dim == 2829
