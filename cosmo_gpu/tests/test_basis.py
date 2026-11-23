from math import comb

import numpy as np

from cosmo_gpu.basis import build_mbasis, build_species_basis
from cosmo_gpu.bits_kernels import find_det, popcount


def test_species_counts(sd_space):
    for n in range(13):
        sb = build_species_basis(sd_space, n)
        assert sb.dets.shape[0] == comb(12, n)
        assert int(sb.sub_dim.sum()) == comb(12, n)
        assert all(popcount(d) == n for d in sb.dets[:50])


def test_subblocks_are_homogeneous(sd_space):
    sb = build_species_basis(sd_space, 4)
    for s in range(sb.n_sub):
        dets = sb.dets[sb.sub_start[s] : sb.sub_start[s + 1]]
        assert np.all(np.diff(dets.astype(np.float64)) > 0)
        for d in dets:
            bits = [(int(d) >> k) & 1 for k in range(12)]
            assert sum(b * m for b, m in zip(bits, sd_space.morb_m2, strict=True)) == sb.sub_m2[s]
            occ = np.bincount(sd_space.morb_orbit, weights=bits, minlength=3)
            np.testing.assert_array_equal(occ, sb.sub_partition[s])


def test_lookup_roundtrip(sd_space):
    sb = build_species_basis(sd_space, 5)
    for i in range(0, sb.dets.shape[0], 37):
        s, loc = find_det(sb.dets[i], sb.sorted_dets, sb.det_sub, sb.det_local)
        assert sb.sub_start[s] + loc == i


def test_mbasis_dimensions(sd_space):
    P, N = build_species_basis(sd_space, 4), build_species_basis(sd_space, 4)
    dims = {}
    for m2 in range(-24, 25, 2):
        b = build_mbasis(sd_space, P, N, m2)
        dims[m2] = b.dim
        assert b.block_offset[-1] == b.dim
    assert dims[0] == 28503  # 24Mg USDB M=0, as in CoSMo
    assert sum(dims.values()) == comb(12, 4) ** 2
    assert dims[2] == dims[-2]
