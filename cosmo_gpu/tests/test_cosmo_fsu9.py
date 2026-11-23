"""psdpf (CoSMo ``spsdpf``) with FSU9 and Nħω truncation, against CoSMo DEXHVJTcsb.

References: ``data/cosmo_fsu9_reference.json`` (``data/extract_cosmo.py``). CoSMo prints 6
decimals and reads interaction files in float32, so inputs are rounded the same way here.
FSU9 includes the Lawson centre-of-mass term (``spsdpfcm.int`` scaled by 10) as a regular file.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from cosmo_gpu.cli import compute_levels

REF = json.loads((Path(__file__).parent / "data" / "cosmo_fsu9_reference.json").read_text())
E_TOL = 6e-7
CONVERGED = 1e-6  # CoSMo eigenvector error below which a state is compared


def _levels(nucleus, parity, hw, nev=10):
    return compute_levels(
        nucleus, "spsdpf", "fsu9", nev=nev, parity=parity, hw=hw, tol=1e-9, cosmo_float32=True
    )


def _compare(lev, ref, min_states):
    assert lev.converged
    assert lev.basis.dim == ref["dim"]
    ok = np.array(ref["vector_error"]) < CONVERGED
    assert ok.sum() >= min_states
    n = len(ref["energies"])
    np.testing.assert_allclose(
        lev.energies[:n][ok], np.array(ref["energies"])[ok], atol=E_TOL, rtol=0
    )
    np.testing.assert_allclose(2 * lev.spins[:n][ok], np.array(ref["j2"])[ok], atol=1e-6)


@pytest.mark.parametrize(
    ("key", "nucleus", "parity", "hw", "min_states"),
    [
        ("33Al_0hw+", "33Al", 1, (0, 0), 10),
        ("33Al_1hw-", "33Al", -1, (1, 1), 8),
        ("20Ne_1hw-", "20Ne", -1, (1, 1), 10),
        pytest.param("22Na_1hw-", "22Na", -1, (1, 1), 2, marks=pytest.mark.slow),
    ],
)
def test_levels_match_cosmo(key, nucleus, parity, hw, min_states):
    _compare(_levels(nucleus, parity, hw), REF[key], min_states)


def test_33al_merged_spectrum():
    """CoSMo Xlevels on both parities: 5/2+ ground state -303.079, 1/2- at 1.947 MeV."""
    pos = _levels("33Al", 1, (0, 0), nev=1)
    neg = _levels("33Al", -1, (1, 1), nev=1)
    assert pos.energies[0] == pytest.approx(-303.079, abs=5e-4)
    assert round(2 * pos.spins[0]) == 5
    assert round(2 * neg.spins[0]) == 1
    assert neg.energies[0] - pos.energies[0] == pytest.approx(1.947, abs=5e-4)


def test_mixed_parity_requires_parity():
    with pytest.raises(ValueError, match="mixes parities"):
        compute_levels("20Ne", "spsdpf", "fsu9", hw=(1, 1))
