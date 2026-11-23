"""Edge cases found in the stage 0-2 code review."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import Levels, load_interaction, main
from cosmo_gpu.interaction import read_interaction, to_mscheme
from cosmo_gpu.modelspace import (
    HamiltonianSpec,
    InteractionFile,
    load_model_space,
    parse_hamiltonian,
)
from cosmo_gpu.operators import MatrixCache, _basis_at
from cosmo_gpu.snt import read_snt
from cosmo_gpu.truncation import make_truncation, min_quanta, parse_hw

from .systems import build_system


def test_charge_exchange_destination_gets_its_own_truncation():
    """17N -> 17O in spsdpf: minimum quanta change by -1, so the window must be rebuilt."""
    space = load_model_space("spsdpf")
    hw = (1, 1)
    src = build_basis(space, 7, 10, parity=1, truncation=make_truncation(space, 7, 10, hw))
    want = build_basis(
        space, 8, 9, m2=src.m2, parity=1, truncation=make_truncation(space, 8, 9, hw)
    )
    got = _basis_at(src, src.m2, 1, want.protons, want.neutrons)
    assert got.truncation is not None
    assert got.truncation.min_quanta_p == min_quanta(space, 8)
    assert got.truncation.min_quanta_n == min_quanta(space, 9)
    assert got.dim == want.dim
    # the source window, as before the fix, keeps a different set of blocks
    wrong = build_basis(space, 8, 9, m2=src.m2, parity=1, truncation=src.truncation)
    assert wrong.dim != want.dim


def test_snt_without_core_and_with_mass_scaling(tmp_path: Path):
    snt = tmp_path / "nocore.snt"
    snt.write_text(
        """  3   3    0   0
  1  0 2 3 -1
  2  0 2 5 -1
  3  1 0 1 -1
  4  0 2 3  1
  5  0 2 5  1
  6  1 0 1  1
  6 0
  1 1 1.0
  2 2 2.0
  3 3 3.0
  4 4 1.0
  5 5 2.0
  6 6 3.0
 1 1 18 -0.3
  2 2 2 2 0 -2.0
"""
    )
    space, zv, nv, inter = load_interaction("4He", None, None, snt)
    assert (space.core_z, space.core_n, zv, nv) == (0, 0, 2, 2)
    factor = (4 / 18) ** -0.3  # A = core + valence = 4
    _, coupled = read_snt(snt, zv + nv)
    assert coupled.pp[(1, 1, 1, 1, 0)] == pytest.approx(-2.0 * factor, rel=1e-12)
    np.testing.assert_allclose(inter.pp_v, to_mscheme(coupled).pp_v, rtol=1e-14)


def _levels(spins: list[float], converged: bool = True) -> Levels:
    s = build_system("20Ne")
    n = len(spins)
    return Levels(
        basis=s.basis,
        energies=np.zeros(n),
        spins=np.array(spins),
        residuals=np.zeros(n),
        converged=converged,
        n_matvec=0,
        timings={},
        vectors=np.zeros((n, s.basis.dim)),
        H=s.H,
    )


def test_states_require_good_j_and_convergence():
    assert [st.j2 for st in _levels([0.0, 2.0]).states()] == [0, 4]
    with pytest.raises(ValueError, match="good J"):
        _levels([0.0, 2.4]).states()
    with pytest.raises(ValueError, match="converged"):
        _levels([0.0, 2.0], converged=False).states()


def test_sf_computes_both_final_parities_in_mixed_spaces(capsys: pytest.CaptureFixture[str]):
    """16O -> 17O in spsdpf at 0-1hw: 0f and 1p additions need negative-parity 17O."""
    code = main(["sf", "16O", "17O", "spsdpf", "fsu9", "--parity", "1", "--hw", "0:1", "-n", "2"])
    out = capsys.readouterr().out
    assert code == 0
    assert " 0d5 " in out  # positive parity final states
    assert "-(1)" in out  # negative-parity final states (p or pf additions)


@pytest.mark.parametrize(("text", "want"), [("2", (0, 2)), ("1:1", (1, 1)), ("0:2", (0, 2))])
def test_parse_hw_valid(text: str, want: tuple[int, int]):
    assert parse_hw(text) == want


@pytest.mark.parametrize("text", ["1:", ":2", "2:1", "-1:1"])
def test_parse_hw_rejects_bad_windows(text: str):
    with pytest.raises(ValueError, match="Nhw window"):
        parse_hw(text)


def _spec(path: Path) -> HamiltonianSpec:
    return HamiltonianSpec("t", "sd", (InteractionFile(path),))


def test_oxbash_header_continuation_and_short_lines(tmp_path: Path, sd_space):
    good = tmp_path / "cont.int"
    good.write_text("! test\n  1  1.0 2.0 -\n     3.0\n  1 1 1 1 0 1 -2.0\n")
    inter = read_interaction(_spec(good), sd_space, 2)
    np.testing.assert_array_equal(
        inter.spe[[sd_space.orbit_position(i) for i in (1, 2, 3)]], [1.0, 2.0, 3.0]
    )
    assert len(inter.tbme) == 1
    bad = tmp_path / "short.int"
    bad.write_text("  2  1.0 2.0 3.0\n  1 1 1 1 0 1 -2.0\n  1 2 1 2 1\n")
    with pytest.raises(ValueError, match="fields"):
        read_interaction(_spec(bad), sd_space, 2)


def test_scale_without_rank_is_an_error():
    h = ET.fromstring(
        '<hamiltonian name="x" valencespace="sd"><file><name>usdb.int</name>'
        "<scale>10</scale></file></hamiltonian>"
    )
    with pytest.raises(ValueError, match="rank 1 or 2"):
        parse_hamiltonian(h)


def test_cache_reuses_jplus_and_keeps_bases_alive():
    s = build_system("20Ne")
    cache = MatrixCache()
    jp1, up1 = cache.jplus(s.basis)
    jp2, up2 = cache.jplus(s.basis)
    assert jp1 is jp2
    assert up1 is up2
