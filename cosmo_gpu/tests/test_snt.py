from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cosmo_gpu.cli import compute_levels
from cosmo_gpu.interaction import _exchange_phase, to_mscheme
from cosmo_gpu.snt import read_snt, write_snt

from .systems import build_system


def _assert_same_mscheme(a, b):
    np.testing.assert_allclose(a.spe_p, b.spe_p, atol=1e-12)
    np.testing.assert_allclose(a.spe_n, b.spe_n, atol=1e-12)
    for name in ("pp", "nn", "pn"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        np.testing.assert_allclose(getattr(a, name + "_v"), getattr(b, name + "_v"), atol=1e-12)


@pytest.mark.parametrize(
    ("nucleus", "space", "ham"), [("20Ne", "sd", "usdb"), ("20Ne", "spsdpf", "fsu9")]
)
def test_write_read_round_trip(tmp_path: Path, nucleus, space, ham):
    s = build_system(
        nucleus, space, ham, parity=1 if space == "sd" else -1, hw=None if space == "sd" else (1, 1)
    )
    pn = s.coupled.to_pn()
    f = tmp_path / "x.snt"
    write_snt(f, pn, comment=f"{space} {ham}")
    space2, back = read_snt(f, 4)
    assert [(o.n, o.l, o.j2) for o in space2.orbits] == [(o.n, o.l, o.j2) for o in s.space.orbits]
    assert (space2.core_z, space2.core_n) == (s.space.core_z, s.space.core_n)
    np.testing.assert_array_equal(back.spe_p, pn.spe_p)
    np.testing.assert_array_equal(back.spe_n, pn.spe_n)
    assert back.pp == pn.pp
    assert back.nn == pn.nn
    assert back.pn == pn.pn
    _assert_same_mscheme(to_mscheme(back), s.inter)


def _tiny_snt(tbme_lines: list[str], method: str = "0") -> str:
    """sd shell, two-body part only given by the caller."""
    head = """! test
  3   3    8   8
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
  4 4 1.5
  5 5 2.5
  6 6 3.5
"""
    return head + f" {len(tbme_lines)} {method}\n" + "\n".join(tbme_lines) + "\n"


def test_mass_scaling_and_species(tmp_path: Path):
    f = tmp_path / "t.snt"
    f.write_text(_tiny_snt(["1 1 1 1 0 -2.0", "4 5 4 5 2 0.5", "1 5 2 4 1 0.25"], "1 18 -0.3"))
    space, inter = read_snt(f, 8)  # A = 24
    factor = (24 / 18) ** -0.3
    np.testing.assert_array_equal(inter.spe_p, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(inter.spe_n, [1.5, 2.5, 3.5])
    assert inter.pp == {(0, 0, 0, 0, 0): pytest.approx(-2.0 * factor)}
    assert inter.nn == {(0, 1, 0, 1, 2): pytest.approx(0.5 * factor)}
    assert inter.pn == {(0, 1, 1, 0, 1): pytest.approx(0.25 * factor)}
    assert space.orbits[1].name == "0d5"


def test_neutron_first_pairs_get_exchange_phase(tmp_path: Path):
    """<d5_n d3_p|V|d3_p d5_n> must equal -(-1)^(ja+jb-J) <d3_p d5_n|V|d3_p d5_n>."""
    J = 2
    a = tmp_path / "a.snt"
    b = tmp_path / "b.snt"
    a.write_text(_tiny_snt([f"1 5 1 5 {J} 0.7"]))
    b.write_text(_tiny_snt([f"5 1 1 5 {J} 0.7"]))
    _, ia = read_snt(a, 2)
    _, ib = read_snt(b, 2)
    phase = -_exchange_phase(3, 5, J, 0)
    assert ib.pn[(0, 1, 0, 1, J)] == pytest.approx(phase * ia.pn[(0, 1, 0, 1, J)])


def test_rejects_off_diagonal_one_body(tmp_path: Path):
    f = tmp_path / "t.snt"
    f.write_text(
        _tiny_snt([]).replace("  2 2 2.0", "  2 2 2.0\n  1 3 0.1").replace("  6 0", "  7 0")
    )
    with pytest.raises(NotImplementedError, match="off-diagonal"):
        read_snt(f, 2)


def test_cli_snt_matches_cosmo_names(tmp_path: Path):
    s = build_system("20Ne")
    f = tmp_path / "usdb20.snt"
    write_snt(f, s.coupled.to_pn())
    via_names = compute_levels("20Ne", "sd", "usdb", nev=5)
    via_snt = compute_levels("20Ne", None, None, nev=5, snt=f)
    np.testing.assert_allclose(via_snt.energies, via_names.energies, atol=1e-10)
    np.testing.assert_allclose(via_snt.spins, via_names.spins, atol=1e-8)
