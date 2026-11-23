from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import compute_levels, main
from cosmo_gpu.cosmo_xml import read_cosmo_xml, system_names
from cosmo_gpu.truncation import make_truncation

DATA = Path(__file__).parent / "data"
XMLS = sorted([*DATA.glob("cosmo_usdb/*.xml"), *DATA.glob("cosmo_fsu9/*.xml")])
SIZED = [p for p in XMLS if read_cosmo_xml(p).cosmo_dimension is not None]


@pytest.mark.parametrize("path", SIZED, ids=lambda p: p.stem)
def test_dimension_matches_cosmo_file(path: Path):
    s = read_cosmo_xml(path)
    trunc = (
        make_truncation(s.space, s.z_val, s.n_val, s.hw, s.occupation)
        if (s.hw or s.occupation)
        else None
    )
    basis = build_basis(s.space, s.z_val, s.n_val, s.m2, s.parity, trunc)
    assert basis.dim == s.cosmo_dimension


def test_reads_system_attributes():
    s = read_cosmo_xml(DATA / "cosmo_fsu9" / "33Al_1hw-.xml")
    assert (s.z_val, s.n_val, s.m2, s.parity, s.hw) == (13, 20, 1, -1, (1, 1))
    assert s.space.name == "spsdpf"
    assert s.hamiltonian.name == "fsu9"
    assert len(s.hamiltonian.files) == 6
    assert s.hamiltonian.files[-1].tbme_scale == 10.0
    assert s.cosmo_dimension == 23842
    assert s.cosmo_states[0].parity == -1


def test_reads_rejection():
    s = read_cosmo_xml(DATA / "cosmo_usdb" / "24Mg_rej+.xml")
    (rule,) = s.occupation
    names = [o.name for o in s.space.orbits]
    chosen = {names[i] for i, on in enumerate(rule.proton_orbits) if on}
    assert chosen == {"1s1", "0d3"}
    assert rule.neutron_orbits == rule.proton_orbits
    assert (rule.nmin, rule.nmax) == (0, 2)
    assert s.hw is None


@pytest.mark.parametrize("name", ["cosmo_usdb/18O_usdb+.xml", "cosmo_fsu9/33Al_0hw+.xml"])
def test_energies_match_states_in_file(name: str):
    s = read_cosmo_xml(DATA / name)
    lev = compute_levels("", None, None, nev=len(s.cosmo_states), system=s)
    cosmo = np.array([st.energy for st in s.cosmo_states])
    np.testing.assert_allclose(lev.energies, cosmo, atol=1.1e-3)
    np.testing.assert_array_equal(
        np.rint(2 * lev.spins).astype(int), [st.j2 for st in s.cosmo_states]
    )


def _write_variant(tmp_path: Path, src: Path, old: str, new: str) -> Path:
    out = tmp_path / src.name
    text = src.read_text(encoding="utf-8-sig")
    assert old in text
    out.write_text(text.replace(old, new, 1), encoding="utf-8")
    return out


def test_several_systems_need_a_name(tmp_path: Path):
    src = DATA / "cosmo_fsu9" / "33Al_0hw+.xml"
    text = src.read_text(encoding="utf-8-sig")
    start, end = text.index("<system "), text.index("</system>") + len("</system>")
    copy = text[start:end].replace('name="33Al_fsu9_M1T7+"', 'name="other"', 1)
    path = tmp_path / "two.xml"
    path.write_text(text[:end] + copy + text[end:], encoding="utf-8")
    assert system_names(path) == ["33Al_fsu9_M1T7+", "other"]
    with pytest.raises(ValueError, match="choose one system"):
        read_cosmo_xml(path)
    assert read_cosmo_xml(path, "other").name == "other"


def test_jj_shift_is_rejected(tmp_path: Path):
    src = DATA / "cosmo_usdb" / "24Mg_rej+.xml"
    path = _write_variant(tmp_path, src, "<rejection ", '<matrix JJ="0.5"/><rejection ')
    with pytest.raises(NotImplementedError, match="JJ"):
        read_cosmo_xml(path)


def test_system_excludes_other_sources():
    s = read_cosmo_xml(DATA / "cosmo_fsu9" / "33Al_0hw+.xml")
    with pytest.raises(ValueError, match="system file"):
        compute_levels("33Al", "spsdpf", "fsu9", system=s)


def test_cli_runs_on_xml(capsys: pytest.CaptureFixture[str]):
    assert main(["levels", str(DATA / "cosmo_fsu9" / "33Al_0hw+.xml"), "-n", "3"]) == 0
    out = capsys.readouterr().out
    assert "matches CoSMo 119" in out
    assert "CoSMo E" in out
