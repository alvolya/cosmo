from __future__ import annotations

import pytest

from cosmo_gpu.cli import main


def test_transitions_command(capsys: pytest.CaptureFixture[str]):
    assert main(["transitions", "20Ne", "sd", "usdb", "-n", "3"]) == 0
    out = capsys.readouterr().out
    assert "0+(1) -> 2+(1)" in out
    assert "30.3598" in out


def test_sf_command(capsys: pytest.CaptureFixture[str]):
    assert main(["sf", "18O", "19O", "sd", "usdb", "-n", "2"]) == 0
    out = capsys.readouterr().out
    assert "neutron added" in out
    assert "0+(1) -> 5/2+(1)   0d5" in out


def test_strength_command(capsys: pytest.CaptureFixture[str]):
    assert main(["strength", "20Ne", "sd", "usdb", "--steps", "40"]) == 0
    out = capsys.readouterr().out
    assert "B(E2) to J = 2" in out
    assert "30.35986" in out  # converged first peak equals B(E2; 0+ -> 2+_1)
