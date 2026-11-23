"""Frozen random-vector matvec references (``tools/freeze_matvec.py``)."""

from __future__ import annotations

import importlib.util
from functools import cache
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from cosmo_gpu.reference_csr import hamiltonian_csr

ROOT = Path(__file__).resolve().parents[1]
REFS = ROOT / "tests" / "data" / "matvec_refs"
TOL = 1e-13


@cache
def _freeze() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "freeze_matvec", ROOT / "tools" / "freeze_matvec.py"
    )
    assert spec is not None
    assert spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@cache
def _case(name: str) -> tuple[Any, Any, dict[str, Any]]:
    basis, inter = _freeze().build_case(name)
    with np.load(REFS / f"{name}.npz") as f:
        ref = {k: f[k] for k in f.files}
    return basis, inter, ref


def rel_error(Y: np.ndarray, Yref: np.ndarray, scale: np.ndarray) -> float:
    """Max elementwise error over a per-row scale (|H| |X|, so cancellations are allowed)."""
    return float(np.max(np.abs(Y - Yref) / np.maximum(scale, 1e-300)))


CASES = ["20Ne_sd", "33Al_1hw-", "44Ti_kb3g"]


@pytest.mark.parametrize("name", CASES)
def test_layout_unchanged(name: str):
    basis, _, ref = _case(name)
    assert basis.dim == int(ref["dim"])
    assert _freeze().layout_hash(basis) == str(ref["layout"])


@pytest.mark.parametrize("name", CASES)
def test_stored_h_reproduces_reference(name: str):
    basis, inter, ref = _case(name)
    H = hamiltonian_csr(basis, inter)
    Y = np.asarray((H @ ref["X"].T).T)
    scale = np.asarray((abs(H) @ np.abs(ref["X"]).T).T)
    assert rel_error(Y, ref["Y"], scale) < TOL
