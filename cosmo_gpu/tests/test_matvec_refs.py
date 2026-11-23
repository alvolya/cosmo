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


@cache
def _operator(name: str) -> Any:
    from cosmo_gpu.onthefly import hamiltonian_operator

    basis, inter, _ = _case(name)
    return hamiltonian_operator(basis, inter)


@cache
def _scale(name: str) -> np.ndarray:
    basis, inter, ref = _case(name)
    H = hamiltonian_csr(basis, inter)
    return np.asarray((abs(H) @ np.abs(ref["X"]).T).T)


@pytest.mark.parametrize("name", CASES)
def test_kernel_one_call(name: str):
    _, _, ref = _case(name)
    Y = _operator(name).matvec(ref["X"])
    assert rel_error(Y, ref["Y"], _scale(name)) < TOL


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("n_src_chunks", [2, 3])
def test_kernel_source_chunks(name: str, n_src_chunks: int):
    """A matvec split by source block range: first call overwrites, the rest accumulate."""
    basis, _, ref = _case(name)
    op = _operator(name)
    X = ref["X"]
    cuts = np.linspace(0, basis.n_blocks, n_src_chunks + 1).astype(int)
    Y = np.full_like(ref["Y"], np.nan)
    for k in range(n_src_chunks):
        s0, s1 = int(cuts[k]), int(cuts[k + 1])
        lo, hi = basis.block_offset[s0], basis.block_offset[s1]
        op.apply(np.ascontiguousarray(X[:, lo:hi]), Y, src_range=(s0, s1), accumulate=k > 0)
    assert rel_error(Y, ref["Y"], _scale(name)) < TOL


@pytest.mark.parametrize("name", CASES)
def test_kernel_destination_range_leaves_the_rest(name: str):
    basis, _, ref = _case(name)
    op = _operator(name)
    d0, d1 = basis.n_blocks // 3, 2 * basis.n_blocks // 3
    Y = np.full_like(ref["Y"], 7.0)
    op.apply(ref["X"], Y, dst_range=(d0, d1))
    lo, hi = basis.block_offset[d0], basis.block_offset[d1]
    assert rel_error(Y[:, lo:hi], ref["Y"][:, lo:hi], _scale(name)[:, lo:hi]) < TOL
    assert np.all(Y[:, :lo] == 7.0)
    assert np.all(Y[:, hi:] == 7.0)
