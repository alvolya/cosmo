"""GPU checks; skipped unless CuPy is installed and sees a device (``uv sync --extra gpu``).

Run on GPU nodes with ``tools/gpu/gpu_tests.sb``.
"""

from __future__ import annotations

import numpy as np
import pytest

cp = pytest.importorskip("cupy")
try:
    HAS_DEVICE = cp.cuda.runtime.getDeviceCount() > 0
except Exception:  # no driver or no device
    HAS_DEVICE = False

pytestmark = pytest.mark.skipif(not HAS_DEVICE, reason="no CUDA device")

from .systems import build_system  # noqa: E402
from .test_matvec_refs import CASES as REF_CASES  # noqa: E402
from .test_matvec_refs import TOL, _case, _scale, rel_error  # noqa: E402

LAYOUTS = ["rows", "interleaved"]


def test_gather_kernel_matches_numpy():
    from cosmo_gpu.gpu_smoke import MAX_REL_ERR, run_smoke

    report = run_smoke(n_rows=50_000, nnz_per_row=16, nvec=3)
    assert report.max_rel_err <= MAX_REL_ERR
    assert report.lanczos_ok


def _cuda_op(name: str, layout: str):
    from cosmo_gpu.onthefly_cuda import hamiltonian_cuda

    basis, inter, _ = _case(name)
    return hamiltonian_cuda(basis, inter, layout)


@pytest.mark.parametrize("layout", LAYOUTS)
@pytest.mark.parametrize("name", REF_CASES)
def test_cuda_matches_frozen_references(name: str, layout: str):
    _, _, ref = _case(name)
    Y = cp.asnumpy(_cuda_op(name, layout).matvec(cp.asarray(ref["X"])))
    assert rel_error(Y, ref["Y"], _scale(name)) < TOL


@pytest.mark.parametrize("name", REF_CASES)
def test_cuda_source_chunks_and_destination_range(name: str):
    basis, _, ref = _case(name)
    op = _cuda_op(name, "rows")
    X = cp.asarray(ref["X"])
    Y = cp.full(ref["Y"].shape, np.nan)
    cuts = np.linspace(0, basis.n_blocks, 4).astype(int)
    for k in range(3):
        s0, s1 = int(cuts[k]), int(cuts[k + 1])
        lo, hi = basis.block_offset[s0], basis.block_offset[s1]
        op.apply(cp.ascontiguousarray(X[:, lo:hi]), Y, src_range=(s0, s1), accumulate=k > 0)
    assert rel_error(cp.asnumpy(Y), ref["Y"], _scale(name)) < TOL
    d0, d1 = basis.n_blocks // 3, 2 * basis.n_blocks // 3
    Y2 = cp.full(ref["Y"].shape, 7.0)
    op.apply(X, Y2, dst_range=(d0, d1))
    Y2h = cp.asnumpy(Y2)
    lo, hi = basis.block_offset[d0], basis.block_offset[d1]
    assert rel_error(Y2h[:, lo:hi], ref["Y"][:, lo:hi], _scale(name)[:, lo:hi]) < TOL
    assert np.all(Y2h[:, :lo] == 7.0)
    assert np.all(Y2h[:, hi:] == 7.0)


@pytest.mark.parametrize("layout", LAYOUTS)
@pytest.mark.parametrize(
    "spec",
    [("25Mg", "sd", "usdb", None, None), ("20Ne", "spsdpf", "fsu9", -1, (1, 1))],
    ids=["25Mg", "20Ne_1hw-"],
)
def test_cuda_matches_cpu_kernel_and_is_deterministic(spec, layout: str):
    from cosmo_gpu.onthefly import hamiltonian_operator
    from cosmo_gpu.onthefly_cuda import hamiltonian_cuda

    nucleus, space, ham, parity, hw = spec
    s = build_system(nucleus, space, ham, parity=parity, hw=hw)
    X = np.random.default_rng(11).standard_normal((4, s.basis.dim))
    cpu = hamiltonian_operator(s.basis, s.inter).matvec(X)
    op = hamiltonian_cuda(s.basis, s.inter, layout)
    g1 = op.matvec(cp.asarray(X))
    g2 = op.matvec(cp.asarray(X))
    assert bool(cp.all(g1 == g2))
    assert float(np.max(np.abs(cp.asnumpy(g1) - cpu)) / np.max(np.abs(cpu))) < 1e-13


def test_cuda_jplus_matches_stored():
    from cosmo_gpu.onthefly_cuda import jplus_cuda
    from cosmo_gpu.reference_csr import jplus_csr

    s = build_system("25Mg")
    X = np.random.default_rng(12).standard_normal((2, s.basis.dim))
    ref, _ = jplus_csr(s.basis)
    got = cp.asnumpy(jplus_cuda(s.basis).matvec(cp.asarray(X)))
    want = np.asarray((ref @ X.T).T)
    assert float(np.max(np.abs(got - want)) / np.max(np.abs(want))) < 1e-13
