"""GPU checks; skipped unless CuPy is installed and sees a device (``uv sync --extra gpu``)."""

from __future__ import annotations

import pytest

cp = pytest.importorskip("cupy")
try:
    HAS_DEVICE = cp.cuda.runtime.getDeviceCount() > 0
except Exception:  # no driver or no device
    HAS_DEVICE = False

pytestmark = pytest.mark.skipif(not HAS_DEVICE, reason="no CUDA device")


def test_gather_kernel_matches_numpy():
    from cosmo_gpu.gpu_smoke import MAX_REL_ERR, run_smoke

    report = run_smoke(n_rows=50_000, nnz_per_row=16, nvec=3)
    assert report.max_rel_err <= MAX_REL_ERR
    assert report.lanczos_ok
