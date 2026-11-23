"""Stage-0b GPU smoke test: a CUDA C gather kernel and the solver on CuPy.

``python -m cosmo_gpu.gpu_smoke`` prints one JSON line with the device, the CUDA versions, the
kernel's error against NumPy, its achieved bandwidth, and a block-Lanczos run of 20Ne sd USDB
on the device compared with the CPU.

The kernel is a multi-vector CSR gather, ``y[v, i] = sum_k val[k] * x[v, src[k]]`` over the row
``i`` of a CSR matrix: FP64, gather only (one thread per output, no atomics), which is the
access pattern of the jump-table gather of the on-the-fly matvec.
"""

from __future__ import annotations

import json
import platform
import time
from dataclasses import asdict, dataclass

import cupy as cp
import cupyx.scipy.sparse as csp
import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.backend_cupy import CupyBackend
from cosmo_gpu.basis import build_basis
from cosmo_gpu.interaction import read_interaction, to_mscheme
from cosmo_gpu.lanczos import block_lanczos
from cosmo_gpu.modelspace import load_hamiltonian_spec, load_model_space
from cosmo_gpu.reference_csr import hamiltonian_csr

GATHER_SOURCE = r"""
extern "C" __global__
void csr_gather(const long long* __restrict__ row_ptr,
                const int* __restrict__ src,
                const double* __restrict__ val,
                const double* __restrict__ x,
                double* __restrict__ y,
                const long long n_rows,
                const long long n_cols,
                const int nvec)
{
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n_rows * nvec) return;
    long long i = t % n_rows;          // consecutive threads take consecutive rows
    int v = (int)(t / n_rows);
    const double* xv = x + (long long)v * n_cols;
    double acc = 0.0;
    for (long long k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
        acc += val[k] * xv[src[k]];
    }
    y[(long long)v * n_rows + i] = acc;
}
"""

MAX_REL_ERR = 1e-13
LANCZOS_TOL = 1e-10


@dataclass
class SmokeReport:
    """Result of one smoke run on one device."""

    host: str
    device: str
    compute_capability: str
    driver_version: int
    runtime_version: int
    cupy_version: str
    n_rows: int
    nnz: int
    nvec: int
    max_rel_err: float
    kernel_ms: float
    gb_per_s: float
    kernel_ok: bool
    lanczos_dim: int = 0
    lanczos_max_de: float = float("nan")
    lanczos_ok: bool = False


def gather_kernel() -> cp.RawKernel:
    """Compile (and cache) the gather kernel."""
    return cp.RawKernel(GATHER_SOURCE, "csr_gather")


def run_gather(
    mat: sp.csr_matrix, x: NDArray[np.float64], repeats: int = 20
) -> tuple[NDArray[np.float64], float]:
    """Return ``(x @ mat.T)`` computed on the device, and the mean kernel time in ms."""
    n_rows, n_cols = mat.shape
    nvec = x.shape[0]
    d_ptr = cp.asarray(mat.indptr.astype(np.int64))
    d_src = cp.asarray(mat.indices.astype(np.int32))
    d_val = cp.asarray(mat.data.astype(np.float64))
    d_x = cp.asarray(np.ascontiguousarray(x, dtype=np.float64))
    d_y = cp.zeros((nvec, n_rows), dtype=cp.float64)
    kernel = gather_kernel()
    threads = 256
    blocks = (n_rows * nvec + threads - 1) // threads
    args = (d_ptr, d_src, d_val, d_x, d_y, np.int64(n_rows), np.int64(n_cols), np.int32(nvec))
    kernel((blocks,), (threads,), args)  # warm-up and compile
    cp.cuda.Device().synchronize()
    start, stop = cp.cuda.Event(), cp.cuda.Event()
    start.record()
    for _ in range(repeats):
        kernel((blocks,), (threads,), args)
    stop.record()
    stop.synchronize()
    ms = cp.cuda.get_elapsed_time(start, stop) / repeats
    return np.asarray(cp.asnumpy(d_y), dtype=np.float64), float(ms)


def lanczos_check(nev: int = 5) -> tuple[int, float]:
    """Run block Lanczos on 20Ne sd USDB on CPU and on the device; return (dim, max |dE|)."""
    space = load_model_space("sd")
    inter = to_mscheme(read_interaction(load_hamiltonian_spec("usdb", "sd"), space, 4))
    basis = build_basis(space, 2, 2, parity=1)
    H = hamiltonian_csr(basis, inter)
    cpu = block_lanczos(
        lambda X: np.ascontiguousarray((H @ X.T).T), basis.dim, nev, NumpyBackend(), tol=1e-12
    )
    Hd = csp.csr_matrix(H)
    gpu = block_lanczos(
        lambda X: cp.ascontiguousarray((Hd @ X.T).T), basis.dim, nev, CupyBackend(), tol=1e-12
    )
    return basis.dim, float(np.max(np.abs(gpu.eigenvalues - cpu.eigenvalues)))


def run_smoke(
    n_rows: int = 2_000_000, nnz_per_row: int = 32, nvec: int = 4, seed: int = 7
) -> SmokeReport:
    """Run the gather kernel on a random CSR matrix and the device Lanczos check."""
    rng = np.random.default_rng(seed)
    n_cols = n_rows
    indptr = np.arange(0, (n_rows + 1) * nnz_per_row, nnz_per_row, dtype=np.int64)
    indices = rng.integers(0, n_cols, size=n_rows * nnz_per_row, dtype=np.int64)
    data = rng.standard_normal(n_rows * nnz_per_row)
    mat = sp.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
    x = rng.standard_normal((nvec, n_cols))
    y, ms = run_gather(mat, x)
    ref = np.asarray((mat @ x.T).T)
    scale = np.abs(mat) @ np.abs(x.T)
    err = float(np.max(np.abs(y - ref) / np.maximum(np.asarray(scale).T, 1e-300)))
    # bytes: row_ptr + src + val once, x gathered per nonzero per vector, y written
    nnz = mat.nnz
    moved = 8 * (n_rows + 1) + 12 * nnz + 8 * nvec * (nnz + n_rows)
    props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
    report = SmokeReport(
        host=platform.node(),
        device=props["name"].decode(),
        compute_capability=f"{props['major']}.{props['minor']}",
        driver_version=cp.cuda.runtime.driverGetVersion(),
        runtime_version=cp.cuda.runtime.runtimeGetVersion(),
        cupy_version=cp.__version__,
        n_rows=n_rows,
        nnz=nnz,
        nvec=nvec,
        max_rel_err=err,
        kernel_ms=round(ms, 3),
        gb_per_s=round(moved / (ms * 1e-3) / 1e9, 1),
        kernel_ok=err <= MAX_REL_ERR,
    )
    report.lanczos_dim, report.lanczos_max_de = lanczos_check()
    report.lanczos_ok = report.lanczos_max_de <= LANCZOS_TOL
    return report


def main() -> int:
    """Run the smoke test and print one JSON line."""
    t0 = time.perf_counter()
    report = run_smoke()
    out = asdict(report) | {"seconds": round(time.perf_counter() - t0, 1)}
    print(json.dumps(out))
    return 0 if report.kernel_ok and report.lanczos_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
