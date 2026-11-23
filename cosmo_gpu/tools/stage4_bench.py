"""Stage-4 benchmark on one GPU: kernel time per layout and block size, end-to-end solves.

For each case: CPU table build, upload and device memory; the kernel alone (CUDA events) for
the "rows" (nvec, D) and "interleaved" (D, nvec) layouts at 1, 4, 8 and 16 vectors; and, for
the solve cases, a 10-state block-Lanczos run on the device (block 4, tol 1e-8) timed end to
end from the interaction to J. Writes ``tools/stage4_bench_<gpu>.json``.

usage (from cosmo_gpu/, on a GPU node; see stage4_bench.sb):
    UV_PROJECT_ENVIRONMENT=.venv-gpu uv run --no-sync python tools/stage4_bench.py
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import cupy as cp
import numpy as np

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import compute_levels, load_interaction
from cosmo_gpu.distribution import block_costs
from cosmo_gpu.onthefly_cuda import CudaOperator
from cosmo_gpu.operator_tables import hamiltonian_tables
from cosmo_gpu.truncation import make_truncation

ROOT = Path(__file__).resolve().parents[1]
SNT = ROOT / "tests" / "data" / "kshell" / "snt"
NVECS = (1, 4, 8, 16)
LAYOUTS = ("rows", "interleaved")


@dataclass(frozen=True)
class Bench:
    """One benchmark case."""

    name: str
    nucleus: str
    snt: str
    parity: int
    hw: tuple[int, int] | None
    solve: bool


CASES = [
    Bench("24Mg_sd", "24Mg", "usdb_from_int.snt", 1, None, False),
    Bench("28Si_sd", "28Si", "usdb_from_int.snt", 1, None, True),
    Bench("48Cr_kb3g", "48Cr", "kb3g_from_int.snt", 1, None, True),
    Bench("24Mg_1hw-", "24Mg", "fsu9_A24.snt", -1, (1, 1), True),
    Bench("28Si_1hw-", "28Si", "fsu9_A28.snt", -1, (1, 1), True),
]


def kernel_ms(op: CudaOperator, nvec: int, reps: int = 5) -> float:
    """Median kernel time (ms) of one ``apply`` on random device vectors."""
    d = op.dim
    shape = (d, nvec) if op.layout == "interleaved" else (nvec, d)
    x = cp.random.default_rng(1).standard_normal(shape)
    y = cp.empty(shape)
    op.apply(x, y)  # warm-up and compile
    times = []
    for _ in range(reps):
        a, b = cp.cuda.Event(), cp.cuda.Event()
        a.record()
        op.apply(x, y)
        b.record()
        b.synchronize()
        times.append(cp.cuda.get_elapsed_time(a, b))
    return float(np.median(times))


def bench(b: Bench) -> dict[str, object]:
    """All measurements for one case."""
    space, zv, nv, inter = load_interaction(b.nucleus, None, None, SNT / b.snt)
    trunc = make_truncation(space, zv, nv, b.hw) if b.hw else None
    basis = build_basis(space, zv, nv, parity=b.parity, truncation=trunc)
    t0 = time.perf_counter()
    tables = hamiltonian_tables(basis, inter)
    t_tables = time.perf_counter() - t0
    rec: dict[str, object] = {"dimension": basis.dim, "tables_s": round(t_tables, 3)}
    madds = float(block_costs(tables).sum())
    for layout in LAYOUTS:
        op = CudaOperator(tables, layout)
        rec[f"upload_s_{layout}"] = round(op.upload_seconds, 3)
        rec["device_table_mb"] = round(op.device_bytes / 1e6, 1)
        for nvec in NVECS:
            rec[f"kernel_ms_{layout}_{nvec}"] = round(kernel_ms(op, nvec), 3)
        del op
        cp.get_default_memory_pool().free_all_blocks()
    rec["madds_per_vector"] = madds
    k4 = {lay: float(rec[f"kernel_ms_{lay}_4"]) for lay in LAYOUTS}  # type: ignore[arg-type]
    best = min(k4, key=lambda lay: k4[lay])
    rec["best_layout_block4"] = best
    rec["gmadds_per_s_block4"] = round(4 * madds / (k4[best] * 1e-3) / 1e9, 1)
    if b.solve:
        cp.get_default_memory_pool().free_all_blocks()
        t0 = time.perf_counter()
        lev = compute_levels(
            b.nucleus,
            None,
            None,
            nev=10,
            parity=b.parity,
            hw=b.hw,
            tol=1e-8,
            snt=SNT / b.snt,
            matvec="gpu",
        )
        rec |= {
            "solve_total_s": round(time.perf_counter() - t0, 2),
            "solve_timings": {k: round(v, 3) for k, v in lev.timings.items()},
            "solve_matvecs": lev.n_matvec,
            "solve_converged": lev.converged,
            "energies": lev.energies.tolist(),
            "peak_device_mb": round(cp.get_default_memory_pool().total_bytes() / 1e6, 1),
        }
    return rec


def main(names: list[str]) -> int:
    """Run the cases (all by default) and write the JSON for this GPU."""
    props = cp.cuda.runtime.getDeviceProperties(0)
    gpu = props["name"].decode()
    tag = gpu.split()[-1].split("-")[0].lower()
    out: dict[str, object] = {"gpu": gpu, "memory_gb": round(props["totalGlobalMem"] / 1e9, 1)}
    for b in CASES:
        if names and b.name not in names:
            continue
        out[b.name] = bench(b)
        print(b.name, json.dumps(out[b.name])[:600], flush=True)
        (ROOT / "tools" / f"stage4_bench_{tag}.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
