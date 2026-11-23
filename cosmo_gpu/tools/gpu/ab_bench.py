"""Interleaved A/B timing of the original and the default kernel (guards against clock ramp).

After a warm-up of at least WARM_S seconds of kernel launches, the two kernels alternate for
ROUNDS rounds (A B, then B A), REPS calls each; medians and per-round ratios are reported.
usage (on a GPU node): uv run --no-sync python tools/gpu/ab_bench.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import cupy as cp
import numpy as np

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import load_interaction
from cosmo_gpu.onthefly_cuda import CudaOperator
from cosmo_gpu.operator_tables import hamiltonian_tables
from cosmo_gpu.truncation import make_truncation

ROOT = Path(__file__).resolve().parents[2]
SNT = ROOT / "tests" / "data" / "kshell" / "snt"
CASES = {
    "48Cr_kb3g": ("48Cr", "kb3g_from_int.snt", 1, None),
    "24Mg_1hw-": ("24Mg", "fsu9_A24.snt", -1, (1, 1)),
    "28Si_1hw-": ("28Si", "fsu9_A28.snt", -1, (1, 1)),
}
ORIGINAL = {"lean": False, "packed": False, "threads": 128, "min_blocks": 1}
WARM_S, ROUNDS, REPS = 10.0, 10, 5


def run(op: CudaOperator, x: cp.ndarray, y: cp.ndarray, reps: int) -> list[float]:
    """Return the times (ms) of ``reps`` back-to-back calls."""
    out = []
    for _ in range(reps):
        a, b = cp.cuda.Event(), cp.cuda.Event()
        a.record()
        op.apply(x, y)
        b.record()
        b.synchronize()
        out.append(cp.cuda.get_elapsed_time(a, b))
    return out


def main() -> None:
    """Warm up, alternate the two kernels, print and save the medians."""
    gpu = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    res: dict[str, object] = {"gpu": gpu}
    for name, (nucleus, snt, parity, hw) in CASES.items():
        space, zv, nv, inter = load_interaction(nucleus, None, None, SNT / snt)
        trunc = make_truncation(space, zv, nv, hw) if hw else None
        basis = build_basis(space, zv, nv, parity=parity, truncation=trunc)
        tables = hamiltonian_tables(basis, inter)
        ops = {"original": CudaOperator(tables, **ORIGINAL), "default": CudaOperator(tables)}
        x = cp.asarray(np.random.default_rng(1).standard_normal((4, basis.dim)))
        y = cp.empty_like(x)
        t0 = time.perf_counter()
        while time.perf_counter() - t0 < WARM_S:
            for op in ops.values():
                run(op, x, y, 1)
        times: dict[str, list[float]] = {k: [] for k in ops}
        ratios = []
        for r in range(ROUNDS):
            order = list(ops) if r % 2 == 0 else list(ops)[::-1]
            med = {}
            for k in order:
                ts = run(ops[k], x, y, REPS)
                times[k] += ts
                med[k] = float(np.median(ts))
            ratios.append(med["original"] / med["default"])
        rec = {
            "original_ms": round(float(np.median(times["original"])), 3),
            "default_ms": round(float(np.median(times["default"])), 3),
            "speedup_median": round(float(np.median(ratios)), 3),
            "speedup_min_max": [round(min(ratios), 3), round(max(ratios), 3)],
        }
        res[name] = rec
        print(name, json.dumps(rec), flush=True)
        del ops
        cp.get_default_memory_pool().free_all_blocks()
    tag = gpu.split()[-1].split("-")[0].lower()
    (ROOT / "tools" / f"ab_bench_{tag}.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
