"""Kernel time per variant (block of 4 vectors) on one GPU; writes variant_bench_<gpu>.json.

usage (on a GPU node): uv run --no-sync python tools/gpu/variant_bench.py
"""

from __future__ import annotations

import json
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
    "24Mg_sd": ("24Mg", "usdb_from_int.snt", 1, None),
    "28Si_sd": ("28Si", "usdb_from_int.snt", 1, None),
    "48Cr_kb3g": ("48Cr", "kb3g_from_int.snt", 1, None),
    "24Mg_1hw-": ("24Mg", "fsu9_A24.snt", -1, (1, 1)),
    "28Si_1hw-": ("28Si", "fsu9_A28.snt", -1, (1, 1)),
}
VARIANTS: dict[str, dict[str, object]] = {
    "baseline": {"lean": False, "packed": False, "threads": 128, "min_blocks": 1},
    "lean": {"lean": True, "packed": False, "threads": 128, "min_blocks": 1},
    "lean_packed": {"lean": True, "packed": True, "threads": 128, "min_blocks": 1},
    **{
        f"lean_packed_t{t}_mb{m}": {"threads": t, "min_blocks": m}
        for t, m in ((128, 6), (128, 8), (256, 3), (256, 4))
    },
    "default": {},
}


def time_ms(op: CudaOperator, nvec: int = 4, reps: int = 7) -> float:
    """Median CUDA-event time of one full matvec."""
    x = cp.asarray(np.random.default_rng(1).standard_normal((nvec, op.dim)))
    y = cp.empty_like(x)
    op.apply(x, y)
    out = []
    for _ in range(reps):
        a, b = cp.cuda.Event(), cp.cuda.Event()
        a.record()
        op.apply(x, y)
        b.record()
        b.synchronize()
        out.append(cp.cuda.get_elapsed_time(a, b))
    return float(np.median(out))


def main() -> None:
    """Time every variant on every case and save the table."""
    gpu = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    res: dict[str, object] = {"gpu": gpu}
    for name, (nucleus, snt, parity, hw) in CASES.items():
        space, zv, nv, inter = load_interaction(nucleus, None, None, SNT / snt)
        trunc = make_truncation(space, zv, nv, hw) if hw else None
        basis = build_basis(space, zv, nv, parity=parity, truncation=trunc)
        tables = hamiltonian_tables(basis, inter)
        row = {}
        for vname, kw in VARIANTS.items():
            op = CudaOperator(tables, **kw)  # type: ignore[arg-type]
            row[vname] = round(time_ms(op), 3)
            del op
            cp.get_default_memory_pool().free_all_blocks()
        res[name] = row
        print(name, json.dumps(row), flush=True)
    tag = gpu.split()[-1].split("-")[0].lower()
    (ROOT / "tools" / f"variant_bench_{tag}.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
