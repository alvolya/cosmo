"""End-to-end GPU solves at Lanczos block sizes 4, 8 and 16 (10 states, tol 1e-8).

usage (from cosmo_gpu/, on a GPU node): uv run --no-sync python tools/gpu/block_sweep.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import cupy as cp

from cosmo_gpu.cli import compute_levels

ROOT = Path(__file__).resolve().parents[2]
SNT = ROOT / "tests" / "data" / "kshell" / "snt"
CASES = [
    ("48Cr_kb3g", "48Cr", "kb3g_from_int.snt", 1, None),
    ("24Mg_1hw-", "24Mg", "fsu9_A24.snt", -1, (1, 1)),
    ("28Si_1hw-", "28Si", "fsu9_A28.snt", -1, (1, 1)),
]


def main() -> None:
    """Print and save matvec counts and times per block size."""
    gpu = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
    out: dict[str, object] = {"gpu": gpu}
    for name, nucleus, snt, parity, hw in CASES:
        rows = {}
        for block in (4, 8, 16):
            t0 = time.perf_counter()
            lev = compute_levels(
                nucleus,
                None,
                None,
                nev=10,
                parity=parity,
                hw=hw,
                tol=1e-8,
                block=block,
                snt=SNT / snt,
                matvec="gpu",
            )
            rows[block] = {
                "total_s": round(time.perf_counter() - t0, 2),
                "solve_s": round(lev.timings["solve"], 2),
                "matvecs": lev.n_matvec,
                "converged": lev.converged,
                "e0": lev.energies[0],
            }
            print(name, block, rows[block], flush=True)
        out[name] = rows
    tag = gpu.split()[-1].split("-")[0].lower()
    (ROOT / "tools" / f"block_sweep_{tag}.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
