"""Stage-3 benchmark: on-the-fly matvec timings and a solver comparison with KSHELL.

Both codes read the same ``.snt`` files and use thick-restart block Lanczos with block 4,
``max_lanc_vec`` (our ``ncv``) 60 and ``n_restart_vec`` (our ``n_keep``) 15. KSHELL stops on the
change of the eigenvalues, we stop on the residual norm, both at 1e-8, so the matvec counts are
compared side by side with the convergence histories rather than as a pass/fail test.

usage (from cosmo_gpu/, in a Slurm job; see stage3_bench.sb):
    KSHELL_DIR=~/kshell uv run python tools/stage3_bench.py
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import numba
import numpy as np

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import load_interaction
from cosmo_gpu.lanczos import block_lanczos
from cosmo_gpu.onthefly import hamiltonian_operator
from cosmo_gpu.truncation import make_truncation

ROOT = Path(__file__).resolve().parents[1]
SNT = ROOT / "tests" / "data" / "kshell" / "snt"
OUT = ROOT / "tools"
BLOCK, NCV, KEEP, TOL, NEV = 4, 60, 15, 1e-8, 10
# CoSMo stored-H matvec (one vector, optimized build, 16 threads, AMD EPYC 9654), stage 0
COSMO_MATVEC_S = {"24Mg_1hw-": 0.073, "28Si_1hw-": 0.365}


@dataclass(frozen=True)
class Bench:
    """One benchmark case."""

    name: str
    nucleus: str
    snt: str
    z_val: int
    n_val: int
    parity: int
    hw: tuple[int, int] | None
    solve: bool  # full solver comparison (else timing only)


CASES = [
    Bench("24Mg_sd", "24Mg", "usdb_from_int.snt", 4, 4, 1, None, False),
    Bench("28Si_sd", "28Si", "usdb_from_int.snt", 6, 6, 1, None, True),
    Bench("48Cr_kb3g", "48Cr", "kb3g_from_int.snt", 4, 4, 1, None, True),
    Bench("24Mg_1hw-", "24Mg", "fsu9_A24.snt", 12, 12, -1, (1, 1), False),
    Bench("28Si_1hw-", "28Si", "fsu9_A28.snt", 14, 14, -1, (1, 1), False),
]


def _kshell_driver() -> ModuleType:
    path = ROOT / "tests" / "data" / "kshell" / "run_kshell.py"
    spec = importlib.util.spec_from_file_location("run_kshell", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses look their module up here
    spec.loader.exec_module(mod)
    return mod


TAG = os.environ.get("BENCH_TAG", "")
SOLVE = {n for n in os.environ.get("BENCH_SOLVE", "").replace(":", ",").split(",") if n}


def _solve(b: Bench) -> bool:
    return b.solve or b.name in SOLVE


def ours(b: Bench) -> dict[str, object]:
    """Table setup, matvec timing and (if ``b.solve``) a block-Lanczos run."""
    space, zv, nv, inter = load_interaction(b.nucleus, None, None, SNT / b.snt)
    assert (zv, nv) == (b.z_val, b.n_val)
    trunc = make_truncation(space, zv, nv, b.hw) if b.hw else None
    t0 = time.perf_counter()
    basis = build_basis(space, zv, nv, parity=b.parity, truncation=trunc)
    t1 = time.perf_counter()
    op = hamiltonian_operator(basis, inter)
    X = np.random.default_rng(1).standard_normal((BLOCK, basis.dim))
    op.matvec(X)  # warm-up (JIT)
    reps = []
    for _ in range(3):
        s = time.perf_counter()
        op.matvec(X)
        reps.append(time.perf_counter() - s)
    rec: dict[str, object] = {
        "dimension": basis.dim,
        "basis_s": round(t1 - t0, 3),
        "tables_s": round(op.setup_seconds, 3),
        "table_mb": round(sum(op.tables.nbytes().values()) / 1e6, 1),
        "madds_per_vector": op.cost,
        "matvec_block4_s": round(min(reps), 4),
        "gmadds_per_s": round(BLOCK * op.cost / min(reps) / 1e9, 2),
        "threads": numba.get_num_threads(),
    }
    if _solve(b):
        history: list[list[object]] = []
        t0 = time.perf_counter()
        res = block_lanczos(
            op.matvec,
            basis.dim,
            NEV,
            NumpyBackend(),
            block=BLOCK,
            ncv=NCV,
            n_keep=KEEP,
            tol=TOL,
            history=lambda n, th, r: history.append([n, th.tolist(), r.tolist()]),
        )
        rec |= {
            "solve_s": round(time.perf_counter() - t0, 2),
            "matvecs": res.n_matvec,
            "converged": res.converged,
            "energies": res.eigenvalues.tolist(),
            "history": history,
        }
    return rec


def kshell(b: Bench, drv: ModuleType, kshell_bin: Path) -> dict[str, object]:
    """KSHELL block Lanczos with the same settings (timing-only cases stop after 25 iterations)."""
    case = drv.Case(b.name, b.nucleus, b.snt, b.z_val, b.n_val, b.parity, b.hw, n_eigen=NEV)
    settings: dict[str, object] = {
        "n_block": BLOCK,
        "max_lanc_vec": NCV,
        "n_restart_vec": KEEP,
        "tol": TOL,
        "maxiter": 2000 if _solve(b) else 25,
    }
    rec: dict[str, object] = drv.run_solver(kshell_bin, case, settings, f"bench{TAG}")
    calls = rec.get("operate_calls")
    if _solve(b) and calls is not None:
        rec["matvecs"] = BLOCK * len(rec["history"])  # one operate (block of 4) per iteration
    return rec


def report(results: dict[str, dict[str, dict[str, object]]]) -> str:
    """Markdown tables for the design plan."""
    lines = [
        "| Case | D | tables (s, MB) | ours: s per matvec of 4 | Gmadd/s "
        "| KSHELL: s per matvec of 4 | CoSMo stored, 1 vector |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for name, r in results.items():
        o, k = r["ours"], r["kshell"]
        cosmo = COSMO_MATVEC_S.get(name)
        lines.append(
            f"| {name} | {o['dimension']:,} | {o['tables_s']} s, {o['table_mb']} MB | "
            f"{o['matvec_block4_s']} | {o['gmadds_per_s']} | {k.get('seconds_per_call')} | "
            f"{cosmo if cosmo is not None else '-'} |"
        )
    lines += [
        "",
        "| Case | ours: matvecs, solve s | KSHELL: matvecs, solve s | max abs dE (MeV) |",
        "| --- | --- | --- | --- |",
    ]
    for name, r in results.items():
        o, k = r["ours"], r["kshell"]
        if "matvecs" not in o:
            continue
        de = float(np.max(np.abs(np.array(o["energies"]) - np.array(k["energies"]))))
        ours_s = f"{o['matvecs']}, {o['solve_s']}"
        kshell_s = f"{k.get('matvecs')}, {k['seconds']}"
        lines.append(f"| {name} | {ours_s} | {kshell_s} | {de:.1e} |")
    return "\n".join(lines) + "\n"


def main(names: list[str]) -> int:
    """Run the benchmark and write ``stage3_bench.json`` and ``stage3_bench.md``."""
    kshell_bin = Path(os.environ.get("KSHELL_DIR", Path.home() / "kshell")) / "bin"
    drv = _kshell_driver()
    results: dict[str, dict[str, dict[str, object]]] = {}
    for b in CASES:
        if names and b.name not in names:
            continue
        results[b.name] = {"ours": ours(b), "kshell": kshell(b, drv, kshell_bin)}
        print(
            b.name,
            json.dumps({k: v for k, v in results[b.name]["ours"].items() if k != "history"}),
            flush=True,
        )
    suffix = f"_{TAG}" if TAG else ""
    (OUT / f"stage3_bench{suffix}.json").write_text(json.dumps(results, indent=1) + "\n")
    (OUT / f"stage3_bench{suffix}.md").write_text(report(results))
    print(report(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
