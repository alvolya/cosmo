"""Launch the GPU matvec kernel for one case, for Nsight Compute / Systems.

usage (on a GPU node): uv run --no-sync python tools/gpu/profile_kernel.py CASE [NVEC] [REPS]
CASE is one of 28Si_1hw-, 24Mg_1hw-, 48Cr_kb3g. Prints the CUDA-event kernel time per call.
"""

from __future__ import annotations

import sys
from pathlib import Path

import cupy as cp
import numpy as np

from cosmo_gpu.basis import build_basis
from cosmo_gpu.cli import load_interaction
from cosmo_gpu.onthefly_cuda import CudaOperator
from cosmo_gpu.operator_tables import hamiltonian_tables
from cosmo_gpu.truncation import make_truncation

SNT = Path(__file__).resolve().parents[2] / "tests" / "data" / "kshell" / "snt"
CASES = {
    "28Si_1hw-": ("28Si", "fsu9_A28.snt", -1, (1, 1)),
    "24Mg_1hw-": ("24Mg", "fsu9_A24.snt", -1, (1, 1)),
    "48Cr_kb3g": ("48Cr", "kb3g_from_int.snt", 1, None),
}


def main() -> None:
    """Build, warm up, then time ``reps`` kernel calls."""
    name = sys.argv[1]
    nvec = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    reps = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    nucleus, snt, parity, hw = CASES[name]
    space, zv, nv, inter = load_interaction(nucleus, None, None, SNT / snt)
    trunc = make_truncation(space, zv, nv, hw) if hw else None
    basis = build_basis(space, zv, nv, parity=parity, truncation=trunc)
    op = CudaOperator(hamiltonian_tables(basis, inter))
    x = cp.asarray(np.random.default_rng(1).standard_normal((nvec, basis.dim)))
    y = cp.empty_like(x)
    op.apply(x, y)
    cp.cuda.Device().synchronize()
    for _ in range(reps):
        a, b = cp.cuda.Event(), cp.cuda.Event()
        a.record()
        op.apply(x, y)
        b.record()
        b.synchronize()
        print(f"{name} nvec={nvec} kernel {cp.cuda.get_elapsed_time(a, b):.3f} ms", flush=True)


if __name__ == "__main__":
    main()
