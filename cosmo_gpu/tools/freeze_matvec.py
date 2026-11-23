"""Freeze random-vector matvec references from the stage-1 stored H (stage-3 fixtures).

For each case it saves ``X`` (3 random vectors), ``Y = X H^T`` from ``hamiltonian_csr``, the
dimension, the seed and a hash of the block layout to ``tests/data/matvec_refs/<case>.npz``.
``X`` is stored, not regenerated, so a change in NumPy's generator cannot break the fixtures.

usage (from cosmo_gpu/): uv run python tools/freeze_matvec.py
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from cosmo_gpu.basis import MBasis, build_basis
from cosmo_gpu.cli import load_interaction
from cosmo_gpu.interaction import MSchemeInteraction
from cosmo_gpu.reference_csr import hamiltonian_csr
from cosmo_gpu.truncation import make_truncation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "data" / "matvec_refs"
SNT = ROOT / "tests" / "data" / "kshell" / "snt"
NVEC = 3

# name: (nucleus, space, hamiltonian, snt, parity, hw, seed)
CASES: dict[
    str, tuple[str, str | None, str | None, str | None, int, tuple[int, int] | None, int]
] = {
    "20Ne_sd": ("20Ne", "sd", "usdb", None, 1, None, 101),
    "33Al_1hw-": ("33Al", "spsdpf", "fsu9", None, -1, (1, 1), 102),
    "44Ti_kb3g": ("44Ti", None, None, "kb3g_from_int.snt", 1, None, 103),
}


def layout_hash(basis: MBasis) -> str:
    """SHA-256 of the block layout; changes if the basis order changes."""
    h = hashlib.sha256()
    for a in (basis.block_offset, basis.block_psub, basis.block_nsub):
        h.update(np.ascontiguousarray(a, dtype=np.int64).tobytes())
    return h.hexdigest()


def build_case(name: str) -> tuple[MBasis, MSchemeInteraction]:
    """Basis and m-scheme interaction of a reference case."""
    nucleus, space_name, ham, snt, parity, hw, _ = CASES[name]
    space, zv, nv, inter = load_interaction(nucleus, space_name, ham, SNT / snt if snt else None)
    trunc = make_truncation(space, zv, nv, hw) if hw else None
    return build_basis(space, zv, nv, parity=parity, truncation=trunc), inter


def main() -> None:
    """Write one ``.npz`` per case."""
    OUT.mkdir(parents=True, exist_ok=True)
    for name, spec in CASES.items():
        basis, inter = build_case(name)
        H = hamiltonian_csr(basis, inter)
        X = np.random.default_rng(spec[-1]).standard_normal((NVEC, basis.dim))
        Y = np.asarray((H @ X.T).T)
        np.savez_compressed(
            OUT / f"{name}.npz",
            X=X,
            Y=Y,
            dim=basis.dim,
            seed=spec[-1],
            layout=layout_hash(basis),
        )
        print(f"{name}: D = {basis.dim}, nnz(H) = {H.nnz}")


if __name__ == "__main__":
    main()
