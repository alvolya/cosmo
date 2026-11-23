"""Table sizes and remote fractions of H for psdpf FSU9 1hw cases (stage-2 report).

usage (from cosmo_gpu/, in a Slurm job for 28Si):
    uv run python tools/remote_fraction.py [24Mg 28Si ...] > tools/remote_fraction.txt
"""

from __future__ import annotations

import sys
import time

from cosmo_gpu.basis import build_basis
from cosmo_gpu.distribution import block_costs, remote_fraction
from cosmo_gpu.interaction import read_interaction, to_mscheme
from cosmo_gpu.modelspace import load_hamiltonian_spec, load_model_space, parse_nucleus
from cosmo_gpu.operator_tables import hamiltonian_tables
from cosmo_gpu.truncation import make_truncation

RANKS = (2, 4, 8)


def report(nucleus: str) -> None:
    """Print table sizes and remote fractions for ``nucleus`` in spsdpf FSU9 at 1hw-."""
    space = load_model_space("spsdpf")
    z, n = parse_nucleus(nucleus)
    inter = to_mscheme(read_interaction(load_hamiltonian_spec("fsu9", "spsdpf"), space, z + n))
    trunc = make_truncation(space, z, n, (1, 1))
    t0 = time.perf_counter()
    basis = build_basis(space, z, n, parity=-1, truncation=trunc)
    t1 = time.perf_counter()
    tables = hamiltonian_tables(basis, inter)
    t2 = time.perf_counter()
    sizes = tables.nbytes()
    cost = block_costs(tables).sum()
    print(f"## {nucleus} spsdpf fsu9 1hw-: D = {basis.dim:,}, blocks = {basis.n_blocks:,}")
    print(f"species dims p = {len(basis.protons.dets):,}, n = {len(basis.neutrons.dets):,}")
    print(f"basis {t1 - t0:.1f} s, tables {t2 - t1:.1f} s, tasks = {tables.n_tasks:,}")
    print("table bytes: " + ", ".join(f"{k} {v / 1e6:.1f} MB" for k, v in sizes.items()))
    print(
        f"total tables {sum(sizes.values()) / 1e6:.1f} MB; gather multiply-adds per vector "
        f"{cost:.3e} (stored H at 12 B per nonzero would be about {12 * cost / 1e9:.1f} GB)"
    )
    print(
        "| ranks | max received states | max remote fraction | mean remote fraction | "
        "max cost share |"
    )
    print("| --- | --- | --- | --- | --- |")
    for n_ranks in RANKS:
        r = remote_fraction(tables, n_ranks)
        print(
            f"| {n_ranks} | {int(r.remote_states.max()):,} | {r.remote_fraction.max():.3f} | "
            f"{r.remote_fraction.mean():.3f} | {r.cost_share.max():.3f} |"
        )
    print(flush=True)


if __name__ == "__main__":
    for nuc in sys.argv[1:] or ["24Mg", "28Si"]:
        report(nuc)
