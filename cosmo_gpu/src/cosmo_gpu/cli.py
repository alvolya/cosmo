"""Command line: ``python -m cosmo_gpu levels 24Mg sd usdb -n 10``."""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla
from numpy.typing import NDArray

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.basis import MBasis, build_mbasis, build_species_basis
from cosmo_gpu.interaction import read_interaction, to_mscheme
from cosmo_gpu.lanczos import block_lanczos
from cosmo_gpu.modelspace import load_hamiltonian_spec, load_model_space, parse_nucleus
from cosmo_gpu.reference_csr import hamiltonian_csr, j_expectation, jplus_csr

DENSE_LIMIT = 2000
BAD_SPIN_TOL = 1e-4


@dataclass
class Levels:
    """Lowest eigenvalues of one system with their J and residuals."""

    basis: MBasis
    energies: NDArray[np.float64]
    spins: NDArray[np.float64]
    residuals: NDArray[np.float64]
    converged: bool
    n_matvec: int
    timings: dict[str, float]


def compute_levels(
    nucleus: str,
    space_name: str,
    hamiltonian: str,
    nev: int = 10,
    parity: int | None = None,
    m2: int | None = None,
    block: int = 4,
    tol: float = 1e-8,
    cosmo_float32: bool = False,
) -> Levels:
    """Build the basis and stored H, then diagonalize (dense below ``DENSE_LIMIT``)."""
    t: dict[str, float] = {}
    t0 = time.perf_counter()
    space = load_model_space(space_name)
    spec = load_hamiltonian_spec(hamiltonian, space_name)
    Z, N = parse_nucleus(nucleus)
    zv, nv = Z - space.core_z, N - space.core_n
    inter = to_mscheme(read_interaction(spec, space, zv + nv, cosmo_float32))
    t["interaction"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    P = build_species_basis(space, zv)
    Nb = build_species_basis(space, nv)
    basis = build_mbasis(space, P, Nb, (zv + nv) % 2 if m2 is None else m2, parity)
    t["basis"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    H = hamiltonian_csr(basis, inter)
    jplus, _ = jplus_csr(basis)
    t["build_H"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    nev = min(nev, basis.dim)
    energies: NDArray[np.float64]
    vectors: NDArray[np.float64]
    residuals: NDArray[np.float64]
    if basis.dim <= DENSE_LIMIT:
        w, v = sla.eigh(H.toarray())
        energies, vectors = w[:nev], np.ascontiguousarray(v[:, :nev].T)
        residuals = np.zeros(nev)
        converged, n_matvec = True, 0
    else:
        res = block_lanczos(
            lambda X: np.ascontiguousarray((H @ X.T).T),
            basis.dim,
            nev,
            NumpyBackend(),
            block=block,
            tol=tol,
        )
        energies, vectors = res.eigenvalues, res.eigenvectors
        residuals, converged, n_matvec = res.residuals, res.converged, res.n_matvec
    t["solve"] = time.perf_counter() - t0
    spins = j_expectation(vectors, jplus, basis.m2)
    return Levels(basis, energies, spins, residuals, converged, n_matvec, t)


def _format_j(j: float) -> str:
    j2 = round(2 * j)
    if abs(2 * j - j2) > 2 * BAD_SPIN_TOL:
        return f"{j:.3f}?"
    return str(j2 // 2) if j2 % 2 == 0 else f"{j2}/2"


def main(argv: list[str] | None = None) -> int:
    """Entry point of the ``cosmo-gpu`` command."""
    parser = argparse.ArgumentParser(prog="cosmo-gpu")
    sub = parser.add_subparsers(dest="command", required=True)
    lv = sub.add_parser("levels", help="lowest levels of a nucleus")
    lv.add_argument("nucleus", help="e.g. 24Mg")
    lv.add_argument("space", help="CoSMo valence space, e.g. sd")
    lv.add_argument("hamiltonian", help="CoSMo Hamiltonian name, e.g. usdb")
    lv.add_argument("-n", "--nev", type=int, default=10)
    lv.add_argument("--parity", type=int, choices=(-1, 1), default=None)
    lv.add_argument("--m2", type=int, default=None, help="doubled total M (default 0 or 1)")
    lv.add_argument("--block", type=int, default=4)
    lv.add_argument("--tol", type=float, default=1e-8)
    lv.add_argument(
        "--cosmo-float32",
        action="store_true",
        help="round interaction-file numbers to single precision, as CoSMo does",
    )
    args = parser.parse_args(argv)

    lev = compute_levels(
        args.nucleus,
        args.space,
        args.hamiltonian,
        args.nev,
        args.parity,
        args.m2,
        args.block,
        args.tol,
        args.cosmo_float32,
    )
    b = lev.basis
    sign = "+" if b.parity > 0 else "-"
    print(f"# {args.nucleus} {args.space} {args.hamiltonian}  2M={b.m2} parity={sign}")
    print(f"# dimension {b.dim}  blocks {b.n_blocks}  matvecs {lev.n_matvec}")
    print("# timings " + "  ".join(f"{k}={v:.2f}s" for k, v in lev.timings.items()))
    print(f"# {'k':>3} {'E (MeV)':>14} {'Ex':>10} {'J':>6} {'residual':>10}")
    for i, (e, j, r) in enumerate(zip(lev.energies, lev.spins, lev.residuals, strict=True)):
        ex = e - lev.energies[0]
        print(f"  {i + 1:>3} {e:14.6f} {ex:10.4f} {_format_j(j) + sign:>6} {r:10.2e}")
    if not lev.converged:
        print("# WARNING: not all requested states converged", file=sys.stderr)
        return 1
    return 0
