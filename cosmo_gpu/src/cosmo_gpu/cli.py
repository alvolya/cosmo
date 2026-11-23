"""Command line: ``cosmo-gpu levels 24Mg sd usdb -n 10``, a ``.snt`` file or a CoSMo XML."""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as sla
from numpy.typing import NDArray

from cosmo_gpu.backend import NumpyBackend
from cosmo_gpu.basis import MBasis, build_basis
from cosmo_gpu.cosmo_xml import CosmoSystem, read_cosmo_xml
from cosmo_gpu.interaction import MSchemeInteraction, read_interaction, to_mscheme
from cosmo_gpu.lanczos import block_lanczos
from cosmo_gpu.modelspace import (
    ModelSpace,
    load_hamiltonian_spec,
    load_model_space,
    parse_nucleus,
)
from cosmo_gpu.reference_csr import hamiltonian_csr, j_expectation, jplus_csr
from cosmo_gpu.snt import read_snt
from cosmo_gpu.truncation import OccupationRule, make_truncation, parse_hw, parse_rejection

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


def load_interaction(
    nucleus: str,
    space_name: str | None,
    hamiltonian: str | None,
    snt: Path | None = None,
    cosmo_float32: bool = False,
) -> tuple[ModelSpace, int, int, MSchemeInteraction]:
    """Return the space, valence (Z, N) and m-scheme interaction, from CoSMo names or a .snt."""
    Z, N = parse_nucleus(nucleus)
    if snt is not None:
        if space_name is not None or hamiltonian is not None:
            msg = "give either a CoSMo space and Hamiltonian or a .snt file, not both"
            raise ValueError(msg)
        space, _ = read_snt(snt, 0)
        zv, nv = Z - space.core_z, N - space.core_n
        space, from_core = read_snt(snt, zv + nv)
        return space, zv, nv, to_mscheme(from_core)
    if space_name is None or hamiltonian is None:
        msg = "give a CoSMo space and Hamiltonian, or a .snt file"
        raise ValueError(msg)
    space = load_model_space(space_name)
    spec = load_hamiltonian_spec(hamiltonian, space_name)
    zv, nv = Z - space.core_z, N - space.core_n
    return space, zv, nv, to_mscheme(read_interaction(spec, space, zv + nv, cosmo_float32))


def compute_levels(
    nucleus: str,
    space_name: str | None,
    hamiltonian: str | None,
    nev: int = 10,
    parity: int | None = None,
    m2: int | None = None,
    block: int = 4,
    tol: float = 1e-8,
    cosmo_float32: bool = False,
    hw: tuple[int, int] | None = None,
    rejections: tuple[str, ...] = (),
    snt: Path | None = None,
    system: CosmoSystem | None = None,
) -> Levels:
    """Build the basis and stored H, then diagonalize (dense below ``DENSE_LIMIT``).

    ``hw`` is the Nħω window above the lowest configuration and ``rejections`` are occupation
    windows in CoSMo's ``ORBS:TYPE:NMIN:NMAX`` form. ``parity`` is required in spaces that mix
    parities. With ``snt`` the space and interaction come from a KSHELL file instead of the CoSMo
    catalogue (``space_name`` and ``hamiltonian`` must then be None). With ``system`` (from
    ``read_cosmo_xml``) the space, Hamiltonian, valence numbers, truncations and, unless given,
    parity and M come from a CoSMo system file; ``nucleus`` is then ignored.
    """
    t: dict[str, float] = {}
    t0 = time.perf_counter()
    rules: tuple[OccupationRule, ...]
    if system is not None:
        if space_name or hamiltonian or snt or hw or rejections:
            msg = "a CoSMo system file sets the space, Hamiltonian and truncations itself"
            raise ValueError(msg)
        space, zv, nv = system.space, system.z_val, system.n_val
        coupled = read_interaction(system.hamiltonian, space, zv + nv, cosmo_float32)
        inter = to_mscheme(coupled)
        parity = system.parity if parity is None else parity
        m2 = system.m2 if m2 is None else m2
        hw, rules = system.hw, system.occupation
    else:
        space, zv, nv, inter = load_interaction(
            nucleus, space_name, hamiltonian, snt, cosmo_float32
        )
        rules = tuple(parse_rejection(r, space) for r in rejections)
    t["interaction"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    trunc = make_truncation(space, zv, nv, hw, rules) if (hw or rules) else None
    if parity is None and len(set(space.orbit_parity.tolist())) > 1:
        msg = f"valence space {space.name} mixes parities: give parity +1 or -1 (--parity)"
        raise ValueError(msg)
    basis = build_basis(space, zv, nv, m2, parity, trunc)
    if basis.dim == 0:
        msg = "the basis is empty: this parity/truncation admits no configuration"
        raise ValueError(msg)
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
    lv.add_argument("nucleus", help="e.g. 24Mg, or a CoSMo system file ending in .xml")
    lv.add_argument("space", nargs="?", help="CoSMo valence space, e.g. sd")
    lv.add_argument("hamiltonian", nargs="?", help="CoSMo Hamiltonian name, e.g. usdb")
    lv.add_argument("--snt", type=Path, default=None, help="KSHELL .snt file instead")
    lv.add_argument("--system", default=None, help="system name in a multi-system CoSMo XML")
    lv.add_argument("-n", "--nev", type=int, default=10)
    lv.add_argument("--parity", type=int, choices=(-1, 1), default=None)
    lv.add_argument("--m2", type=int, default=None, help="doubled total M (default 0 or 1)")
    lv.add_argument(
        "--hw",
        type=parse_hw,
        default=None,
        metavar="MIN:MAX",
        help="oscillator quanta above the lowest configuration, e.g. 1:1",
    )
    lv.add_argument(
        "--rejection",
        action="append",
        default=[],
        metavar="ORBS:TYPE:NMIN:NMAX",
        help="occupation window, e.g. 0f5,1p3,1p1:pn:0:4 (repeatable)",
    )
    lv.add_argument("--block", type=int, default=4)
    lv.add_argument("--tol", type=float, default=1e-8)
    lv.add_argument(
        "--cosmo-float32",
        action="store_true",
        help="round interaction-file numbers to single precision, as CoSMo does",
    )
    args = parser.parse_args(argv)

    system = None
    if args.nucleus.endswith(".xml"):
        system = read_cosmo_xml(Path(args.nucleus), args.system)
    elif args.system is not None:
        parser.error("--system needs a CoSMo XML file")
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
        args.hw,
        tuple(args.rejection),
        args.snt,
        system,
    )
    b = lev.basis
    sign = "+" if b.parity > 0 else "-"
    trunc = b.truncation.describe() if b.truncation else "none"
    if system is not None:
        source = f"{system.name} ({system.space.name} {system.hamiltonian.name})"
    elif args.snt:
        source = f"{args.nucleus} {args.snt}"
    else:
        source = f"{args.nucleus} {args.space} {args.hamiltonian}"
    print(f"# {source}  2M={b.m2} parity={sign}")
    print(f"# truncation {trunc}")
    dim_note = ""
    if system is not None and system.cosmo_dimension is not None:
        same = "matches" if system.cosmo_dimension == b.dim else "DIFFERS FROM"
        dim_note = f"  ({same} CoSMo {system.cosmo_dimension})"
    print(f"# dimension {b.dim}{dim_note}  blocks {b.n_blocks}  matvecs {lev.n_matvec}")
    print("# timings " + "  ".join(f"{k}={v:.2f}s" for k, v in lev.timings.items()))
    cosmo_e = [st.energy for st in system.cosmo_states] if system is not None else []
    extra = f" {'CoSMo E':>11}" if cosmo_e else ""
    print(f"# {'k':>3} {'E (MeV)':>14} {'Ex':>10} {'J':>6} {'residual':>10}{extra}")
    for i, (e, j, r) in enumerate(zip(lev.energies, lev.spins, lev.residuals, strict=True)):
        ex = e - lev.energies[0]
        col = f" {cosmo_e[i]:11.3f}" if i < len(cosmo_e) else ""
        print(f"  {i + 1:>3} {e:14.6f} {ex:10.4f} {_format_j(j) + sign:>6} {r:10.2e}{col}")
    if not lev.converged:
        print("# WARNING: not all requested states converged", file=sys.stderr)
        return 1
    return 0
