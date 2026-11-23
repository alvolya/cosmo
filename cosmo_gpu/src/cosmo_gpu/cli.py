"""Command line: ``cosmo-gpu levels 24Mg sd usdb -n 10``, a ``.snt`` file or a CoSMo XML."""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
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
from cosmo_gpu.onthefly import (
    OnTheFlyOperator,
    hamiltonian_operator,
    j_expectation_otf,
    jplus_operator,
)
from cosmo_gpu.operators import (
    MatrixCache,
    State,
    Tensor,
    e2,
    m1,
    spectroscopic_factor,
    tensor_matrix,
    transition_strength,
)
from cosmo_gpu.reference_csr import hamiltonian_csr, j_expectation, jplus_csr
from cosmo_gpu.snt import read_snt, read_snt_space
from cosmo_gpu.strength import lanczos_strength
from cosmo_gpu.truncation import OccupationRule, make_truncation, parse_hw, parse_rejection

DENSE_LIMIT = 2000
OTF_LIMIT = 200_000  # 'auto' uses the on-the-fly matvec above this dimension
GPU_LAYOUT = "rows"  # vector layout of the GPU kernel (see onthefly_cuda)
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
    vectors: NDArray[np.float64] | None = None
    H: sp.csr_matrix | OnTheFlyOperator | Any | None = None  # Any: CudaOperator (needs CuPy)

    def matvec(self, X: NDArray[np.float64]) -> NDArray[np.float64]:
        """``X H^T`` for a block of row vectors, with whichever H was used."""
        if self.H is None:
            msg = "no Hamiltonian stored"
            raise ValueError(msg)
        if isinstance(self.H, OnTheFlyOperator):
            return self.H.matvec(X)
        if not isinstance(self.H, sp.csr_matrix):  # CudaOperator: device round trip
            import cupy as cp

            return np.asarray(cp.asnumpy(self.H.matvec(cp.asarray(np.atleast_2d(X)))))
        return np.ascontiguousarray((self.H @ np.atleast_2d(X).T).T)

    def states(self) -> list[State]:
        """Return the eigenvectors as ``operators.State`` objects.

        Raises
        ------
        ValueError
            If the solver did not converge or a level's <J> is not a half-integer (within
            ``BAD_SPIN_TOL``): reduced matrix elements need good J.
        """
        if self.vectors is None:
            msg = "no eigenvectors stored"
            raise ValueError(msg)
        if not self.converged:
            msg = "not all requested states converged; raise --nev or check the residuals"
            raise ValueError(msg)
        j2 = np.rint(2 * self.spins).astype(int)
        bad = np.flatnonzero(np.abs(2 * self.spins - j2) > 2 * BAD_SPIN_TOL)
        if bad.size:
            js = ", ".join(f"level {k + 1}: <J> = {self.spins[k]:.4f}" for k in bad)
            msg = f"levels without good J ({js}); near-degenerate states of different J?"
            raise ValueError(msg)
        return [State(self.vectors[k], self.basis, int(j2[k])) for k in range(len(j2))]


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
        space = read_snt_space(snt)
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


def _levels_gpu(
    basis: MBasis,
    inter: MSchemeInteraction,
    nev: int,
    block: int,
    tol: float,
    t: dict[str, float],
) -> Levels:
    """Tables on the GPU, block Lanczos on the device; eigenvectors returned to the host."""
    import cupy as cp

    from cosmo_gpu.backend_cupy import CupyBackend
    from cosmo_gpu.onthefly_cuda import (
        hamiltonian_cuda,
        j_expectation_cuda,
        jplus_cuda,
    )

    t0 = time.perf_counter()
    H = hamiltonian_cuda(basis, inter, GPU_LAYOUT)
    jop = jplus_cuda(basis, GPU_LAYOUT)
    cp.cuda.Device().synchronize()
    t["build_tables"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    res = block_lanczos(
        H.matvec, basis.dim, min(nev, basis.dim), CupyBackend(), block=block, tol=tol
    )
    spins = j_expectation_cuda(res.eigenvectors, jop, basis.m2)
    vectors = np.asarray(cp.asnumpy(res.eigenvectors), dtype=np.float64)
    t["solve"] = time.perf_counter() - t0
    return Levels(
        basis, res.eigenvalues, spins, res.residuals, res.converged, res.n_matvec, t, vectors, H
    )


Solution = tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], bool, int]


def _diagonalize(
    H: sp.csr_matrix | OnTheFlyOperator, dim: int, nev: int, block: int, tol: float
) -> Solution:
    """Lowest ``nev`` pairs: dense for a small stored H, block Lanczos otherwise."""
    if isinstance(H, sp.csr_matrix) and dim <= DENSE_LIMIT:
        w, v = sla.eigh(H.toarray())
        return w[:nev], np.ascontiguousarray(v[:, :nev].T), np.zeros(nev), True, 0
    if isinstance(H, OnTheFlyOperator):
        mv = H.matvec
    else:

        def mv(X: NDArray[np.float64]) -> NDArray[np.float64]:
            return np.ascontiguousarray((H @ X.T).T)

    res = block_lanczos(mv, dim, nev, NumpyBackend(), block=block, tol=tol)
    return res.eigenvalues, res.eigenvectors, res.residuals, res.converged, res.n_matvec


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
    matvec: str = "auto",
) -> Levels:
    """Build the basis and H, then diagonalize (dense below ``DENSE_LIMIT``).

    ``matvec`` is "stored" (CSR reference path), "otf" (operator tables and the on-the-fly
    kernel; H is never stored), "gpu" (the same tables on the GPU, Lanczos on the device; needs
    the ``gpu`` extra) or "auto" (on the fly above ``OTF_LIMIT`` states).

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

    if matvec not in ("auto", "stored", "otf", "gpu"):
        msg = f"matvec must be 'auto', 'stored', 'otf' or 'gpu', not {matvec!r}"
        raise ValueError(msg)
    if matvec == "gpu" and basis.dim > DENSE_LIMIT:
        return _levels_gpu(basis, inter, nev, block, tol, t)
    otf = basis.dim > DENSE_LIMIT and (
        matvec == "otf" or (matvec == "auto" and basis.dim > OTF_LIMIT)
    )
    t0 = time.perf_counter()
    H: sp.csr_matrix | OnTheFlyOperator
    if otf:
        H = hamiltonian_operator(basis, inter)
        jop = jplus_operator(basis)
        t["build_tables"] = time.perf_counter() - t0
    else:
        H = hamiltonian_csr(basis, inter)
        jplus, _ = jplus_csr(basis)
        t["build_H"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    energies, vectors, residuals, converged, n_matvec = _diagonalize(
        H, basis.dim, min(nev, basis.dim), block, tol
    )
    t["solve"] = time.perf_counter() - t0
    spins = (
        j_expectation_otf(vectors, jop, basis.m2)
        if otf
        else j_expectation(vectors, jplus, basis.m2)
    )
    return Levels(basis, energies, spins, residuals, converged, n_matvec, t, vectors, H)


def _format_j(j: float) -> str:
    j2 = round(2 * j)
    if abs(2 * j - j2) > 2 * BAD_SPIN_TOL:
        return f"{j:.3f}?"
    return str(j2 // 2) if j2 % 2 == 0 else f"{j2}/2"


def _system_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("space", nargs="?", help="CoSMo valence space, e.g. sd")
    p.add_argument("hamiltonian", nargs="?", help="CoSMo Hamiltonian name, e.g. usdb")
    p.add_argument("--snt", type=Path, default=None, help="KSHELL .snt file instead")
    p.add_argument("-n", "--nev", type=int, default=10)
    p.add_argument("--parity", type=int, choices=(-1, 1), default=None)
    p.add_argument(
        "--hw",
        type=parse_hw,
        default=None,
        metavar="MIN:MAX",
        help="oscillator quanta above the lowest configuration, e.g. 1:1",
    )
    p.add_argument(
        "--rejection",
        action="append",
        default=[],
        metavar="ORBS:TYPE:NMIN:NMAX",
        help="occupation window, e.g. 0f5,1p3,1p1:pn:0:4 (repeatable)",
    )
    p.add_argument(
        "--cosmo-float32",
        action="store_true",
        help="round interaction-file numbers to single precision, as CoSMo does",
    )


def _operator_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--op", choices=("E2", "M1"), default="E2")
    p.add_argument(
        "--charges", default="1.5,0.5", help="E2 effective charges e_p,e_n (default 1.5,0.5)"
    )
    p.add_argument(
        "--g",
        default="1,5.5857,0,-3.82608",
        help="M1 g factors gl_p,gs_p,gl_n,gs_n (default free nucleons)",
    )
    p.add_argument("--b", type=float, default=1.0, help="oscillator length in fm (E2 units)")


def _levels_for(args: argparse.Namespace, nucleus: str, parity: int | None = None) -> Levels:
    return compute_levels(
        nucleus,
        args.space,
        args.hamiltonian,
        args.nev,
        args.parity if parity is None else parity,
        hw=args.hw,
        rejections=tuple(args.rejection),
        snt=args.snt,
        cosmo_float32=args.cosmo_float32,
    )


def _operator(args: argparse.Namespace, lev: Levels) -> tuple[Tensor, float, str]:
    space = lev.basis.space
    if args.op == "E2":
        ep, en = (float(x) for x in args.charges.split(","))
        unit = "e^2 b^4" if args.b == 1.0 else "e^2 fm^4"
        return e2(space, ep, en), args.b**4, unit
    glp, gsp, gln, gsn = (float(x) for x in args.g.split(","))
    return m1(space, glp, gsp, gln, gsn), 1.0, "mu_N^2"


def _label(s: State, k: int, states: list[State], parity: int) -> str:
    idx = sum(1 for q in states[: k + 1] if q.j2 == s.j2)
    js = str(s.j2 // 2) if s.j2 % 2 == 0 else f"{s.j2}/2"
    return f"{js}{'+' if parity > 0 else '-'}({idx})"


def _states_or_report(lev: Levels) -> list[State] | None:
    """Return the levels as states, or None after printing why they cannot be used."""
    try:
        return lev.states()
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return None


def _cmd_transitions(args: argparse.Namespace) -> int:
    lev = _levels_for(args, args.nucleus)
    op, scale, unit = _operator(args, lev)
    states = _states_or_report(lev)
    if states is None:
        return 1
    cache = MatrixCache()
    print(f"# {args.nucleus} B({args.op}) in {unit}, i -> f among the lowest {len(states)} levels")
    for a, si in enumerate(states):
        for b, sf_ in enumerate(states):
            bval = transition_strength(op, si, sf_, cache) * scale
            if bval > 0.0:
                li = _label(si, a, states, lev.basis.parity)
                lf = _label(sf_, b, states, lev.basis.parity)
                print(f"{li:>8} -> {lf:<8} {bval:14.6f}")
    return 0


def _final_parities(args: argparse.Namespace, li: Levels) -> list[int]:
    """Return the final parities: ``--final-parity``, both if mixed, else the initial."""
    if args.final_parity is not None:
        return [args.final_parity]
    if len(set(li.basis.space.orbit_parity.tolist())) > 1:
        return [1, -1]
    return [li.basis.parity]


def _cmd_sf(args: argparse.Namespace) -> int:
    li = _levels_for(args, args.nucleus)
    si = _states_or_report(li)
    if si is None:
        return 1
    space = li.basis.space
    finals: list[tuple[Levels, list[State]]] = []
    for par in _final_parities(args, li):
        try:
            lf = _levels_for(args, args.final, par)
        except ValueError as err:  # empty basis for this parity
            print(f"# final parity {par:+d}: {err}")
            continue
        sfs = _states_or_report(lf)
        if sfs is None:
            return 1
        finals.append((lf, sfs))
    if not finals:
        print("error: no final states", file=sys.stderr)
        return 1
    lf0 = finals[0][0]
    dz = lf0.basis.protons.n_particles - li.basis.protons.n_particles
    dn = lf0.basis.neutrons.n_particles - li.basis.neutrons.n_particles
    if (dz, dn) not in ((1, 0), (0, 1)):
        print("final nucleus must have one proton or one neutron more", file=sys.stderr)
        return 2
    proton = dz == 1
    print(f"# S = |<f||a+_j||i>|^2/(2J_f+1), {'proton' if proton else 'neutron'} added")
    caches = [MatrixCache() for _ in range(space.n_orbits)]
    for a, i in enumerate(si):
        for lf, sfs in finals:
            for b, f in enumerate(sfs):
                for orb in range(space.n_orbits):
                    val = spectroscopic_factor(space, orb, proton, i, f, caches[orb])
                    if val > 1e-12:
                        print(
                            f"{_label(i, a, si, li.basis.parity):>8} -> "
                            f"{_label(f, b, sfs, lf.basis.parity):<8} "
                            f"{space.orbits[orb].name:>4} {val:10.6f}"
                        )
    return 0


def _cmd_strength(args: argparse.Namespace) -> int:
    lev = _levels_for(args, args.nucleus)
    op, scale, unit = _operator(args, lev)
    states = _states_or_report(lev)
    if states is None:
        return 1
    g = states[0]
    pivot = tensor_matrix(op, 0, lev.basis, lev.basis) @ g.vector
    st = lanczos_strength(
        lev.matvec,
        np.ascontiguousarray(pivot[None, :]),
        args.steps,
        NumpyBackend(),
    )
    if g.j2 == 0:
        # one final J (= k); by rotational invariance B = (2k + 1) |<f 0|T_0|i 0>|^2
        scale *= op.k2 + 1
        note = f"B({args.op}) to J = {op.k2 // 2}"
    else:
        note = "|<f M|T_0|i M>|^2 at the run's M, summed over final J"
    print(f"# {args.nucleus} {args.op} strength from the ground state: {note}, {unit}")
    print(f"# total {st.norm2 * scale:.6f} over {st.n_steps} Lanczos steps")
    print(f"# {'E (MeV)':>12} {'Ex':>10} {'weight':>14}")
    for e, w in zip(st.energies, st.weights, strict=True):
        print(f"  {e:12.4f} {e - lev.energies[0]:10.4f} {w * scale:14.6f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point of the ``cosmo-gpu`` command."""
    parser = argparse.ArgumentParser(prog="cosmo-gpu")
    sub = parser.add_subparsers(dest="command", required=True)
    lv = sub.add_parser("levels", help="lowest levels of a nucleus")
    lv.add_argument("nucleus", help="e.g. 24Mg, or a CoSMo system file ending in .xml")
    _system_args(lv)
    lv.add_argument("--system", default=None, help="system name in a multi-system CoSMo XML")
    lv.add_argument("--m2", type=int, default=None, help="doubled total M (default 0 or 1)")
    lv.add_argument("--block", type=int, default=4)
    lv.add_argument(
        "--matvec",
        choices=("auto", "stored", "otf", "gpu"),
        default="auto",
        help="stored CSR H, on the fly from tables (CPU or gpu), or auto (otf above 2e5 states)",
    )
    lv.add_argument("--tol", type=float, default=1e-8)

    tr = sub.add_parser("transitions", help="B(E2) or B(M1) among the lowest levels")
    tr.add_argument("nucleus")
    _system_args(tr)
    _operator_args(tr)

    sf = sub.add_parser("sf", help="spectroscopic factors for adding one nucleon")
    sf.add_argument("nucleus", help="initial nucleus (A)")
    sf.add_argument("final", help="final nucleus (A + 1)")
    _system_args(sf)
    sf.add_argument(
        "--final-parity",
        type=int,
        choices=(-1, 1),
        default=None,
        help="parity of the A+1 states (default: both in mixed-parity spaces, else --parity)",
    )

    stp = sub.add_parser("strength", help="Lanczos strength function from the ground state")
    stp.add_argument("nucleus")
    _system_args(stp)
    _operator_args(stp)
    stp.add_argument("--steps", type=int, default=100)
    args = parser.parse_args(argv)
    if args.command == "transitions":
        return _cmd_transitions(args)
    if args.command == "sf":
        return _cmd_sf(args)
    if args.command == "strength":
        return _cmd_strength(args)
    return _cmd_levels(args, parser)


def _cmd_levels(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
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
        args.matvec,
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
