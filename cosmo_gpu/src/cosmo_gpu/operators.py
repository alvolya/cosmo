"""One-body spherical tensor operators, reduced matrix elements and transition strengths.

Conventions
-----------
* Single-particle states couple ``j = l + s`` (l first) and harmonic-oscillator radial functions
  are positive at the origin, as in KSHELL (``sp_matrix_element.f90``, ``harmonic_oscillator.f90``).
  Single-particle m-scheme matrices are built directly from Clebsch-Gordan sums over
  ``|l m_l> |s m_s>``, so no reduced-matrix-element phase formula is involved.
* Wigner-Eckart (Edmonds): ``<J_f M_f|T^k_q|J_i M_i> = (-1)^(J_f - M_f) (J_f k J_i; -M_f q M_i)
  <f||T^k||i>``. ``B(k; i -> f) = |<f||T^k||i>|^2 / (2 J_i + 1)``.
* Spectroscopic factor for adding a nucleon in orbit ``j``:
  ``S = |<f||a+_j||i>|^2 / (2 J_f + 1)``.
* Lengths are in units of the oscillator length ``b`` (E2 in ``e^2 b^4``), magnetic moments in
  ``mu_N``. All angular momenta are doubled (``j2``, ``m2``, ``k2``, ``q2``).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from functools import cache
from math import pi, sqrt

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy import integrate
from scipy.special import eval_genlaguerre, gammaln

from cosmo_gpu.basis import MBasis, SpeciesBasis, build_mbasis
from cosmo_gpu.interaction import clebsch_gordan
from cosmo_gpu.modelspace import ModelSpace
from cosmo_gpu.reference_csr import (
    ANNIHILATE_N,
    ANNIHILATE_P,
    CREATE_N,
    CREATE_P,
    CX_NP,
    CX_PN,
    OperatorTerms,
    TransferTerms,
    build_operator,
    jplus_csr,
    transfer_operator,
)
from cosmo_gpu.truncation import make_truncation

THREEJ_MIN = 1e-8
ZERO = 1e-14


# ------------------------------------------------------------------ angular algebra


def threej(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol with doubled arguments."""
    if m1 + m2 + m3 != 0:
        return 0.0
    phase = -1.0 if ((j1 - j2 - m3) // 2) % 2 else 1.0
    return phase / sqrt(j3 + 1) * clebsch_gordan(j1, m1, j2, m2, j3, -m3)


@cache
def ho_radial(n1: int, l1: int, n2: int, l2: int, power: int) -> float:
    """``<n1 l1| r^power |n2 l2>`` for HO radial functions with b = 1, positive at the origin."""

    def norm(n: int, l: int) -> float:
        return sqrt(2.0 * np.exp(gammaln(n + 1) - gammaln(n + l + 1.5)))

    def rad(r: float, n: int, l: int) -> float:
        lag = float(eval_genlaguerre(n, l + 0.5, r * r))
        return float(norm(n, l) * r**l * lag * np.exp(-r * r / 2))

    val, _ = integrate.quad(
        lambda r: rad(r, n1, l1) * rad(r, n2, l2) * r ** (power + 2),
        0.0,
        np.inf,
        epsabs=1e-14,
        epsrel=1e-13,
        limit=200,
    )
    return float(val)


def gaunt(l1: int, m1: int, k: int, q: int, l2: int, m2: int) -> float:
    """``<l1 m1|Y_kq|l2 m2>`` (integer arguments, not doubled)."""
    if (l1 + k + l2) % 2:
        return 0.0
    phase = -1.0 if m1 % 2 else 1.0
    pref = sqrt((2 * l1 + 1) * (2 * k + 1) * (2 * l2 + 1) / (4 * pi))
    return (
        phase
        * pref
        * threej(2 * l1, 2 * k, 2 * l2, 0, 0, 0)
        * threej(2 * l1, 2 * k, 2 * l2, -2 * m1, 2 * q, 2 * m2)
    )


def _ls_matrix(
    space: ModelSpace,
    orbital: Callable[[int, int, int, int], float],
    spin: Callable[[int, int], float],
    radial: Callable[[int, int], float],
) -> NDArray[np.float64]:
    """m-scheme matrix ``<a m_a| R * O_l x O_s |b m_b>`` from the l and s parts.

    ``orbital(l_a, ml_a, l_b, ml_b)`` and ``spin(ms2_a, ms2_b)`` act in the uncoupled basis (``ml``
    integer, ``ms2`` doubled); ``radial(a, b)`` takes orbit positions.
    """
    nm = space.n_morbits
    orb = space.orbits
    M = np.zeros((nm, nm))
    for x in range(nm):
        a = int(space.morb_orbit[x])
        ja, ma = orb[a].j2, int(space.morb_m2[x])
        for y in range(nm):
            b = int(space.morb_orbit[y])
            r = radial(a, b)
            if r == 0.0:
                continue
            jb, mb = orb[b].j2, int(space.morb_m2[y])
            la, lb = orb[a].l, orb[b].l
            tot = 0.0
            for msa in (-1, 1):
                mla2 = ma - msa
                if abs(mla2) > 2 * la:
                    continue
                cga = clebsch_gordan(2 * la, mla2, 1, msa, ja, ma)
                if cga == 0.0:
                    continue
                for msb in (-1, 1):
                    mlb2 = mb - msb
                    if abs(mlb2) > 2 * lb:
                        continue
                    s_me = spin(msa, msb)
                    if s_me == 0.0:
                        continue
                    o_me = orbital(la, mla2 // 2, lb, mlb2 // 2)
                    if o_me == 0.0:
                        continue
                    tot += cga * clebsch_gordan(2 * lb, mlb2, 1, msb, jb, mb) * o_me * s_me
            M[x, y] = r * tot
    return M


def _ang_component(j: float, m_bra: float, m_ket: float, q: int) -> float:
    """Spherical component ``q`` of an angular momentum operator: ``<j m_bra|J_q|j m_ket>``."""
    if q == 0:
        return m_ket if m_bra == m_ket else 0.0
    if m_bra != m_ket + q:
        return 0.0
    ladder = sqrt(j * (j + 1) - m_ket * (m_ket + q))
    return -q * ladder / sqrt(2.0)  # J_{+1} = -J_+/sqrt2, J_{-1} = +J_-/sqrt2


# ------------------------------------------------------------------ operator definitions


@dataclass(frozen=True)
class Tensor:
    """A one-body spherical tensor of rank ``k2/2`` and parity ``parity``.

    ``matrices(q2)`` returns the single-particle m-scheme matrices of component ``q2`` as
    ``(proton part, neutron part)`` for particle-conserving operators, or the single matrix
    ``M[c, a]`` of ``a+_c a_a`` across species for charge exchange (``kind`` "cx_pn": proton
    created, neutron annihilated; "cx_np" the reverse).
    """

    name: str
    k2: int
    parity: int
    kind: str  # "conserving", "cx_pn" or "cx_np"
    matrices: Callable[[int], tuple[NDArray[np.float64], NDArray[np.float64]]]


def e2(space: ModelSpace, e_p: float, e_n: float) -> Tensor:
    """Electric quadrupole ``e r^2 Y_2`` with effective charges, lengths in units of b."""

    def mats(q2: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        q = q2 // 2
        orb = space.orbits

        def radial(a: int, b: int) -> float:
            if (orb[a].l + orb[b].l) % 2:
                return 0.0
            return ho_radial(orb[a].n, orb[a].l, orb[b].n, orb[b].l, 2)

        base = _ls_matrix(
            space,
            lambda la, ma, lb, mb: gaunt(la, ma, 2, q, lb, mb),
            lambda sa, sb: 1.0 if sa == sb else 0.0,
            radial,
        )
        return e_p * base, e_n * base

    return Tensor("E2", 4, 1, "conserving", mats)


def m1(space: ModelSpace, gl_p: float, gs_p: float, gl_n: float, gs_n: float) -> Tensor:
    """Magnetic dipole ``sqrt(3/4pi) (g_l l + g_s s)`` in mu_N."""

    def mats(q2: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        q = q2 // 2
        orb = space.orbits

        def same(a: int, b: int) -> float:
            return 1.0 if (orb[a].n, orb[a].l) == (orb[b].n, orb[b].l) else 0.0

        lmat = _ls_matrix(
            space,
            lambda la, ma, lb, mb: _ang_component(la, ma, mb, q) if la == lb else 0.0,
            lambda sa, sb: 1.0 if sa == sb else 0.0,
            same,
        )
        smat = _ls_matrix(
            space,
            lambda la, ma, lb, mb: 1.0 if (la == lb and ma == mb) else 0.0,
            lambda sa, sb: _ang_component(0.5, sa / 2, sb / 2, q),
            same,
        )
        f = sqrt(3 / (4 * pi))
        return f * (gl_p * lmat + gs_p * smat), f * (gl_n * lmat + gs_n * smat)

    return Tensor("M1", 2, 1, "conserving", mats)


def gamow_teller(space: ModelSpace, direction: str) -> Tensor:
    """``sigma tau``: "n->p" (beta-minus, ``a+_p a_n``) or "p->n" (beta-plus, ``a+_n a_p``)."""
    if direction not in ("n->p", "p->n"):
        msg = "direction must be 'n->p' or 'p->n'"
        raise ValueError(msg)

    def mats(q2: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        q = q2 // 2
        orb = space.orbits
        sigma = 2.0 * _ls_matrix(
            space,
            lambda la, ma, lb, mb: 1.0 if (la == lb and ma == mb) else 0.0,
            lambda sa, sb: _ang_component(0.5, sa / 2, sb / 2, q),
            lambda a, b: 1.0 if (orb[a].n, orb[a].l) == (orb[b].n, orb[b].l) else 0.0,
        )
        return sigma, sigma

    return Tensor("GT", 2, 1, "cx_pn" if direction == "n->p" else "cx_np", mats)


def number(space: ModelSpace, orbit: int, proton: bool) -> Tensor:
    """Occupation of orbit ``orbit`` (position in ``space.orbits``) for one species."""

    def mats(q2: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        d = (space.morb_orbit == orbit).astype(np.float64)
        z = np.zeros(space.n_morbits)
        return (np.diag(d), np.diag(z)) if proton else (np.diag(z), np.diag(d))

    return Tensor(f"N{orbit}{'p' if proton else 'n'}", 0, 1, "conserving", mats)


def coupled_density(space: ModelSpace, a: int, b: int, k2: int, proton: bool) -> Tensor:
    """Return the unit tensor ``[a+_a ~a_b]^k`` with ``~a_(b m) = (-1)^(j_b + m) a_(b, -m)``."""
    orb = space.orbits
    ja, jb = orb[a].j2, orb[b].j2
    parity = orb[a].parity * orb[b].parity

    def mats(q2: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        M = np.zeros((space.n_morbits, space.n_morbits))
        for x in np.flatnonzero(space.morb_orbit == a):
            ma = int(space.morb_m2[x])
            for y in np.flatnonzero(space.morb_orbit == b):
                m_annih = int(space.morb_m2[y])  # a_(b, m') with m' = -mb
                mb = -m_annih
                phase = -1.0 if ((jb + mb) // 2) % 2 else 1.0
                M[x, y] = clebsch_gordan(ja, ma, jb, mb, k2, q2) * phase
        Z = np.zeros_like(M)
        return (M, Z) if proton else (Z, M)

    return Tensor(f"[{a}x{b}]^{k2 / 2:g}", k2, parity, "conserving", mats)


# ------------------------------------------------------------------ many-body matrices


@dataclass(frozen=True)
class State:
    """An eigenvector with its basis and doubled angular momentum."""

    vector: NDArray[np.float64]
    basis: MBasis
    j2: int


def _basis_at(
    template: MBasis,
    m2: int,
    parity: int,
    protons: SpeciesBasis | None = None,
    neutrons: SpeciesBasis | None = None,
) -> MBasis:
    """Basis like ``template`` at ``m2``/``parity``, optionally with other species bases.

    A truncation is rebuilt for the destination (Z, N): its minimum quanta depend on the
    particle numbers, so the source's cannot be reused across charge exchange or transfer.
    """
    P = protons or template.protons
    N = neutrons or template.neutrons
    trunc = template.truncation
    if trunc is not None:
        trunc = make_truncation(
            template.space, P.n_particles, N.n_particles, trunc.hw, trunc.occupation
        )
    return build_mbasis(template.space, P, N, m2, parity, trunc)


def tensor_matrix(op: Tensor, q2: int, dst: MBasis, src: MBasis) -> sp.csr_matrix:
    """Return the stored matrix of component ``q2`` of ``op`` from ``src`` to ``dst``."""
    mp, mn = op.matrices(q2)
    if op.kind == "conserving":
        rows_p = np.argwhere(np.abs(mp) > ZERO).astype(np.int32)
        rows_n = np.argwhere(np.abs(mn) > ZERO).astype(np.int32)
        terms = OperatorTerms(
            ob_p=rows_p,
            ob_p_v=mp[rows_p[:, 0], rows_p[:, 1]],
            ob_n=rows_n,
            ob_n_v=mn[rows_n[:, 0], rows_n[:, 1]],
        )
        return build_operator(dst, src, terms)
    kind = CX_PN if op.kind == "cx_pn" else CX_NP
    idx = np.argwhere(np.abs(mp) > ZERO).astype(np.int32)
    return transfer_operator(
        dst, src, TransferTerms(kind, idx[:, 0], idx[:, 1], mp[idx[:, 0], idx[:, 1]])
    )


def creation_matrix(
    space: ModelSpace, orbit: int, m2: int, proton: bool, dst: MBasis, src: MBasis
) -> sp.csr_matrix:
    """``a+_(orbit, m2)`` for one species from ``src`` to ``dst``."""
    k = int(np.flatnonzero((space.morb_orbit == orbit) & (space.morb_m2 == m2))[0])
    terms = TransferTerms(
        CREATE_P if proton else CREATE_N, np.array([k], np.int32), np.zeros(1, np.int32), np.ones(1)
    )
    return transfer_operator(dst, src, terms)


def annihilation_matrix(
    space: ModelSpace, orbit: int, m2: int, proton: bool, dst: MBasis, src: MBasis
) -> sp.csr_matrix:
    """``a_(orbit, m2)`` for one species from ``src`` to ``dst``."""
    k = int(np.flatnonzero((space.morb_orbit == orbit) & (space.morb_m2 == m2))[0])
    terms = TransferTerms(
        ANNIHILATE_P if proton else ANNIHILATE_N,
        np.zeros(1, np.int32),
        np.array([k], np.int32),
        np.ones(1),
    )
    return transfer_operator(dst, src, terms)


def ladder_up(state: State, cache: MatrixCache | None = None) -> State:
    """``|J M+1> = J+ |J M> / sqrt(J(J+1) - M(M+1))`` (J+ reused from ``cache`` if given)."""
    jp, up = cache.jplus(state.basis) if cache is not None else jplus_csr(state.basis)
    j2, m2 = state.j2, state.basis.m2
    norm = sqrt((j2 * (j2 + 2) - m2 * (m2 + 2)) / 4.0)
    if norm == 0.0:
        msg = f"cannot raise M = {m2 / 2:g} for J = {j2 / 2:g}"
        raise ValueError(msg)
    return State(np.asarray(jp @ state.vector) / norm, up, j2)


def _raise(state: State, steps: int, cache: MatrixCache | None = None) -> State:
    for _ in range(steps):
        state = ladder_up(state, cache)
    return state


def _choose_components(k2: int, i: State, f: State) -> tuple[int, int, int, float]:
    """Pick ladder steps (i, f) and q2 with the largest 3j symbol; return (di, df, q2, 3j)."""
    best: tuple[int, int, int, float] | None = None
    for total in range(5):
        for di in range(total + 1):
            df = total - di
            mi, mf = i.basis.m2 + 2 * di, f.basis.m2 + 2 * df
            q2 = mf - mi
            if abs(q2) > k2 or abs(mi) > i.j2 or abs(mf) > f.j2 or (k2 - q2) % 2:
                continue
            w = threej(f.j2, k2, i.j2, -mf, q2, mi)
            if abs(w) > THREEJ_MIN and (best is None or abs(w) > abs(best[3])):
                best = (di, df, q2, w)
        if best is not None:
            return best
    msg = "no non-zero Wigner-Eckart component found (selection rule forbids it?)"
    raise ValueError(msg)


BasisKey = tuple[int, int, int, int, int, int]


def _basis_key(b: MBasis) -> BasisKey:
    return (b.m2, b.parity, id(b.protons), id(b.neutrons), id(b.truncation), b.dim)


class MatrixCache:
    """Stored operator components and J+ matrices, reused across reduced matrix elements.

    Keys identify a basis by its M, parity, species bases, truncation and dimension (bases
    rebuilt from the same parts are identical). Each entry keeps its bases alive, so the
    ``id`` values in its key cannot be reused by other objects. Use one cache per operator.
    """

    def __init__(self) -> None:
        self._ops: dict[tuple[int, ...], tuple[MBasis, MBasis, sp.csr_matrix]] = {}
        self._jp: dict[BasisKey, tuple[MBasis, sp.csr_matrix, MBasis]] = {}

    def operator(
        self,
        q2: int,
        dst: MBasis,
        src: MBasis,
        build: Callable[[int, MBasis, MBasis], sp.csr_matrix],
    ) -> sp.csr_matrix:
        """Component ``q2`` from ``src`` to ``dst``, built once."""
        key = (q2, *_basis_key(dst), *_basis_key(src))
        hit = self._ops.get(key)
        if hit is None:
            hit = (dst, src, build(q2, dst, src))
            self._ops[key] = hit
        return hit[2]

    def jplus(self, basis: MBasis) -> tuple[sp.csr_matrix, MBasis]:
        """J+ from ``basis`` and the 2M + 2 basis, built once per basis."""
        key = _basis_key(basis)
        hit = self._jp.get(key)
        if hit is None:
            jp, up = jplus_csr(basis)
            hit = (basis, jp, up)
            self._jp[key] = hit
        return hit[1], hit[2]


def reduced_matrix_element(
    i: State,
    f: State,
    component: Callable[[int, MBasis, MBasis], sp.csr_matrix],
    k2: int,
    cache: MatrixCache | None = None,
) -> float:
    """``<f||T^k||i>`` from one m-scheme matrix element, using ladders when needed.

    ``component(q2, dst, src)`` returns the stored matrix of ``T^k_q`` between two bases;
    ``cache`` (optional, one per operator) keeps those matrices and the J+ ladders.
    """
    di, df, q2, w = _choose_components(k2, i, f)
    ii, ff = _raise(i, di, cache), _raise(f, df, cache)
    if cache is None:
        mat = component(q2, ff.basis, ii.basis)
    else:
        mat = cache.operator(q2, ff.basis, ii.basis, component)
    me = float(ff.vector @ (mat @ ii.vector))
    phase = -1.0 if ((f.j2 - ff.basis.m2) // 2) % 2 else 1.0
    return me / (phase * w)


def tensor_rme(op: Tensor, i: State, f: State, cache: MatrixCache | None = None) -> float:
    """``<f||op||i>`` for a one-body tensor (conserving or charge exchange)."""
    return reduced_matrix_element(i, f, lambda q2, d, s: tensor_matrix(op, q2, d, s), op.k2, cache)


def transition_strength(op: Tensor, i: State, f: State, cache: MatrixCache | None = None) -> float:
    """``B(op; i -> f) = |<f||op||i>|^2 / (2 J_i + 1)``; 0 if a selection rule forbids it."""
    if abs(f.j2 - i.j2) > op.k2 or f.j2 + i.j2 < op.k2:
        return 0.0
    if f.basis.parity * i.basis.parity != op.parity:
        return 0.0
    return tensor_rme(op, i, f, cache) ** 2 / (i.j2 + 1)


def spectroscopic_factor(
    space: ModelSpace,
    orbit: int,
    proton: bool,
    i: State,
    f: State,
    cache: MatrixCache | None = None,
) -> float:
    """``S = |<f||a+_j||i>|^2 / (2 J_f + 1)`` for adding one nucleon in ``orbit``."""
    j2 = space.orbits[orbit].j2
    if abs(f.j2 - i.j2) > j2 or f.j2 + i.j2 < j2:
        return 0.0
    if f.basis.parity != i.basis.parity * space.orbits[orbit].parity:
        return 0.0

    def comp(q2: int, dst: MBasis, src: MBasis) -> sp.csr_matrix:
        return creation_matrix(space, orbit, q2, proton, dst, src)

    rme = reduced_matrix_element(i, f, comp, j2, cache)
    return rme**2 / (f.j2 + 1)


def total_strength(
    op: Tensor, i: State, protons: SpeciesBasis | None = None, neutrons: SpeciesBasis | None = None
) -> float:
    """``sum_f B(op; i -> f) = sum_q ||T_q psi_i||^2`` over the destination space.

    For charge exchange pass the destination species bases (``protons``, ``neutrons``).
    """
    tot = 0.0
    parity = i.basis.parity * op.parity
    for q2 in range(-op.k2, op.k2 + 1, 2):
        dst = _basis_at(i.basis, i.basis.m2 + q2, parity, protons, neutrons)
        if dst.dim == 0:
            continue
        v = tensor_matrix(op, q2, dst, i.basis) @ i.vector
        tot += float(v @ v)
    return tot


def expectation(op: Tensor, s: State) -> float:
    """``<psi|T^0_0|psi>`` for a rank-0 conserving operator (e.g. an occupation number)."""
    if op.k2 != 0 or op.kind != "conserving":
        msg = "expectation needs a rank-0 particle-conserving operator"
        raise ValueError(msg)
    return float(s.vector @ (tensor_matrix(op, 0, s.basis, s.basis) @ s.vector))


def sp_reduced(space: ModelSpace, op: Tensor, a: int, b: int, proton: bool = True) -> float:
    """Single-particle ``<a||t||b>`` of a conserving tensor from its m-scheme matrices."""
    orb = space.orbits
    ja, jb = orb[a].j2, orb[b].j2
    for q2 in range(-op.k2, op.k2 + 1, 2):
        mp, mn = op.matrices(q2)
        M = mp if proton else mn
        for x in np.flatnonzero(space.morb_orbit == a):
            ma = int(space.morb_m2[x])
            for y in np.flatnonzero(space.morb_orbit == b):
                mb = int(space.morb_m2[y])
                w = threej(ja, op.k2, jb, -ma, q2, mb)
                if abs(w) > THREEJ_MIN:
                    phase = -1.0 if ((ja - ma) // 2) % 2 else 1.0
                    return float(M[x, y]) / (phase * w)
    return 0.0
