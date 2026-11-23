"""Operator tables: the factorized form of H, J+ and other operators that stage-3 kernels apply.

Layout (see ``docs/design-plan.md``, "Data layouts" and "Kernel ABI"):

* **Moves.** A one-body move ``a+_c a_a`` of one species has a class ``(d2m, dparity)``. Within
  a class, moves are numbered by a label ``q``.
* **Jump tables.** For each (destination sub-block, source sub-block) pair of one species that a
  single move connects: CSR over the destination local index; each entry holds the source local
  index (``int32``) and the label with the fermion sign in the top bit (``int16``). A sub-block
  pair fixes ``d2m`` and parity, so a table belongs to one class. All tables (both species) share
  one row-pointer array: table ``t`` has its row ``r`` at ``jump_row_ptr[jump_table_start[t] + r]``.
* **Coupling tables.** Dense ``V[q, r]`` per (proton class, neutron class) pair. A one-body
  operator acting on one species uses the identity on the other, with a single label.
* **Same-species tables.** pp (or nn) part of H plus single-particle energies, as CSR per
  (destination sub-block, source sub-block) pair of one species; table ``u`` has its row ``r``
  at ``ss_row_ptr[ss_start[u] + r]``.
* **Task list.** CSR over destination blocks; each task names its source block, proton and
  neutron jump tables (-1 = identity), coupling or same-species table, and kind
  (0 = pn product, 1 = pp with neutron spectators, 2 = nn with proton spectators).

``expand_to_csr`` turns the tables back into a stored matrix; it exists to test the tables
against ``reference_csr``, not to be fast.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray

from cosmo_gpu.basis import MBasis, SpeciesBasis, build_mbasis
from cosmo_gpu.interaction import MSchemeInteraction
from cosmo_gpu.modelspace import ModelSpace
from cosmo_gpu.operator_tables_kernels import (
    build_tasks,
    expand_tasks,
    jump_entries,
    species_entries,
)
from cosmo_gpu.reference_csr import _max_bucket, _two_body

LABEL_MAX = 0x7FFF
SIGN_BIT = 0x8000


@dataclass(frozen=True)
class MoveSet:
    """One-body moves of one species, grouped into classes ``(d2m, dparity)``."""

    c: NDArray[np.int32]
    a: NDArray[np.int32]
    cls: NDArray[np.int32]
    label: NDArray[np.int32]
    class_keys: tuple[tuple[int, int], ...]

    @property
    def n_class(self) -> int:
        """Number of classes."""
        return len(self.class_keys)

    def class_size(self, k: int) -> int:
        """Return the number of moves (labels) in class ``k``."""
        return int(np.count_nonzero(self.cls == k))


def make_moves(space: ModelSpace, pairs: NDArray[np.int32]) -> MoveSet:
    """Group distinct moves ``(c, a)`` into classes; labels follow sorted ``(c, a)`` order."""
    pairs = np.unique(np.asarray(pairs, dtype=np.int32).reshape(-1, 2), axis=0)
    m2 = space.morb_m2
    par = space.orbit_parity[space.morb_orbit]
    d2m = m2[pairs[:, 0]] - m2[pairs[:, 1]]
    dpar = par[pairs[:, 0]] * par[pairs[:, 1]]
    keys = sorted({(int(x), int(y)) for x, y in zip(d2m, dpar, strict=True)})
    index = {k: i for i, k in enumerate(keys)}
    cls = np.array([index[(int(x), int(y))] for x, y in zip(d2m, dpar, strict=True)], np.int32)
    label = np.zeros(len(pairs), dtype=np.int32)
    for k in range(len(keys)):
        sel = np.flatnonzero(cls == k)
        if len(sel) > LABEL_MAX:
            msg = f"class {keys[k]} has {len(sel)} moves, more than {LABEL_MAX} labels"
            raise ValueError(msg)
        label[sel] = np.arange(len(sel), dtype=np.int32)
    return MoveSet(pairs[:, 0].copy(), pairs[:, 1].copy(), cls, label, tuple(keys))


@dataclass(frozen=True)
class JumpTables:
    """Jump tables of one species between two species bases."""

    table_dst_sub: NDArray[np.int32]
    table_src_sub: NDArray[np.int32]
    table_class: NDArray[np.int32]
    row_start: NDArray[np.int64]  # first row of each table in row_ptr
    row_ptr: NDArray[np.int64]
    src: NDArray[np.int32]
    op: NDArray[np.int16]


def _group_by_pair(
    dst_sub: NDArray[np.int32],
    dst_local: NDArray[np.int32],
    src_sub: NDArray[np.int32],
    src_local: NDArray[np.int32],
    extra: NDArray[np.int64],
    dst_dims: NDArray[np.int64],
) -> tuple[
    NDArray[np.int64], NDArray[np.int32], NDArray[np.int32], NDArray[np.int64], NDArray[np.int64]
]:
    """Sort entries into per-(dst_sub, src_sub) CSR tables.

    Returns the sort order, each table's (dst_sub, src_sub), each table's first row and the
    shared row pointer.
    """
    order = np.lexsort((extra, src_local, dst_local, src_sub, dst_sub))
    ds, ss = dst_sub[order], src_sub[order]
    new = np.ones(len(order), dtype=bool)
    new[1:] = (ds[1:] != ds[:-1]) | (ss[1:] != ss[:-1])
    first = np.flatnonzero(new)
    t_dst = ds[first].astype(np.int32)
    t_src = ss[first].astype(np.int32)
    n_rows = dst_dims[t_dst]
    row_start = np.concatenate([[0], np.cumsum(n_rows)]).astype(np.int64)
    tid = np.cumsum(new) - 1
    g = row_start[tid] + dst_local[order]
    counts = np.bincount(g, minlength=int(row_start[-1]))
    row_ptr = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    return order, t_dst, t_src, row_start[:-1].copy(), row_ptr


def jump_tables(
    dst: SpeciesBasis, src: SpeciesBasis, space: ModelSpace, moves: MoveSet
) -> JumpTables:
    """All single-move connections from ``src`` to ``dst`` determinants of one species."""
    nm = space.n_morbits
    order_a = np.argsort(moves.a, kind="stable")
    move_ptr = np.zeros(nm + 1, dtype=np.int64)
    np.add.at(move_ptr, moves.a.astype(np.int64) + 1, 1)
    move_ptr = np.cumsum(move_ptr)
    e_src, e_dsub, e_dloc, e_cls, e_lab, e_sgn = jump_entries(
        src.dets,
        nm,
        move_ptr,
        moves.c[order_a],
        moves.cls[order_a],
        moves.label[order_a],
        dst.sorted_dets,
        dst.det_sub,
        dst.det_local,
    )
    s_sub = (np.searchsorted(src.sub_start, e_src, side="right") - 1).astype(np.int32)
    s_loc = (e_src - src.sub_start[s_sub]).astype(np.int32)
    order, t_dst, t_src, row_start, row_ptr = _group_by_pair(
        e_dsub, e_dloc, s_sub, s_loc, e_lab.astype(np.int64), dst.sub_dim
    )
    first = row_ptr[row_start]
    t_cls = e_cls[order][first] if len(order) else np.zeros(0, np.int32)
    op = e_lab[order].astype(np.int32) | np.where(e_sgn[order] < 0, SIGN_BIT, 0)
    return JumpTables(
        table_dst_sub=t_dst,
        table_src_sub=t_src,
        table_class=t_cls.astype(np.int32),
        row_start=row_start,
        row_ptr=row_ptr,
        src=s_loc[order].astype(np.int32),
        op=op.astype(np.uint16).view(np.int16),
    )


@dataclass(frozen=True)
class SpeciesTables:
    """Same-species operator of one species as CSR tables per sub-block pair."""

    table_dst_sub: NDArray[np.int32]
    table_src_sub: NDArray[np.int32]
    row_start: NDArray[np.int64]
    row_ptr: NDArray[np.int64]
    col: NDArray[np.int32]
    val: NDArray[np.float64]


def species_tables(
    sb: SpeciesBasis,
    space: ModelSpace,
    diag: NDArray[np.float64],
    tb: NDArray[np.int32],
    tb_v: NDArray[np.float64],
) -> SpeciesTables:
    """Single-particle energies ``diag`` plus like-particle two-body terms on one species."""
    nm = space.n_morbits
    ptr, cc, vv = _two_body(tb, tb_v, nm)
    z = sb.n_particles
    buf = 1 + z * (z - 1) // 2 * _max_bucket((ptr, cc, vv))
    col_ptr, rows, vals = species_entries(
        sb.dets,
        nm,
        buf,
        np.ascontiguousarray(diag, dtype=np.float64),
        ptr,
        cc,
        vv,
        sb.sorted_dets,
        sb.det_sub,
        sb.det_local,
        sb.sub_start,
    )
    cols = np.repeat(np.arange(len(sb.dets), dtype=np.int64), np.diff(col_ptr))
    d_sub = (np.searchsorted(sb.sub_start, rows, side="right") - 1).astype(np.int32)
    s_sub = (np.searchsorted(sb.sub_start, cols, side="right") - 1).astype(np.int32)
    d_loc = (rows - sb.sub_start[d_sub]).astype(np.int32)
    s_loc = (cols - sb.sub_start[s_sub]).astype(np.int32)
    order, t_dst, t_src, row_start, row_ptr = _group_by_pair(
        d_sub, d_loc, s_sub, s_loc, np.zeros(len(rows), np.int64), sb.sub_dim
    )
    return SpeciesTables(t_dst, t_src, row_start, row_ptr, s_loc[order], vals[order].copy())


@dataclass(frozen=True)
class OperatorTables:
    """All tables of one operator from basis ``src`` to basis ``dst``, flattened for kernels."""

    dst: MBasis
    src: MBasis
    # tasks, CSR over destination blocks
    task_ptr: NDArray[np.int64]
    task_src_block: NDArray[np.int32]
    task_pjump: NDArray[np.int32]
    task_njump: NDArray[np.int32]
    task_coupling: NDArray[np.int32]
    task_kind: NDArray[np.int8]
    # jump tables of both species, proton tables first
    jump_table_start: NDArray[np.int64]
    jump_row_ptr: NDArray[np.int64]
    jump_src: NDArray[np.int32]
    jump_op: NDArray[np.int16]
    # coupling tables, dense row-major
    coup_start: NDArray[np.int64]
    coup_ncol: NDArray[np.int32]
    coup_data: NDArray[np.float64]
    # same-species tables, proton tables first
    ss_start: NDArray[np.int64]
    ss_row_ptr: NDArray[np.int64]
    ss_col: NDArray[np.int32]
    ss_val: NDArray[np.float64]

    @property
    def n_tasks(self) -> int:
        """Number of tasks."""
        return int(self.task_ptr[-1])

    def nbytes(self) -> dict[str, int]:
        """Bytes per table family."""
        fam = {
            "tasks": (
                "task_ptr",
                "task_src_block",
                "task_pjump",
                "task_njump",
                "task_coupling",
                "task_kind",
            ),
            "jumps": ("jump_table_start", "jump_row_ptr", "jump_src", "jump_op"),
            "couplings": ("coup_start", "coup_ncol", "coup_data"),
            "same_species": ("ss_start", "ss_row_ptr", "ss_col", "ss_val"),
        }
        return {k: sum(int(getattr(self, a).nbytes) for a in v) for k, v in fam.items()}


def _by_dst(
    table_dst_sub: NDArray[np.int32], n_sub: int, offset: int = 0
) -> tuple[NDArray[np.int64], NDArray[np.int32]]:
    """CSR from destination sub-block to table ids (``offset`` added to the ids)."""
    order = np.argsort(table_dst_sub, kind="stable")
    ptr = np.zeros(n_sub + 1, dtype=np.int64)
    np.add.at(ptr, table_dst_sub.astype(np.int64) + 1, 1)
    return np.cumsum(ptr), (order + offset).astype(np.int32)


_EMPTY_SS = SpeciesTables(
    np.zeros(0, np.int32),
    np.zeros(0, np.int32),
    np.zeros(0, np.int64),
    np.zeros(1, np.int64),
    np.zeros(0, np.int32),
    np.zeros(0, np.float64),
)


def assemble(
    dst: MBasis,
    src: MBasis,
    jp: JumpTables | None,
    jn: JumpTables | None,
    n_pclass: int,
    n_nclass: int,
    couplings: dict[tuple[int, int], NDArray[np.float64]],
    ssp: SpeciesTables | None = None,
    ssn: SpeciesTables | None = None,
) -> OperatorTables:
    """Flatten species tables and build the task list.

    ``couplings`` maps (proton class, neutron class) to a dense ``V[q, r]``; class ``-1`` on
    either side means the identity on that species (a single label 0).
    """
    empty_j = JumpTables(
        *(np.zeros(0, np.int32),) * 3,
        np.zeros(0, np.int64),
        np.zeros(1, np.int64),
        np.zeros(0, np.int32),
        np.zeros(0, np.int16),
    )
    jp = jp or empty_j
    jn = jn or empty_j
    ssp = ssp or _EMPTY_SS
    ssn = ssn or _EMPTY_SS
    n_jp = len(jp.table_dst_sub)
    n_ssp = len(ssp.table_dst_sub)
    # jump tables: concatenate, shifting neutron row pointers and row starts
    jump_row_ptr = np.concatenate([jp.row_ptr[:-1], jn.row_ptr + jp.row_ptr[-1]])
    jump_table_start = np.concatenate([jp.row_start, jn.row_start + len(jp.row_ptr) - 1])
    jump_src = np.concatenate([jp.src, jn.src])
    jump_op = np.concatenate([jp.op, jn.op])
    ss_row_ptr = np.concatenate([ssp.row_ptr[:-1], ssn.row_ptr + ssp.row_ptr[-1]])
    ss_start = np.concatenate([ssp.row_start, ssn.row_start + len(ssp.row_ptr) - 1])
    ss_col = np.concatenate([ssp.col, ssn.col])
    ss_val = np.concatenate([ssp.val, ssn.val])
    # couplings
    lookup = np.full((n_pclass + 1, n_nclass + 1), -1, dtype=np.int32)
    coup_start = [0]
    coup_ncol: list[int] = []
    data: list[NDArray[np.float64]] = []
    for cid, ((pc, nc), V0) in enumerate(sorted(couplings.items())):
        lookup[pc if pc >= 0 else n_pclass, nc if nc >= 0 else n_nclass] = cid
        V = np.ascontiguousarray(V0, dtype=np.float64)
        coup_ncol.append(V.shape[1])
        data.append(V.ravel())
        coup_start.append(coup_start[-1] + V.size)
    P, N = dst.protons, dst.neutrons
    pj_ptr, pj_ids = _by_dst(jp.table_dst_sub, P.n_sub)
    nj_ptr, nj_ids = _by_dst(jn.table_dst_sub, N.n_sub)
    pss_ptr, pss_ids = _by_dst(ssp.table_dst_sub, P.n_sub)
    nss_ptr, nss_ids = _by_dst(ssn.table_dst_sub, N.n_sub, offset=n_ssp)
    nj_src = np.concatenate([np.zeros(n_jp, np.int32), jn.table_src_sub])  # indexed by global id
    nj_cls = np.concatenate([np.zeros(n_jp, np.int32), jn.table_class])
    nj_ids = (nj_ids + n_jp).astype(np.int32)
    nss_src = np.concatenate([np.zeros(n_ssp, np.int32), ssn.table_src_sub])
    max_opts = 1 + max(int(np.diff(pj_ptr).max(initial=0)), int(np.diff(nj_ptr).max(initial=0)))
    args = (
        dst.block_psub,
        dst.block_nsub,
        src.block_of,
        pj_ptr,
        pj_ids,
        jp.table_src_sub,
        jp.table_class,
        nj_ptr,
        nj_ids,
        nj_src,
        nj_cls,
        lookup,
        pss_ptr,
        pss_ids,
        ssp.table_src_sub,
        nss_ptr,
        nss_ids,
        nss_src,
        max_opts,
    )
    nb = dst.n_blocks
    dummy64 = np.zeros(nb + 1, np.int64)
    d32 = np.zeros(0, np.int32)
    n_tasks = build_tasks(*args, True, dummy64, d32, d32, d32, d32, np.zeros(0, np.int8))
    t_ptr = np.zeros(nb + 1, np.int64)
    t_src = np.empty(n_tasks, np.int32)
    t_pj = np.empty(n_tasks, np.int32)
    t_nj = np.empty(n_tasks, np.int32)
    t_coup = np.empty(n_tasks, np.int32)
    t_kind = np.empty(n_tasks, np.int8)
    build_tasks(*args, False, t_ptr, t_src, t_pj, t_nj, t_coup, t_kind)
    return OperatorTables(
        dst=dst,
        src=src,
        task_ptr=t_ptr,
        task_src_block=t_src,
        task_pjump=t_pj,
        task_njump=t_nj,
        task_coupling=t_coup,
        task_kind=t_kind,
        jump_table_start=jump_table_start.astype(np.int64),
        jump_row_ptr=jump_row_ptr.astype(np.int64),
        jump_src=jump_src.astype(np.int32),
        jump_op=jump_op.astype(np.int16),
        coup_start=np.array(coup_start, np.int64),
        coup_ncol=np.array(coup_ncol, np.int32),
        coup_data=np.concatenate(data) if data else np.zeros(0),
        ss_start=ss_start.astype(np.int64),
        ss_row_ptr=ss_row_ptr.astype(np.int64),
        ss_col=ss_col.astype(np.int32),
        ss_val=ss_val.astype(np.float64),
    )


def hamiltonian_tables(basis: MBasis, inter: MSchemeInteraction) -> OperatorTables:
    """Tables of H on ``basis``: pn coupling tables, pp and nn same-species tables."""
    space = basis.space
    pn = inter.pn.reshape(-1, 4)  # (cp, ap, cn, an)
    mp = make_moves(space, pn[:, [0, 1]])
    mn = make_moves(space, pn[:, [2, 3]])
    jp = jump_tables(basis.protons, basis.protons, space, mp)
    jn = jump_tables(basis.neutrons, basis.neutrons, space, mn)
    p_index = {
        (int(c), int(a)): (int(k), int(q))
        for c, a, k, q in zip(mp.c, mp.a, mp.cls, mp.label, strict=True)
    }
    n_index = {
        (int(c), int(a)): (int(k), int(q))
        for c, a, k, q in zip(mn.c, mn.a, mn.cls, mn.label, strict=True)
    }
    couplings: dict[tuple[int, int], NDArray[np.float64]] = {}
    for (cp, ap, cn, an), v in zip(pn, inter.pn_v, strict=True):
        pk, q = p_index[(int(cp), int(ap))]
        nk, r = n_index[(int(cn), int(an))]
        V = couplings.get((pk, nk))
        if V is None:
            V = couplings[(pk, nk)] = np.zeros((mp.class_size(pk), mn.class_size(nk)))
        V[q, r] += v
    ssp = species_tables(basis.protons, space, inter.spe_p, inter.pp, inter.pp_v)
    ssn = species_tables(basis.neutrons, space, inter.spe_n, inter.nn, inter.nn_v)
    return assemble(basis, basis, jp, jn, mp.n_class, mn.n_class, couplings, ssp, ssn)


def one_body_tables(
    dst: MBasis,
    src: MBasis,
    p_terms: list[tuple[int, int, float]],
    n_terms: list[tuple[int, int, float]],
) -> OperatorTables:
    """Tables of a particle-number-conserving one-body operator ``sum v a+_c a_a``.

    ``dst`` and ``src`` share their species bases (different total M is allowed, as for J+).
    """
    space = src.space
    couplings: dict[tuple[int, int], NDArray[np.float64]] = {}
    jp = jn = None
    mp = mn = None
    if p_terms:
        mp = make_moves(space, np.array([(c, a) for c, a, _ in p_terms]))
        jp = jump_tables(dst.protons, src.protons, space, mp)
        for c, a, v in p_terms:
            k = int(np.flatnonzero((mp.c == c) & (mp.a == a))[0])
            key = (int(mp.cls[k]), -1)
            couplings.setdefault(key, np.zeros((mp.class_size(key[0]), 1)))[mp.label[k], 0] += v
    if n_terms:
        mn = make_moves(space, np.array([(c, a) for c, a, _ in n_terms]))
        jn = jump_tables(dst.neutrons, src.neutrons, space, mn)
        for c, a, v in n_terms:
            k = int(np.flatnonzero((mn.c == c) & (mn.a == a))[0])
            key = (-1, int(mn.cls[k]))
            couplings.setdefault(key, np.zeros((1, mn.class_size(key[1]))))[0, mn.label[k]] += v
    return assemble(dst, src, jp, jn, mp.n_class if mp else 0, mn.n_class if mn else 0, couplings)


def jplus_tables(basis: MBasis) -> tuple[OperatorTables, MBasis]:
    """Tables of J+ from ``basis`` (2M) to the basis with 2M + 2, which is also returned."""
    space = basis.space
    up = build_mbasis(
        space, basis.protons, basis.neutrons, basis.m2 + 2, basis.parity, basis.truncation
    )
    terms: list[tuple[int, int, float]] = []
    for k in range(space.n_morbits - 1):
        if space.morb_orbit[k] == space.morb_orbit[k + 1]:
            j2, m2 = int(space.morb_j2[k]), int(space.morb_m2[k])
            terms.append((k + 1, k, sqrt((j2 * (j2 + 2) - m2 * (m2 + 2)) / 4.0)))
    return one_body_tables(up, basis, terms, terms), up


def expand_to_csr(t: OperatorTables) -> sp.csr_matrix:
    """Return the stored matrix of the tables (tests only)."""
    dst, src = t.dst, t.src
    args = (
        dst.n_blocks,
        dst.block_offset,
        dst.block_psub,
        dst.block_nsub,
        src.block_offset,
        src.block_psub,
        src.block_nsub,
        dst.protons.sub_dim,
        dst.neutrons.sub_dim,
        src.protons.sub_dim,
        src.neutrons.sub_dim,
        t.task_ptr,
        t.task_src_block,
        t.task_pjump,
        t.task_njump,
        t.task_coupling,
        t.task_kind,
        t.jump_table_start,
        t.jump_row_ptr,
        t.jump_src,
        t.jump_op,
        t.coup_start,
        t.coup_ncol,
        t.coup_data,
        t.ss_start,
        t.ss_row_ptr,
        t.ss_col,
        t.ss_val,
    )
    e64 = np.zeros(0, np.int64)
    n = expand_tasks(*args, True, e64, e64, np.zeros(0))
    rows = np.empty(n, np.int64)
    cols = np.empty(n, np.int64)
    vals = np.empty(n, np.float64)
    expand_tasks(*args, False, rows, cols, vals)
    mat = sp.coo_matrix((vals, (rows, cols)), shape=(dst.dim, src.dim)).tocsr()
    mat.sum_duplicates()
    mat.eliminate_zeros()
    return mat
