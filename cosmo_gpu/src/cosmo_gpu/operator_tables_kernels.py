"""Numba kernels that enumerate operator-table entries from bitstrings, and the test expander.

Entries are produced in a deterministic order (source sub-block, source determinant, then the
order of the move or term lists), and sorted into tables on the Python side.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray

from cosmo_gpu.bits_kernels import ONE, annihilate, create, find_det

MERGE_TOL = 1e-14


@njit(cache=True)
def _occupied(s: np.uint64, nm: int, occ: NDArray[np.int64]) -> int:
    z = 0
    for b in range(nm):
        if (s >> np.uint64(b)) & ONE:
            occ[z] = b
            z += 1
    return z


@njit(cache=True, parallel=True)
def jump_entries(
    src_dets: NDArray[np.uint64],
    nm: int,
    move_ptr: NDArray[np.int64],
    move_c: NDArray[np.int32],
    move_class: NDArray[np.int32],
    move_label: NDArray[np.int32],
    dst_sorted: NDArray[np.uint64],
    dst_sub: NDArray[np.int32],
    dst_local: NDArray[np.int32],
) -> tuple[
    NDArray[np.int64],
    NDArray[np.int32],
    NDArray[np.int32],
    NDArray[np.int32],
    NDArray[np.int32],
    NDArray[np.int8],
]:
    """All single moves ``a+_c a_a`` applied to every source determinant.

    Moves are grouped by annihilated m-orbit ``a`` (``move_ptr``). Returns, per entry, the
    source determinant index, destination sub-block and local index, class, label and sign.
    """
    n = src_dets.shape[0]
    counts = np.zeros(n + 1, dtype=np.int64)
    for i in prange(n):
        occ = np.empty(nm, dtype=np.int64)
        s = src_dets[i]
        z = _occupied(s, nm, occ)
        c_ = 0
        for k in range(z):
            a = occ[k]
            s1, g1 = annihilate(s, a)
            for t in range(move_ptr[a], move_ptr[a + 1]):
                s2, g2 = create(s1, move_c[t])
                if g2 != 0:
                    sub, _ = find_det(s2, dst_sorted, dst_sub, dst_local)
                    if sub >= 0:
                        c_ += 1
        counts[i + 1] = c_
    ptr = np.cumsum(counts)
    total = ptr[-1]
    e_src = np.empty(total, dtype=np.int64)
    e_dsub = np.empty(total, dtype=np.int32)
    e_dloc = np.empty(total, dtype=np.int32)
    e_cls = np.empty(total, dtype=np.int32)
    e_lab = np.empty(total, dtype=np.int32)
    e_sgn = np.empty(total, dtype=np.int8)
    for i in prange(n):
        occ = np.empty(nm, dtype=np.int64)
        s = src_dets[i]
        z = _occupied(s, nm, occ)
        w = ptr[i]
        for k in range(z):
            a = occ[k]
            s1, g1 = annihilate(s, a)
            for t in range(move_ptr[a], move_ptr[a + 1]):
                s2, g2 = create(s1, move_c[t])
                if g2 != 0:
                    sub, loc = find_det(s2, dst_sorted, dst_sub, dst_local)
                    if sub >= 0:
                        e_src[w] = i
                        e_dsub[w] = sub
                        e_dloc[w] = loc
                        e_cls[w] = move_class[t]
                        e_lab[w] = move_label[t]
                        e_sgn[w] = g1 * g2
                        w += 1
    return e_src, e_dsub, e_dloc, e_cls, e_lab, e_sgn


@njit(cache=True)
def _merge_column(k: int, idx: NDArray[np.int64], val: NDArray[np.float64]) -> int:
    """Sort entries by index, sum duplicates and drop zeros; return the new count."""
    if k == 0:
        return 0
    order = np.argsort(idx[:k], kind="mergesort")
    ti = idx[:k][order].copy()
    tv = val[:k][order].copy()
    m = 0
    cur = ti[0]
    acc = 0.0
    for q in range(k):
        if ti[q] != cur:
            if abs(acc) > MERGE_TOL:
                idx[m] = cur
                val[m] = acc
                m += 1
            cur = ti[q]
            acc = 0.0
        acc += tv[q]
    if abs(acc) > MERGE_TOL:
        idx[m] = cur
        val[m] = acc
        m += 1
    return m


@njit(cache=True)
def _species_column(
    s: np.uint64,
    nm: int,
    diag: NDArray[np.float64],
    tb_ptr: NDArray[np.int64],
    tb_c: NDArray[np.int32],
    tb_v: NDArray[np.float64],
    dst_sorted: NDArray[np.uint64],
    dst_sub: NDArray[np.int32],
    dst_local: NDArray[np.int32],
    sub_start: NDArray[np.int64],
    occ: NDArray[np.int64],
    out_idx: NDArray[np.int64],
    out_val: NDArray[np.float64],
) -> int:
    """One column of the single-species operator (diagonal plus like-particle two-body)."""
    z = _occupied(s, nm, occ)
    k = 0
    e = 0.0
    for i in range(z):
        e += diag[occ[i]]
    sub, loc = find_det(s, dst_sorted, dst_sub, dst_local)
    if sub >= 0:
        out_idx[k] = sub_start[sub] + loc
        out_val[k] = e
        k += 1
    for i in range(z):
        a1 = occ[i]
        s1, g1 = annihilate(s, a1)
        for j in range(i + 1, z):
            a2 = occ[j]
            key = a1 * nm + a2
            if tb_ptr[key] == tb_ptr[key + 1]:
                continue
            s2, g2 = annihilate(s1, a2)
            for t in range(tb_ptr[key], tb_ptr[key + 1]):
                s3, g3 = create(s2, tb_c[t, 1])
                if g3 == 0:
                    continue
                s4, g4 = create(s3, tb_c[t, 0])
                if g4 == 0:
                    continue
                sub, loc = find_det(s4, dst_sorted, dst_sub, dst_local)
                if sub >= 0:
                    out_idx[k] = sub_start[sub] + loc
                    out_val[k] = g1 * g2 * g3 * g4 * tb_v[t]
                    k += 1
    return _merge_column(k, out_idx, out_val)


@njit(cache=True, parallel=True)
def species_entries(
    dets: NDArray[np.uint64],
    nm: int,
    buf: int,
    diag: NDArray[np.float64],
    tb_ptr: NDArray[np.int64],
    tb_c: NDArray[np.int32],
    tb_v: NDArray[np.float64],
    sorted_dets: NDArray[np.uint64],
    det_sub: NDArray[np.int32],
    det_local: NDArray[np.int32],
    sub_start: NDArray[np.int64],
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """Single-species operator on one species basis, as CSC ``(col_ptr, row, value)``.

    Rows and columns are global determinant indices of the species (sub-block order).
    """
    n = dets.shape[0]
    counts = np.zeros(n + 1, dtype=np.int64)
    for c in prange(n):
        occ = np.empty(nm, dtype=np.int64)
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        counts[c + 1] = _species_column(
            dets[c],
            nm,
            diag,
            tb_ptr,
            tb_c,
            tb_v,
            sorted_dets,
            det_sub,
            det_local,
            sub_start,
            occ,
            bi,
            bv,
        )
    ptr = np.cumsum(counts)
    rows = np.empty(ptr[-1], dtype=np.int64)
    vals = np.empty(ptr[-1], dtype=np.float64)
    for c in prange(n):
        occ = np.empty(nm, dtype=np.int64)
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        m = _species_column(
            dets[c],
            nm,
            diag,
            tb_ptr,
            tb_c,
            tb_v,
            sorted_dets,
            det_sub,
            det_local,
            sub_start,
            occ,
            bi,
            bv,
        )
        rows[ptr[c] : ptr[c] + m] = bi[:m]
        vals[ptr[c] : ptr[c] + m] = bv[:m]
    return ptr, rows, vals


@njit(cache=True)
def expand_tasks(
    n_dst_blocks: int,
    dst_block_offset: NDArray[np.int64],
    dst_block_psub: NDArray[np.int32],
    dst_block_nsub: NDArray[np.int32],
    src_block_offset: NDArray[np.int64],
    src_block_psub: NDArray[np.int32],
    src_block_nsub: NDArray[np.int32],
    dst_p_dim: NDArray[np.int64],
    dst_n_dim: NDArray[np.int64],
    src_p_dim: NDArray[np.int64],
    src_n_dim: NDArray[np.int64],
    task_ptr: NDArray[np.int64],
    task_src_block: NDArray[np.int32],
    task_pjump: NDArray[np.int32],
    task_njump: NDArray[np.int32],
    task_coupling: NDArray[np.int32],
    task_kind: NDArray[np.int8],
    jump_table_start: NDArray[np.int64],
    jump_row_ptr: NDArray[np.int64],
    jump_src: NDArray[np.int32],
    jump_op: NDArray[np.int16],
    coup_start: NDArray[np.int64],
    coup_ncol: NDArray[np.int32],
    coup_data: NDArray[np.float64],
    ss_start: NDArray[np.int64],
    ss_row_ptr: NDArray[np.int64],
    ss_col: NDArray[np.int32],
    ss_val: NDArray[np.float64],
    count_only: bool,
    out_row: NDArray[np.int64],
    out_col: NDArray[np.int64],
    out_val: NDArray[np.float64],
) -> int:
    """Expand every task into COO entries (row, col, value); deliberately simple loops.

    Jump table ``t`` has its destination rows at ``jump_row_ptr[jump_table_start[t] + r]``
    (``r`` = destination local index), and its entries index ``jump_src``/``jump_op``.
    Same-species table ``u`` has rows at ``ss_row_ptr[ss_start[u] + r]``.
    """
    k = 0
    for B in range(n_dst_blocks):
        off_d = dst_block_offset[B]
        dpd = dst_p_dim[dst_block_psub[B]]
        dnd = dst_n_dim[dst_block_nsub[B]]
        for t in range(task_ptr[B], task_ptr[B + 1]):
            S = task_src_block[t]
            off_s = src_block_offset[S]
            dns = src_n_dim[src_block_nsub[S]]
            kind = task_kind[t]
            if kind != 0:  # 1 = pp, 2 = nn
                u = task_coupling[t]
                rows = dpd if kind == 1 else dnd
                for r in range(rows):
                    base = ss_start[u] + r
                    for e in range(ss_row_ptr[base], ss_row_ptr[base + 1]):
                        c = ss_col[e]
                        v = ss_val[e]
                        if kind == 1:  # pp: rows are proton indices, neutrons spectators
                            for j in range(dnd):
                                if not count_only:
                                    out_row[k] = off_d + r * dnd + j
                                    out_col[k] = off_s + c * dns + j
                                    out_val[k] = v
                                k += 1
                        else:  # nn
                            for i in range(dpd):
                                if not count_only:
                                    out_row[k] = off_d + i * dnd + r
                                    out_col[k] = off_s + i * dns + c
                                    out_val[k] = v
                                k += 1
                continue
            # kind 0: proton jump x neutron jump x coupling (either jump may be identity)
            pj = task_pjump[t]
            nj = task_njump[t]
            cs = coup_start[task_coupling[t]]
            ncol = coup_ncol[task_coupling[t]]
            for ip in range(dpd):
                if pj >= 0:
                    pbase = jump_table_start[pj] + ip
                    p_lo, p_hi = jump_row_ptr[pbase], jump_row_ptr[pbase + 1]
                else:
                    p_lo, p_hi = 0, 1
                for pe in range(p_lo, p_hi):
                    if pj >= 0:
                        i_src = jump_src[pe]
                        q = int(jump_op[pe]) & 0x7FFF
                        sp = -1.0 if (int(jump_op[pe]) & 0x8000) else 1.0
                    else:
                        i_src = ip
                        q = 0
                        sp = 1.0
                    for jn in range(dnd):
                        if nj >= 0:
                            nbase = jump_table_start[nj] + jn
                            n_lo, n_hi = jump_row_ptr[nbase], jump_row_ptr[nbase + 1]
                        else:
                            n_lo, n_hi = 0, 1
                        for ne in range(n_lo, n_hi):
                            if nj >= 0:
                                j_src = jump_src[ne]
                                r_ = int(jump_op[ne]) & 0x7FFF
                                sn = -1.0 if (int(jump_op[ne]) & 0x8000) else 1.0
                            else:
                                j_src = jn
                                r_ = 0
                                sn = 1.0
                            v = sp * sn * coup_data[cs + q * ncol + r_]
                            if v == 0.0:
                                continue
                            if not count_only:
                                out_row[k] = off_d + ip * dnd + jn
                                out_col[k] = off_s + i_src * dns + j_src
                                out_val[k] = v
                            k += 1
    return k


@njit(cache=True)
def _options(
    sub: int,
    ptr: NDArray[np.int64],
    ids: NDArray[np.int32],
    table_src_sub: NDArray[np.int32],
    table_class: NDArray[np.int32],
    identity_class: int,
    use_identity: bool,
    out_table: NDArray[np.int32],
    out_src: NDArray[np.int32],
    out_class: NDArray[np.int32],
) -> int:
    """Jump options of one species for destination sub-block ``sub`` (identity last)."""
    k = 0
    for e in range(ptr[sub], ptr[sub + 1]):
        t = ids[e]
        out_table[k] = t
        out_src[k] = table_src_sub[t]
        out_class[k] = table_class[t]
        k += 1
    if use_identity:
        out_table[k] = -1
        out_src[k] = sub
        out_class[k] = identity_class
        k += 1
    return k


@njit(cache=True)
def build_tasks(
    dst_block_psub: NDArray[np.int32],
    dst_block_nsub: NDArray[np.int32],
    src_block_of: NDArray[np.int32],
    pj_ptr: NDArray[np.int64],
    pj_ids: NDArray[np.int32],
    pj_src: NDArray[np.int32],
    pj_class: NDArray[np.int32],
    nj_ptr: NDArray[np.int64],
    nj_ids: NDArray[np.int32],
    nj_src: NDArray[np.int32],
    nj_class: NDArray[np.int32],
    coup_lookup: NDArray[np.int32],
    pss_ptr: NDArray[np.int64],
    pss_ids: NDArray[np.int32],
    pss_src: NDArray[np.int32],
    nss_ptr: NDArray[np.int64],
    nss_ids: NDArray[np.int32],
    nss_src: NDArray[np.int32],
    max_opts: int,
    count_only: bool,
    t_ptr: NDArray[np.int64],
    t_src: NDArray[np.int32],
    t_pj: NDArray[np.int32],
    t_nj: NDArray[np.int32],
    t_coup: NDArray[np.int32],
    t_kind: NDArray[np.int8],
) -> int:
    """Task list over destination blocks: pp, then nn, then pn tasks, in table order.

    ``coup_lookup[pc, nc]`` is the coupling-table id of proton class ``pc`` with neutron class
    ``nc`` (-1 if none); the last row and column stand for the identity on that species.
    """
    n_pc = coup_lookup.shape[0] - 1
    n_nc = coup_lookup.shape[1] - 1
    use_id_p = False
    use_id_n = False
    for c in range(n_nc):
        if coup_lookup[n_pc, c] >= 0:
            use_id_p = True
    for c in range(n_pc):
        if coup_lookup[c, n_nc] >= 0:
            use_id_n = True
    p_tab = np.empty(max_opts, dtype=np.int32)
    p_src = np.empty(max_opts, dtype=np.int32)
    p_cls = np.empty(max_opts, dtype=np.int32)
    n_tab = np.empty(max_opts, dtype=np.int32)
    n_src = np.empty(max_opts, dtype=np.int32)
    n_cls = np.empty(max_opts, dtype=np.int32)
    k = 0
    for B in range(dst_block_psub.shape[0]):
        ps = dst_block_psub[B]
        ns = dst_block_nsub[B]
        if not count_only:
            t_ptr[B] = k
        for e in range(pss_ptr[ps], pss_ptr[ps + 1]):
            S = src_block_of[pss_src[pss_ids[e]], ns]
            if S >= 0:
                if not count_only:
                    t_src[k] = S
                    t_pj[k] = -1
                    t_nj[k] = -1
                    t_coup[k] = pss_ids[e]
                    t_kind[k] = 1
                k += 1
        for e in range(nss_ptr[ns], nss_ptr[ns + 1]):
            S = src_block_of[ps, nss_src[nss_ids[e]]]
            if S >= 0:
                if not count_only:
                    t_src[k] = S
                    t_pj[k] = -1
                    t_nj[k] = -1
                    t_coup[k] = nss_ids[e]
                    t_kind[k] = 2
                k += 1
        npo = _options(ps, pj_ptr, pj_ids, pj_src, pj_class, n_pc, use_id_p, p_tab, p_src, p_cls)
        nno = _options(ns, nj_ptr, nj_ids, nj_src, nj_class, n_nc, use_id_n, n_tab, n_src, n_cls)
        for a in range(npo):
            for b in range(nno):
                if p_tab[a] < 0 and n_tab[b] < 0:
                    continue
                cid = coup_lookup[p_cls[a], n_cls[b]]
                if cid < 0:
                    continue
                S = src_block_of[p_src[a], n_src[b]]
                if S < 0:
                    continue
                if not count_only:
                    t_src[k] = S
                    t_pj[k] = p_tab[a]
                    t_nj[k] = n_tab[b]
                    t_coup[k] = cid
                    t_kind[k] = 0
                k += 1
    if not count_only:
        t_ptr[dst_block_psub.shape[0]] = k
    return k
