"""Numba kernels that build a stored sparse operator column by column from bitstrings.

This is the stage-1 oracle: every term of the operator is applied to every source basis state,
with no factorization. It is slow and memory-hungry by design and must stay obviously correct.

Operator terms (m-orbit indices, ``c`` created, ``a`` annihilated):

* ``diag``: one-body diagonal energies per m-orbit, the same for both species.
* ``ob_p`` / ``ob_n``: ``v a+_c a_a`` on protons / neutrons, rows ``(c, a)``.
* ``tb_p`` / ``tb_n``: ``v a+_c1 a+_c2 a_a2 a_a1``, rows ``(c1, c2, a1, a2)``.
* ``tb_pn``: ``v a+_cp a_ap a+_cn a_an``, rows ``(cp, ap, cn, an)``.

Neutron operators always come in pairs, so passing them through the proton creation string gives
no sign.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray

from cosmo_gpu.bits_kernels import ONE, annihilate, create, find_det

DROP_TOL = 1e-14


@njit(cache=True, inline="always")
def _state(
    col: int,
    offset: NDArray[np.int64],
    psub: NDArray[np.int32],
    nsub: NDArray[np.int32],
    p_start: NDArray[np.int64],
    n_start: NDArray[np.int64],
    p_dets: NDArray[np.uint64],
    n_dets: NDArray[np.uint64],
) -> tuple[np.uint64, np.uint64]:
    B = np.searchsorted(offset, col, side="right") - 1
    dn = n_start[nsub[B] + 1] - n_start[nsub[B]]
    r = col - offset[B]
    ip = r // dn
    i_n = r - ip * dn
    return p_dets[p_start[psub[B]] + ip], n_dets[n_start[nsub[B]] + i_n]


@njit(cache=True, inline="always")
def _index(
    p: np.uint64,
    n: np.uint64,
    p_sorted: NDArray[np.uint64],
    p_sub: NDArray[np.int32],
    p_local: NDArray[np.int32],
    n_sorted: NDArray[np.uint64],
    n_sub: NDArray[np.int32],
    n_local: NDArray[np.int32],
    n_start: NDArray[np.int64],
    dst_block_of: NDArray[np.int32],
    dst_offset: NDArray[np.int64],
) -> int:
    sp, lp = find_det(p, p_sorted, p_sub, p_local)
    sn, ln = find_det(n, n_sorted, n_sub, n_local)
    if sp < 0 or sn < 0:
        return -1
    B = dst_block_of[sp, sn]
    if B < 0:
        return -1
    dn = n_start[sn + 1] - n_start[sn]
    return dst_offset[B] + lp * dn + ln


@njit(cache=True)
def _column(
    col: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    diag: NDArray[np.float64],
    ob_p: NDArray[np.int32],
    ob_p_v: NDArray[np.float64],
    ob_n: NDArray[np.int32],
    ob_n_v: NDArray[np.float64],
    tb_p: NDArray[np.int32],
    tb_p_v: NDArray[np.float64],
    tb_n: NDArray[np.int32],
    tb_n_v: NDArray[np.float64],
    tb_pn: NDArray[np.int32],
    tb_pn_v: NDArray[np.float64],
    out_idx: NDArray[np.int64],
    out_val: NDArray[np.float64],
) -> int:
    offset, psub, nsub, p_start, n_start, p_dets, n_dets = src
    dst_block_of, dst_offset = dst
    p_sorted, p_sub, p_local, n_sorted, n_sub, n_local = lookup
    p, n = _state(col, offset, psub, nsub, p_start, n_start, p_dets, n_dets)
    k = 0

    # diagonal one-body energies
    if diag.shape[0] > 0:
        e = 0.0
        for b in range(diag.shape[0]):
            if (p >> np.uint64(b)) & ONE:
                e += diag[b]
            if (n >> np.uint64(b)) & ONE:
                e += diag[b]
        idx = _index(
            p,
            n,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = e
            k += 1

    # one-body, protons
    for t in range(ob_p.shape[0]):
        s, sg1 = annihilate(p, ob_p[t, 1])
        if sg1 == 0:
            continue
        s, sg2 = create(s, ob_p[t, 0])
        if sg2 == 0:
            continue
        idx = _index(
            s,
            n,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = sg1 * sg2 * ob_p_v[t]
            k += 1

    # one-body, neutrons
    for t in range(ob_n.shape[0]):
        s, sg1 = annihilate(n, ob_n[t, 1])
        if sg1 == 0:
            continue
        s, sg2 = create(s, ob_n[t, 0])
        if sg2 == 0:
            continue
        idx = _index(
            p,
            s,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = sg1 * sg2 * ob_n_v[t]
            k += 1

    # two-body, protons
    for t in range(tb_p.shape[0]):
        s, sg1 = annihilate(p, tb_p[t, 2])
        if sg1 == 0:
            continue
        s, sg2 = annihilate(s, tb_p[t, 3])
        if sg2 == 0:
            continue
        s, sg3 = create(s, tb_p[t, 1])
        if sg3 == 0:
            continue
        s, sg4 = create(s, tb_p[t, 0])
        if sg4 == 0:
            continue
        idx = _index(
            s,
            n,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = sg1 * sg2 * sg3 * sg4 * tb_p_v[t]
            k += 1

    # two-body, neutrons
    for t in range(tb_n.shape[0]):
        s, sg1 = annihilate(n, tb_n[t, 2])
        if sg1 == 0:
            continue
        s, sg2 = annihilate(s, tb_n[t, 3])
        if sg2 == 0:
            continue
        s, sg3 = create(s, tb_n[t, 1])
        if sg3 == 0:
            continue
        s, sg4 = create(s, tb_n[t, 0])
        if sg4 == 0:
            continue
        idx = _index(
            p,
            s,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = sg1 * sg2 * sg3 * sg4 * tb_n_v[t]
            k += 1

    # two-body, proton-neutron
    for t in range(tb_pn.shape[0]):
        sp, sg1 = annihilate(p, tb_pn[t, 1])
        if sg1 == 0:
            continue
        sp, sg2 = create(sp, tb_pn[t, 0])
        if sg2 == 0:
            continue
        sn, sg3 = annihilate(n, tb_pn[t, 3])
        if sg3 == 0:
            continue
        sn, sg4 = create(sn, tb_pn[t, 2])
        if sg4 == 0:
            continue
        idx = _index(
            sp,
            sn,
            p_sorted,
            p_sub,
            p_local,
            n_sorted,
            n_sub,
            n_local,
            n_start,
            dst_block_of,
            dst_offset,
        )
        if idx >= 0:
            out_idx[k] = idx
            out_val[k] = sg1 * sg2 * sg3 * sg4 * tb_pn_v[t]
            k += 1

    # merge duplicates
    if k == 0:
        return 0
    order = np.argsort(out_idx[:k])
    tmp_i = out_idx[:k][order].copy()
    tmp_v = out_val[:k][order].copy()
    m = 0
    cur = tmp_i[0]
    acc = 0.0
    for q in range(k):
        if tmp_i[q] != cur:
            if abs(acc) > DROP_TOL:
                out_idx[m] = cur
                out_val[m] = acc
                m += 1
            cur = tmp_i[q]
            acc = 0.0
        acc += tmp_v[q]
    if abs(acc) > DROP_TOL:
        out_idx[m] = cur
        out_val[m] = acc
        m += 1
    return m


@njit(cache=True, parallel=True)
def build_csc(
    n_cols: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    diag: NDArray[np.float64],
    ob_p: NDArray[np.int32],
    ob_p_v: NDArray[np.float64],
    ob_n: NDArray[np.int32],
    ob_n_v: NDArray[np.float64],
    tb_p: NDArray[np.int32],
    tb_p_v: NDArray[np.float64],
    tb_n: NDArray[np.int32],
    tb_n_v: NDArray[np.float64],
    tb_pn: NDArray[np.int32],
    tb_pn_v: NDArray[np.float64],
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """Return CSC arrays ``(indptr, indices, data)``, one column per source state."""
    buf = 1 + ob_p.shape[0] + ob_n.shape[0] + tb_p.shape[0] + tb_n.shape[0] + tb_pn.shape[0]
    counts = np.zeros(n_cols + 1, dtype=np.int64)
    for c in prange(n_cols):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        counts[c + 1] = _column(
            c,
            src,
            dst,
            lookup,
            diag,
            ob_p,
            ob_p_v,
            ob_n,
            ob_n_v,
            tb_p,
            tb_p_v,
            tb_n,
            tb_n_v,
            tb_pn,
            tb_pn_v,
            bi,
            bv,
        )
    indptr = np.cumsum(counts)
    indices = np.empty(indptr[-1], dtype=np.int64)
    data = np.empty(indptr[-1], dtype=np.float64)
    for c in prange(n_cols):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        m = _column(
            c,
            src,
            dst,
            lookup,
            diag,
            ob_p,
            ob_p_v,
            ob_n,
            ob_n_v,
            tb_p,
            tb_p_v,
            tb_n,
            tb_n_v,
            tb_pn,
            tb_pn_v,
            bi,
            bv,
        )
        indices[indptr[c] : indptr[c] + m] = bi[:m]
        data[indptr[c] : indptr[c] + m] = bv[:m]
    return indptr, indices, data
