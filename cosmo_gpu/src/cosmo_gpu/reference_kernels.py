"""Numba kernels that build a stored sparse operator column by column from bitstrings.

This is the stage-1 oracle: each source basis state is acted on by every operator term whose
annihilators it contains, with no factorization. It must stay obviously correct.

Terms arrive indexed by their annihilated m-orbits (CSR over the annihilator key):

* one-body ``v a+_c a_a``: key ``a``; entries ``c``.
* like-particle two-body ``v a+_c1 a+_c2 a_a2 a_a1`` (``a1 < a2``): key ``a1 * nm + a2``;
  entries ``(c1, c2)``.
* proton-neutron ``v a+_cp a_ap a+_cn a_an``: key ``ap * nm + an``; entries ``(cp, cn)``.

So per state only occupied annihilators are visited. Neutron operators always come in pairs, so
moving them past the proton creation string gives no sign.
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
def _index(p: np.uint64, n: np.uint64, lookup: tuple, dst: tuple) -> int:  # type: ignore[type-arg]
    p_sorted, p_sub, p_local, n_sorted, n_sub, n_local, n_start = lookup
    dst_block_of, dst_offset = dst
    sp, lp = find_det(p, p_sorted, p_sub, p_local)
    sn, ln = find_det(n, n_sorted, n_sub, n_local)
    if sp < 0 or sn < 0:
        return -1
    B = dst_block_of[sp, sn]
    if B < 0:
        return -1
    dn = n_start[sn + 1] - n_start[sn]
    return dst_offset[B] + lp * dn + ln


@njit(cache=True, inline="always")
def _emit(
    idx: int, val: float, k: int, out_idx: NDArray[np.int64], out_val: NDArray[np.float64]
) -> int:
    if idx >= 0:
        out_idx[k] = idx
        out_val[k] = val
        return k + 1
    return k


@njit(cache=True)
def _column(
    col: int,
    nm: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    diag_p: NDArray[np.float64],
    diag_n: NDArray[np.float64],
    ob_p: tuple,  # type: ignore[type-arg]
    ob_n: tuple,  # type: ignore[type-arg]
    tb_p: tuple,  # type: ignore[type-arg]
    tb_n: tuple,  # type: ignore[type-arg]
    tb_pn: tuple,  # type: ignore[type-arg]
    out_idx: NDArray[np.int64],
    out_val: NDArray[np.float64],
) -> int:
    offset, psub, nsub, p_start, n_start, p_dets, n_dets = src
    p, n = _state(col, offset, psub, nsub, p_start, n_start, p_dets, n_dets)
    k = 0

    occ_p = np.empty(nm, dtype=np.int64)
    occ_n = np.empty(nm, dtype=np.int64)
    zp = 0
    zn = 0
    for b in range(nm):
        if (p >> np.uint64(b)) & ONE:
            occ_p[zp] = b
            zp += 1
        if (n >> np.uint64(b)) & ONE:
            occ_n[zn] = b
            zn += 1

    if diag_p.shape[0] > 0:
        e = 0.0
        for i in range(zp):
            e += diag_p[occ_p[i]]
        for i in range(zn):
            e += diag_n[occ_n[i]]
        k = _emit(_index(p, n, lookup, dst), e, k, out_idx, out_val)

    # one-body, protons then neutrons
    ptr, cc, vv = ob_p
    for i in range(zp):
        a = occ_p[i]
        s1, g1 = annihilate(p, a)
        for t in range(ptr[a], ptr[a + 1]):
            s2, g2 = create(s1, cc[t])
            if g2 != 0:
                k = _emit(_index(s2, n, lookup, dst), g1 * g2 * vv[t], k, out_idx, out_val)
    ptr, cc, vv = ob_n
    for i in range(zn):
        a = occ_n[i]
        s1, g1 = annihilate(n, a)
        for t in range(ptr[a], ptr[a + 1]):
            s2, g2 = create(s1, cc[t])
            if g2 != 0:
                k = _emit(_index(p, s2, lookup, dst), g1 * g2 * vv[t], k, out_idx, out_val)

    # two-body, protons
    ptr, cc, vv = tb_p
    for i in range(zp):
        a1 = occ_p[i]
        s1, g1 = annihilate(p, a1)
        for j in range(i + 1, zp):
            a2 = occ_p[j]
            key = a1 * nm + a2
            if ptr[key] == ptr[key + 1]:
                continue
            s2, g2 = annihilate(s1, a2)
            for t in range(ptr[key], ptr[key + 1]):
                s3, g3 = create(s2, cc[t, 1])
                if g3 == 0:
                    continue
                s4, g4 = create(s3, cc[t, 0])
                if g4 == 0:
                    continue
                val = g1 * g2 * g3 * g4 * vv[t]
                k = _emit(_index(s4, n, lookup, dst), val, k, out_idx, out_val)

    # two-body, neutrons
    ptr, cc, vv = tb_n
    for i in range(zn):
        a1 = occ_n[i]
        s1, g1 = annihilate(n, a1)
        for j in range(i + 1, zn):
            a2 = occ_n[j]
            key = a1 * nm + a2
            if ptr[key] == ptr[key + 1]:
                continue
            s2, g2 = annihilate(s1, a2)
            for t in range(ptr[key], ptr[key + 1]):
                s3, g3 = create(s2, cc[t, 1])
                if g3 == 0:
                    continue
                s4, g4 = create(s3, cc[t, 0])
                if g4 == 0:
                    continue
                val = g1 * g2 * g3 * g4 * vv[t]
                k = _emit(_index(p, s4, lookup, dst), val, k, out_idx, out_val)

    # two-body, proton-neutron
    ptr, cc, vv = tb_pn
    for i in range(zp):
        ap = occ_p[i]
        sp1, g1 = annihilate(p, ap)
        for j in range(zn):
            an = occ_n[j]
            key = ap * nm + an
            if ptr[key] == ptr[key + 1]:
                continue
            sn1, g3 = annihilate(n, an)
            for t in range(ptr[key], ptr[key + 1]):
                sp2, g2 = create(sp1, cc[t, 0])
                if g2 == 0:
                    continue
                sn2, g4 = create(sn1, cc[t, 1])
                if g4 == 0:
                    continue
                val = g1 * g2 * g3 * g4 * vv[t]
                k = _emit(_index(sp2, sn2, lookup, dst), val, k, out_idx, out_val)

    return _merge(k, out_idx, out_val)


@njit(cache=True)
def _merge(k: int, out_idx: NDArray[np.int64], out_val: NDArray[np.float64]) -> int:
    """Sort entries by row, sum duplicates and drop zeros; return the new count."""
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
    nm: int,
    buf: int,
    n_chunks: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    diag_p: NDArray[np.float64],
    diag_n: NDArray[np.float64],
    ob_p: tuple,  # type: ignore[type-arg]
    ob_n: tuple,  # type: ignore[type-arg]
    tb_p: tuple,  # type: ignore[type-arg]
    tb_n: tuple,  # type: ignore[type-arg]
    tb_pn: tuple,  # type: ignore[type-arg]
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """Return CSC arrays ``(indptr, indices, data)``, one column per source state.

    Columns are processed in ``n_chunks`` parallel chunks with one scratch buffer of ``buf``
    entries per chunk (an upper bound on the entries of one column before merging).
    """
    counts = np.zeros(n_cols + 1, dtype=np.int64)
    step = (n_cols + n_chunks - 1) // n_chunks
    for ch in prange(n_chunks):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        for c in range(ch * step, min(n_cols, (ch + 1) * step)):
            counts[c + 1] = _column(
                c, nm, src, dst, lookup, diag_p, diag_n, ob_p, ob_n, tb_p, tb_n, tb_pn, bi, bv
            )
    indptr = np.cumsum(counts)
    indices = np.empty(indptr[-1], dtype=np.int64)
    data = np.empty(indptr[-1], dtype=np.float64)
    for ch in prange(n_chunks):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        for c in range(ch * step, min(n_cols, (ch + 1) * step)):
            m = _column(
                c, nm, src, dst, lookup, diag_p, diag_n, ob_p, ob_n, tb_p, tb_n, tb_pn, bi, bv
            )
            indices[indptr[c] : indptr[c] + m] = bi[:m]
            data[indptr[c] : indptr[c] + m] = bv[:m]
    return indptr, indices, data


@njit(cache=True)
def _transfer_column(
    col: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    kind: int,
    tc: NDArray[np.int32],
    ta: NDArray[np.int32],
    tv: NDArray[np.float64],
    out_idx: NDArray[np.int64],
    out_val: NDArray[np.float64],
) -> int:
    offset, psub, nsub, p_start, n_start, p_dets, n_dets = src
    p, n = _state(col, offset, psub, nsub, p_start, n_start, p_dets, n_dets)
    zp = 0
    for b in range(64):
        if (p >> np.uint64(b)) & ONE:
            zp += 1
    k = 0
    for t in range(tv.shape[0]):
        p2, n2, g = p, n, 1
        if kind == 0:  # a+_p c
            p2, g = create(p, tc[t])
        elif kind == 1:  # a_p a
            p2, g = annihilate(p, ta[t])
        elif kind == 2:  # a+_n c, passes Z protons
            n2, g = create(n, tc[t])
            g *= 1 - 2 * (zp % 2)
        elif kind == 3:  # a_n a
            n2, g = annihilate(n, ta[t])
            g *= 1 - 2 * (zp % 2)
        elif kind == 4:  # a+_(p,c) a_(n,a): a_n passes Z protons, then a+_p
            n2, g1 = annihilate(n, ta[t])
            p2, g2 = create(p, tc[t])
            g = g1 * g2 * (1 - 2 * (zp % 2))
        else:  # a+_(n,c) a_(p,a): a_p first, then a+_n passes Z - 1 protons
            p2, g1 = annihilate(p, ta[t])
            n2, g2 = create(n, tc[t])
            g = g1 * g2 * (1 - 2 * ((zp - 1) % 2))
        if g != 0:
            k = _emit(_index(p2, n2, lookup, dst), g * tv[t], k, out_idx, out_val)
    return _merge(k, out_idx, out_val)


@njit(cache=True, parallel=True)
def build_csc_transfer(
    n_cols: int,
    n_chunks: int,
    src: tuple,  # type: ignore[type-arg]
    dst: tuple,  # type: ignore[type-arg]
    lookup: tuple,  # type: ignore[type-arg]
    kind: int,
    tc: NDArray[np.int32],
    ta: NDArray[np.int32],
    tv: NDArray[np.float64],
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]:
    """CSC arrays of a single-nucleon or charge-exchange operator, one column per source state."""
    buf = max(1, tv.shape[0])
    counts = np.zeros(n_cols + 1, dtype=np.int64)
    step = (n_cols + n_chunks - 1) // n_chunks
    for ch in prange(n_chunks):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        for c in range(ch * step, min(n_cols, (ch + 1) * step)):
            counts[c + 1] = _transfer_column(c, src, dst, lookup, kind, tc, ta, tv, bi, bv)
    indptr = np.cumsum(counts)
    indices = np.empty(indptr[-1], dtype=np.int64)
    data = np.empty(indptr[-1], dtype=np.float64)
    for ch in prange(n_chunks):
        bi = np.empty(buf, dtype=np.int64)
        bv = np.empty(buf, dtype=np.float64)
        for c in range(ch * step, min(n_cols, (ch + 1) * step)):
            m = _transfer_column(c, src, dst, lookup, kind, tc, ta, tv, bi, bv)
            indices[indptr[c] : indptr[c] + m] = bi[:m]
            data[indptr[c] : indptr[c] + m] = bv[:m]
    return indptr, indices, data
