"""Numba kernel of the on-the-fly matvec ``y = O x`` from operator tables (design plan ABI).

Contract (``docs/design-plan.md``, "Kernel ABI"): gather only, every output amplitude is
written by one thread; tasks whose source block lies outside ``[src_first, src_last)`` are
skipped; ``accumulate=False`` overwrites the owned destination range, ``True`` adds to it; ``y``
outside ``[dst_first, dst_last)`` is untouched. ``x`` holds only the source blocks
``src_first .. src_last - 1`` (its column 0 is ``src_block_offset[src_first]``).
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange
from numpy.typing import NDArray


@njit(cache=True, parallel=True)
def apply_tables(
    chunk_bounds: NDArray[np.int64],
    dst_block_offset: NDArray[np.int64],
    dst_block_psub: NDArray[np.int32],
    dst_block_nsub: NDArray[np.int32],
    src_block_offset: NDArray[np.int64],
    src_block_nsub: NDArray[np.int32],
    dst_p_dim: NDArray[np.int64],
    dst_n_dim: NDArray[np.int64],
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
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    dst_first: int,
    dst_last: int,
    src_first: int,
    src_last: int,
    accumulate: bool,
) -> None:
    """``y[:, owned] (+)= O x`` for the destination blocks ``dst_first .. dst_last - 1``."""
    nvec = x.shape[0]
    x0 = src_block_offset[src_first]
    n_chunks = chunk_bounds.shape[0] - 1
    for c in prange(n_chunks):
        lo = max(chunk_bounds[c], dst_first)
        hi = min(chunk_bounds[c + 1], dst_last)
        for B in range(lo, hi):
            off_d = dst_block_offset[B]
            dpd = dst_p_dim[dst_block_psub[B]]
            dnd = dst_n_dim[dst_block_nsub[B]]
            if not accumulate:
                for v in range(nvec):
                    for k in range(off_d, off_d + dpd * dnd):
                        y[v, k] = 0.0
            for t in range(task_ptr[B], task_ptr[B + 1]):
                sblk = task_src_block[t]
                if sblk < src_first or sblk >= src_last:
                    continue
                off_s = src_block_offset[sblk] - x0
                dns = src_n_dim[src_block_nsub[sblk]]
                kind = task_kind[t]
                if kind == 1:  # pp: proton CSR row times a neutron spectator run
                    base = ss_start[task_coupling[t]]
                    for r in range(dpd):
                        yd = off_d + r * dnd
                        for e in range(ss_row_ptr[base + r], ss_row_ptr[base + r + 1]):
                            xs = off_s + ss_col[e] * dns
                            val = ss_val[e]
                            for j in range(dnd):
                                for v in range(nvec):
                                    y[v, yd + j] += val * x[v, xs + j]
                elif kind == 2:  # nn: neutron CSR row, protons spectators
                    base = ss_start[task_coupling[t]]
                    for r in range(dnd):
                        for e in range(ss_row_ptr[base + r], ss_row_ptr[base + r + 1]):
                            col = ss_col[e]
                            val = ss_val[e]
                            for i in range(dpd):
                                yd = off_d + i * dnd + r
                                xs = off_s + i * dns + col
                                for v in range(nvec):
                                    y[v, yd] += val * x[v, xs]
                else:  # pn: proton jumps x neutron jumps x coupling (either may be identity)
                    pj = task_pjump[t]
                    nj = task_njump[t]
                    cs = coup_start[task_coupling[t]]
                    ncol = coup_ncol[task_coupling[t]]
                    for ip in range(dpd):
                        if pj >= 0:
                            pb = jump_table_start[pj] + ip
                            p_lo = jump_row_ptr[pb]
                            p_hi = jump_row_ptr[pb + 1]
                        else:
                            p_lo = 0
                            p_hi = 1
                        for pe in range(p_lo, p_hi):
                            if pj >= 0:
                                op = int(jump_op[pe])
                                i_src = jump_src[pe]
                                vrow = cs + (op & 0x7FFF) * ncol
                                sp = -1.0 if (op & 0x8000) else 1.0
                            else:
                                i_src = ip
                                vrow = cs
                                sp = 1.0
                            xrow = off_s + i_src * dns
                            for jn in range(dnd):
                                yd = off_d + ip * dnd + jn
                                if nj >= 0:
                                    nb = jump_table_start[nj] + jn
                                    for ne in range(jump_row_ptr[nb], jump_row_ptr[nb + 1]):
                                        opn = int(jump_op[ne])
                                        val = coup_data[vrow + (opn & 0x7FFF)]
                                        if opn & 0x8000:
                                            val = -val
                                        val *= sp
                                        xs = xrow + jump_src[ne]
                                        for v in range(nvec):
                                            y[v, yd] += val * x[v, xs]
                                else:
                                    val = sp * coup_data[vrow]
                                    xs = xrow + jn
                                    for v in range(nvec):
                                        y[v, yd] += val * x[v, xs]
