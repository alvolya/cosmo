"""Block ownership across ranks and the source data each rank needs from the others.

The vector is split by destination block: rank ``r`` owns a contiguous range of blocks of both
``x`` and ``y``. A matvec on rank ``r`` needs every source block that a task of its destination
blocks reads; the ones owned elsewhere must be received. That count, the *remote fraction*, is
computed here from the task list before any GPU code exists.

Task cost is the exact number of multiply-adds of its gather: a pn task costs
``nnz(proton jumps) * nnz(neutron jumps)`` (rows factorize), with a spectator dimension in place
of an identity; a pp (nn) task costs ``nnz(table) * dn`` (``* dp``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.operator_tables import OperatorTables


def _table_nnz(
    start: NDArray[np.int64], row_ptr: NDArray[np.int64], n_rows: NDArray[np.int64]
) -> NDArray[np.int64]:
    return np.asarray(row_ptr[start + n_rows] - row_ptr[start], dtype=np.int64)


def task_costs(t: OperatorTables) -> NDArray[np.float64]:
    """Multiply-adds of each task (per vector)."""
    dst = t.dst
    nb = dst.n_blocks
    task_block = np.repeat(np.arange(nb), np.diff(t.task_ptr))
    dp = dst.protons.sub_dim[dst.block_psub[task_block]].astype(np.float64)
    dn = dst.neutrons.sub_dim[dst.block_nsub[task_block]].astype(np.float64)
    # rows of each jump table: destination sub-block dimension of its species
    n_jump = len(t.jump_table_start)
    jrows = np.diff(np.append(t.jump_table_start, len(t.jump_row_ptr) - 1))
    jnnz = _table_nnz(t.jump_table_start, t.jump_row_ptr, jrows).astype(np.float64)
    n_ss = len(t.ss_start)
    srows = np.diff(np.append(t.ss_start, len(t.ss_row_ptr) - 1))
    snnz = _table_nnz(t.ss_start, t.ss_row_ptr, srows).astype(np.float64)
    cost = np.zeros(t.n_tasks)
    kind = t.task_kind
    pn = kind == 0
    pj, nj = t.task_pjump, t.task_njump
    p_part = np.where(pj >= 0, jnnz[np.clip(pj, 0, max(n_jump - 1, 0))] if n_jump else 0.0, dp)
    n_part = np.where(nj >= 0, jnnz[np.clip(nj, 0, max(n_jump - 1, 0))] if n_jump else 0.0, dn)
    cost[pn] = (p_part * n_part)[pn]
    if n_ss:
        ss = snnz[np.clip(t.task_coupling, 0, n_ss - 1)]
        cost[kind == 1] = (ss * dn)[kind == 1]
        cost[kind == 2] = (ss * dp)[kind == 2]
    return cost


def block_costs(t: OperatorTables) -> NDArray[np.float64]:
    """Multiply-adds per destination block."""
    task_block = np.repeat(np.arange(t.dst.n_blocks), np.diff(t.task_ptr))
    return np.bincount(task_block, weights=task_costs(t), minlength=t.dst.n_blocks)


def partition_blocks(cost: NDArray[np.float64], n_ranks: int) -> NDArray[np.int64]:
    """Split blocks into contiguous ranges of about equal cost.

    Rank ``r`` owns blocks ``[bounds[r], bounds[r+1])``.
    """
    cum = np.concatenate([[0.0], np.cumsum(cost)])
    targets = cum[-1] * np.arange(1, n_ranks) / n_ranks
    cuts = np.searchsorted(cum, targets, side="left")
    return np.concatenate([[0], cuts, [len(cost)]]).astype(np.int64)


@dataclass(frozen=True)
class RemoteReport:
    """Source data each rank receives in one matvec (per vector, FP64)."""

    n_ranks: int
    bounds: NDArray[np.int64]
    owned_states: NDArray[np.int64]
    remote_states: NDArray[np.int64]
    cost_share: NDArray[np.float64]

    @property
    def remote_fraction(self) -> NDArray[np.float64]:
        """Received states over the full dimension, per rank."""
        return self.remote_states / max(1, int(self.owned_states.sum()))

    @property
    def remote_bytes(self) -> NDArray[np.int64]:
        """Received bytes per rank per vector."""
        return self.remote_states * 8


def remote_fraction(t: OperatorTables, n_ranks: int) -> RemoteReport:
    """Partition ``t``'s destination blocks over ``n_ranks`` and count remote source states."""
    if t.dst is not t.src:
        msg = "remote fraction is defined for square operators (x and y split alike)"
        raise ValueError(msg)
    cost = block_costs(t)
    bounds = partition_blocks(cost, n_ranks)
    owner = np.searchsorted(bounds, np.arange(t.dst.n_blocks), side="right") - 1
    sizes = np.diff(t.dst.block_offset)
    owned = np.array([sizes[bounds[r] : bounds[r + 1]].sum() for r in range(n_ranks)], np.int64)
    remote = np.zeros(n_ranks, np.int64)
    share = np.zeros(n_ranks)
    total = max(float(cost.sum()), 1.0)
    for r in range(n_ranks):
        lo, hi = t.task_ptr[bounds[r]], t.task_ptr[bounds[r + 1]]
        need = np.unique(t.task_src_block[lo:hi])
        foreign = need[owner[need] != r]
        remote[r] = int(sizes[foreign].sum())
        share[r] = float(cost[bounds[r] : bounds[r + 1]].sum()) / total
    return RemoteReport(n_ranks, bounds, owned, remote, share)
