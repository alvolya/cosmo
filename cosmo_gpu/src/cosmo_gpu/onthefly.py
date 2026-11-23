"""On-the-fly operators: ``y = O x`` from operator tables, with H never stored.

``OnTheFlyOperator`` wraps an ``OperatorTables`` and the Numba kernel
``onthefly_kernels.apply_tables``. Destination blocks are split into contiguous work chunks of
about equal exact cost (``distribution.block_costs``), several per thread, because Numba
schedules ``prange`` statically.
"""

from __future__ import annotations

import time

import numba
import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.basis import MBasis
from cosmo_gpu.distribution import block_costs, partition_blocks
from cosmo_gpu.interaction import MSchemeInteraction
from cosmo_gpu.onthefly_kernels import apply_tables
from cosmo_gpu.operator_tables import OperatorTables, hamiltonian_tables, jplus_tables
from cosmo_gpu.reference_csr import j_from_jsq

CHUNKS_PER_THREAD = 8


class OnTheFlyOperator:
    """A linear operator from ``tables.src`` to ``tables.dst`` applied from its tables."""

    def __init__(self, tables: OperatorTables, n_chunks: int | None = None) -> None:
        self.tables = tables
        t = tables
        nb = t.dst.n_blocks
        if n_chunks is None:
            n_chunks = CHUNKS_PER_THREAD * numba.get_num_threads()
        n_chunks = max(1, min(n_chunks, nb))
        self.chunk_bounds = (
            partition_blocks(block_costs(t), n_chunks) if nb else np.zeros(1, np.int64)
        )
        self.cost = float(block_costs(t).sum()) if nb else 0.0
        self.setup_seconds = 0.0

    @property
    def shape(self) -> tuple[int, int]:
        """``(D_dst, D_src)``."""
        return self.tables.dst.dim, self.tables.src.dim

    @property
    def dim(self) -> int:
        """Destination dimension (equal to the source one for H)."""
        return self.tables.dst.dim

    def apply(
        self,
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        dst_range: tuple[int, int] | None = None,
        src_range: tuple[int, int] | None = None,
        accumulate: bool = False,
    ) -> None:
        """Raw ABI call: ``x`` holds the source blocks of ``src_range`` only (see the kernel)."""
        t = self.tables
        dst, src = t.dst, t.src
        d0, d1 = dst_range if dst_range is not None else (0, dst.n_blocks)
        s0, s1 = src_range if src_range is not None else (0, src.n_blocks)
        width = int(src.block_offset[s1] - src.block_offset[s0])
        if x.ndim != 2 or x.shape[1] != width or y.shape != (x.shape[0], dst.dim):
            msg = f"x must be (nvec, {width}) and y (nvec, {dst.dim}); got {x.shape}, {y.shape}"
            raise ValueError(msg)
        apply_tables(
            self.chunk_bounds,
            dst.block_offset,
            dst.block_psub,
            dst.block_nsub,
            src.block_offset,
            src.block_nsub,
            dst.protons.sub_dim,
            dst.neutrons.sub_dim,
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
            np.ascontiguousarray(x, dtype=np.float64),
            y,
            d0,
            d1,
            s0,
            s1,
            accumulate,
        )

    def matvec(self, X: NDArray[np.float64]) -> NDArray[np.float64]:
        """``X O^T`` for a block of row vectors (one kernel call over all blocks)."""
        X2 = np.atleast_2d(np.ascontiguousarray(X, dtype=np.float64))
        Y = np.empty((X2.shape[0], self.tables.dst.dim))
        self.apply(X2, Y)
        return Y


def hamiltonian_operator(basis: MBasis, inter: MSchemeInteraction) -> OnTheFlyOperator:
    """H on ``basis`` from its tables (``setup_seconds`` records the table build)."""
    t0 = time.perf_counter()
    op = OnTheFlyOperator(hamiltonian_tables(basis, inter))
    op.setup_seconds = time.perf_counter() - t0
    return op


def jplus_operator(basis: MBasis) -> OnTheFlyOperator:
    """J+ from ``basis`` (2M) to the basis with 2M + 2."""
    tables, _ = jplus_tables(basis)
    return OnTheFlyOperator(tables)


def j_expectation_otf(
    vectors: NDArray[np.float64], jplus: OnTheFlyOperator, m2: int
) -> NDArray[np.float64]:
    """``J`` of each row of ``vectors`` from ``<J^2> = |J+ psi|^2 + M (M + 1)``."""
    M = m2 / 2.0
    jp = jplus.matvec(vectors)
    jsq = np.einsum("ij,ij->i", jp, jp) / np.einsum("ij,ij->i", vectors, vectors) + M * (M + 1)
    return j_from_jsq(np.asarray(jsq, dtype=np.float64))
