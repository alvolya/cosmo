"""GPU on-the-fly operator: the stage-3 kernel ABI as one CUDA C kernel under CuPy.

Mapping: a work item is 32 consecutive output amplitudes of one destination block, one per
lane of a warp (amplitude ``k = ip * dn + jn`` of the block). Every lane walks the same task
list, so branches on the task kind are uniform across the warp; each lane sums its own output
in registers over all tasks (pn, pp, nn) and writes it once. This is gather only: no atomics,
a fixed summation order (task, proton entry, neutron entry), so reruns are bitwise equal.

The default "lean" pn path computes per-vector base pointers once per task, uses 32-bit offsets
inside a block, hoists the identity cases out of the inner loops and applies signs by flipping
the sign bit; jump entries are packed into one 64-bit load on the device; and the launch caps
registers (``default_launch``) to raise occupancy. Nsight Compute showed the original kernel
issue-bound (about 110 warp instructions per pn entry step, 6% of DRAM bandwidth, 41%
occupancy); these changes are 1.4 to 1.6x faster on H200 and A100. ``lean=False,
packed=False, threads=128, min_blocks=1`` is the original kernel, kept as a reference.

The vector layout is a compile-time choice: ``"rows"`` is the ABI's ``(nvec, D)``;
``"interleaved"`` is ``(D, nvec)`` (measured slower on H200; kept for experiments; the lean
path applies to "rows" only).
"""

from __future__ import annotations

import time
from functools import cache

import cupy as cp
import numpy as np
from numpy.typing import NDArray

from cosmo_gpu.backend_cupy import DeviceArray
from cosmo_gpu.basis import MBasis
from cosmo_gpu.interaction import MSchemeInteraction
from cosmo_gpu.operator_tables import OperatorTables, hamiltonian_tables, jplus_tables
from cosmo_gpu.reference_csr import j_from_jsq

WARP = 32
THREADS = 128
CAP_MAX_NVEC = 4  # the register cap of default_launch applies up to this many vectors

KERNEL = r"""
#ifndef PACKED
#define PACKED 0
#endif
#ifndef LEAN
#define LEAN 0
#endif
#ifndef MINBLOCKS
#define MINBLOCKS 1
#endif
#ifndef TPB
#define TPB 128
#endif

// flip the sign of v when bit 15 of the label is set, without a branch
__device__ __forceinline__ double signed_by(double v, int o) {
    return __longlong_as_double(__double_as_longlong(v) ^ ((long long)(o & 0x8000) << 48));
}

// Jump entry e: source local index and label (sign in bit 15).
#if PACKED
#define JGET(e, s, o) { \
    const unsigned long long w_ = jump_pk[e]; s = (int)(w_ >> 16); o = (int)(w_ & 0xFFFFull); }
#else
#define JGET(e, s, o) { s = jump_src[e]; o = (int)(unsigned short)jump_op[e]; }
#endif

#if INTERLEAVED
#define XAT(v, pos) x[(long long)(pos) * NVEC + (v)]
#define YAT(v, pos) y[(long long)(pos) * NVEC + (v)]
#else
#define XAT(v, pos) x[(long long)(v) * x_len + (pos)]
#define YAT(v, pos) y[(long long)(v) * y_len + (pos)]
#endif

extern "C" __global__ void __launch_bounds__(TPB, MINBLOCKS) apply_tables(
    const int* __restrict__ item_block, const int* __restrict__ item_k0, const long long n_items,
    const long long* __restrict__ dst_block_offset, const int* __restrict__ dst_block_psub,
    const int* __restrict__ dst_block_nsub, const long long* __restrict__ src_block_offset,
    const int* __restrict__ src_block_nsub, const long long* __restrict__ dst_p_dim,
    const long long* __restrict__ dst_n_dim, const long long* __restrict__ src_n_dim,
    const long long* __restrict__ task_ptr, const int* __restrict__ task_src_block,
    const int* __restrict__ task_pjump, const int* __restrict__ task_njump,
    const int* __restrict__ task_coupling, const signed char* __restrict__ task_kind,
    const long long* __restrict__ jump_table_start, const long long* __restrict__ jump_row_ptr,
    const int* __restrict__ jump_src, const short* __restrict__ jump_op,
    const unsigned long long* __restrict__ jump_pk,
    const long long* __restrict__ coup_start, const int* __restrict__ coup_ncol,
    const double* __restrict__ coup_data,
    const long long* __restrict__ ss_start, const long long* __restrict__ ss_row_ptr,
    const int* __restrict__ ss_col, const double* __restrict__ ss_val,
    const double* __restrict__ x, double* __restrict__ y,
    const long long x_len, const long long y_len, const long long x0,
    const int src_first, const int src_last, const int accumulate)
{
    const int lane = threadIdx.x & 31;
    const long long warp = ((long long)blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const long long n_warps = ((long long)gridDim.x * blockDim.x) >> 5;
    for (long long it = warp; it < n_items; it += n_warps) {
        const int B = item_block[it];
        const long long dpd = dst_p_dim[dst_block_psub[B]];
        const long long dnd = dst_n_dim[dst_block_nsub[B]];
        const long long k = (long long)item_k0[it] + lane;
        if (k >= dpd * dnd) continue;
        const long long ip = k / dnd;
        const long long jn = k - ip * dnd;
        const long long out = dst_block_offset[B] + k;
        double acc[NVEC];
        #pragma unroll
        for (int v = 0; v < NVEC; ++v) acc[v] = 0.0;
        for (long long t = task_ptr[B]; t < task_ptr[B + 1]; ++t) {
            const int sb = task_src_block[t];
            if (sb < src_first || sb >= src_last) continue;
            const long long off_s = src_block_offset[sb] - x0;
            const long long dns = src_n_dim[src_block_nsub[sb]];
            const int kind = task_kind[t];
            if (kind == 1) {            // pp: proton CSR row ip, neutron spectator jn
                const long long base = ss_start[task_coupling[t]] + ip;
                for (long long e = ss_row_ptr[base]; e < ss_row_ptr[base + 1]; ++e) {
                    const double val = ss_val[e];
                    const long long xs = off_s + (long long)ss_col[e] * dns + jn;
                    #pragma unroll
                    for (int v = 0; v < NVEC; ++v) acc[v] += val * XAT(v, xs);
                }
            } else if (kind == 2) {     // nn: neutron CSR row jn, proton spectator ip
                const long long base = ss_start[task_coupling[t]] + jn;
                for (long long e = ss_row_ptr[base]; e < ss_row_ptr[base + 1]; ++e) {
                    const double val = ss_val[e];
                    const long long xs = off_s + ip * dns + ss_col[e];
                    #pragma unroll
                    for (int v = 0; v < NVEC; ++v) acc[v] += val * XAT(v, xs);
                }
            } else {                    // pn: proton jumps x neutron jumps x coupling
                const int pj = task_pjump[t];
                const int nj = task_njump[t];
                const long long cs = coup_start[task_coupling[t]];
                const long long ncol = coup_ncol[task_coupling[t]];
                long long p_lo = 0, p_hi = 1;
                if (pj >= 0) {
                    const long long pb = jump_table_start[pj] + ip;
                    p_lo = jump_row_ptr[pb];
                    p_hi = jump_row_ptr[pb + 1];
                }
                long long n_lo = 0, n_hi = 1;
                if (nj >= 0) {
                    const long long nb = jump_table_start[nj] + jn;
                    n_lo = jump_row_ptr[nb];
                    n_hi = jump_row_ptr[nb + 1];
                }
#if LEAN && !INTERLEAVED
                // lean path: per-vector base pointers once per task, 32-bit offsets inside the
                // block, identity cases hoisted out of the loops, signs applied by bit flips
                const double* xb[NVEC];
                #pragma unroll
                for (int v = 0; v < NVEC; ++v) xb[v] = x + (long long)v * x_len + off_s;
                const int dnsi = (int)dns;
                if (pj >= 0 && nj >= 0) {
                    for (long long pe = p_lo; pe < p_hi; ++pe) {
                        int i_src, o;
                        JGET(pe, i_src, o);
                        const double* __restrict__ vr =
                            coup_data + cs + (long long)(o & 0x7FFF) * ncol;
                        const int prow = i_src * dnsi;
                        const int sgn = o & 0x8000;
                        for (long long ne = n_lo; ne < n_hi; ++ne) {
                            int j_src, on;
                            JGET(ne, j_src, on);
                            const double val = signed_by(vr[on & 0x7FFF], (on ^ sgn) & 0x8000);
                            const int xs = prow + j_src;
                            #pragma unroll
                            for (int v = 0; v < NVEC; ++v) acc[v] = fma(val, xb[v][xs], acc[v]);
                        }
                    }
                } else if (pj >= 0) {           // neutron identity: j_src = jn, label 0
                    for (long long pe = p_lo; pe < p_hi; ++pe) {
                        int i_src, o;
                        JGET(pe, i_src, o);
                        const double val =
                            signed_by(coup_data[cs + (long long)(o & 0x7FFF) * ncol], o);
                        const int xs = i_src * dnsi + (int)jn;
                        #pragma unroll
                        for (int v = 0; v < NVEC; ++v) acc[v] = fma(val, xb[v][xs], acc[v]);
                    }
                } else if (nj >= 0) {           // proton identity: i_src = ip, label 0
                    const int prow = (int)ip * dnsi;
                    for (long long ne = n_lo; ne < n_hi; ++ne) {
                        int j_src, on;
                        JGET(ne, j_src, on);
                        const double val = signed_by(coup_data[cs + (on & 0x7FFF)], on);
                        const int xs = prow + j_src;
                        #pragma unroll
                        for (int v = 0; v < NVEC; ++v) acc[v] = fma(val, xb[v][xs], acc[v]);
                    }
                }
#else
                for (long long pe = p_lo; pe < p_hi; ++pe) {
                    long long i_src = ip, vrow = cs;
                    double sp = 1.0;
                    if (pj >= 0) {
                        int s_, o;
                        JGET(pe, s_, o);
                        i_src = s_;
                        vrow = cs + (long long)(o & 0x7FFF) * ncol;
                        if (o & 0x8000) sp = -1.0;
                    }
                    const long long xrow = off_s + i_src * dns;
                    for (long long ne = n_lo; ne < n_hi; ++ne) {
                        long long j_src = jn;
                        double val = sp * coup_data[vrow];
                        if (nj >= 0) {
                            int s_, o;
                            JGET(ne, s_, o);
                            j_src = s_;
                            val = sp * coup_data[vrow + (o & 0x7FFF)];
                            if (o & 0x8000) val = -val;
                        }
                        const long long xs = xrow + j_src;
                        #pragma unroll
                        for (int v = 0; v < NVEC; ++v) acc[v] += val * XAT(v, xs);
                    }
                }
#endif
            }
        }
        #pragma unroll
        for (int v = 0; v < NVEC; ++v) {
            if (accumulate) YAT(v, out) += acc[v];
            else YAT(v, out) = acc[v];
        }
    }
}
"""


@cache
def _kernel(
    nvec: int,
    interleaved: bool,
    packed: bool,
    lean: bool = False,
    threads: int = THREADS,
    min_blocks: int = 1,
) -> cp.RawKernel:
    """Compile (once per vector count, layout and variant) and return the kernel."""
    mod = cp.RawModule(
        code=KERNEL,
        options=(
            "-std=c++14",
            f"-DNVEC={nvec}",
            f"-DINTERLEAVED={int(interleaved)}",
            f"-DPACKED={int(packed)}",
            f"-DLEAN={int(lean)}",
            f"-DTPB={threads}",
            f"-DMINBLOCKS={min_blocks}",
        ),
    )
    return mod.get_function("apply_tables")


def default_launch() -> tuple[int, int]:
    """Threads per block and minimum resident blocks (register cap) for this GPU.

    Measured on 28Si spsdpf 1hw, 24Mg 1hw and 48Cr (Nov 2026): 256 threads with at least 3
    resident blocks per SM is fastest on H200 (compute capability 9.0), with 4 on A100 (8.0),
    both 1.5 to 1.6x faster than the unconstrained 128-thread launch.
    """
    cc = int(cp.cuda.Device().compute_capability)  # e.g. "90" for H200, "80" for A100
    return (256, 3) if cc >= 90 else (256, 4)


def _work_items(
    t: OperatorTables,
) -> tuple[NDArray[np.int32], NDArray[np.int32], NDArray[np.int64]]:
    """Items (block, first amplitude) of 32 outputs each, in block order; and items per block."""
    dst = t.dst
    sizes = (dst.protons.sub_dim[dst.block_psub] * dst.neutrons.sub_dim[dst.block_nsub]).astype(
        np.int64
    )
    n_items = (sizes + WARP - 1) // WARP
    block = np.repeat(np.arange(dst.n_blocks, dtype=np.int32), n_items)
    first = np.concatenate([[0], np.cumsum(n_items)]).astype(np.int64)
    k0 = (np.arange(len(block), dtype=np.int64) - first[block]) * WARP
    return block, k0.astype(np.int32), first


class CudaOperator:
    """A linear operator applied on the GPU from its tables (uploaded once)."""

    def __init__(
        self,
        tables: OperatorTables,
        layout: str = "rows",
        packed: bool = True,
        threads: int | None = None,
        lean: bool = True,
        min_blocks: int | None = None,
    ) -> None:
        """Upload ``tables``; the keyword options select kernel variants (see the module)."""
        if layout not in ("rows", "interleaved"):
            msg = f"layout must be 'rows' or 'interleaved', not {layout!r}"
            raise ValueError(msg)
        t0 = time.perf_counter()
        self.tables = tables
        self.layout = layout
        self.packed = packed
        d_threads, d_min_blocks = default_launch()
        self.threads = d_threads if threads is None else threads
        self.lean = lean
        self.min_blocks = d_min_blocks if min_blocks is None else min_blocks
        t = tables
        dst, src = t.dst, t.src
        block, k0, self.items_first = _work_items(t)
        dev = cp.asarray

        def i32(a: DeviceArray) -> DeviceArray:
            return dev(np.ascontiguousarray(a, dtype=np.int32))

        def i64(a: DeviceArray) -> DeviceArray:
            return dev(np.ascontiguousarray(a, dtype=np.int64))

        self._item_block = i32(block)
        self._item_k0 = i32(k0)
        pk = (np.asarray(t.jump_src, dtype=np.uint64) << np.uint64(16)) | np.asarray(
            t.jump_op, dtype=np.int16
        ).view(np.uint16).astype(np.uint64)
        self._jump_pk = dev(pk) if packed else dev(np.zeros(1, dtype=np.uint64))
        self._args = (
            i64(dst.block_offset),
            i32(dst.block_psub),
            i32(dst.block_nsub),
            i64(src.block_offset),
            i32(src.block_nsub),
            i64(dst.protons.sub_dim),
            i64(dst.neutrons.sub_dim),
            i64(src.neutrons.sub_dim),
            i64(t.task_ptr),
            i32(t.task_src_block),
            i32(t.task_pjump),
            i32(t.task_njump),
            i32(t.task_coupling),
            dev(np.ascontiguousarray(t.task_kind, dtype=np.int8)),
            i64(t.jump_table_start),
            i64(t.jump_row_ptr),
            # the packed kernel reads only jump_pk: do not keep the unpacked copy on the device
            i32(t.jump_src if not packed else np.zeros(1, np.int32)),
            dev(np.ascontiguousarray(t.jump_op if not packed else np.zeros(1), dtype=np.int16)),
            self._jump_pk,
            i64(t.coup_start),
            i32(t.coup_ncol),
            dev(np.ascontiguousarray(t.coup_data, dtype=np.float64)),
            i64(t.ss_start),
            i64(t.ss_row_ptr),
            i32(t.ss_col),
            dev(np.ascontiguousarray(t.ss_val, dtype=np.float64)),
        )
        cp.cuda.Device().synchronize()
        self.upload_seconds = time.perf_counter() - t0
        self.setup_seconds = self.upload_seconds
        self.device_bytes = sum(int(a.nbytes) for a in self._args) + int(
            self._item_block.nbytes + self._item_k0.nbytes
        )

    @property
    def shape(self) -> tuple[int, int]:
        """``(D_dst, D_src)``."""
        return self.tables.dst.dim, self.tables.src.dim

    @property
    def dim(self) -> int:
        """Destination dimension."""
        return self.tables.dst.dim

    def apply(
        self,
        x: DeviceArray,
        y: DeviceArray,
        dst_range: tuple[int, int] | None = None,
        src_range: tuple[int, int] | None = None,
        accumulate: bool = False,
    ) -> None:
        """Raw ABI call on device arrays.

        ``x`` holds the source blocks of ``src_range`` only; layout ``(nvec, width)`` for "rows"
        or ``(width, nvec)`` for "interleaved", and ``y`` likewise with the full destination
        dimension. Outputs outside ``dst_range`` are untouched.
        """
        t = self.tables
        dst, src = t.dst, t.src
        d0, d1 = dst_range if dst_range is not None else (0, dst.n_blocks)
        s0, s1 = src_range if src_range is not None else (0, src.n_blocks)
        width = int(src.block_offset[s1] - src.block_offset[s0])
        inter = self.layout == "interleaved"
        nvec = x.shape[1] if inter else x.shape[0]
        want_x = (width, nvec) if inter else (nvec, width)
        want_y = (dst.dim, nvec) if inter else (nvec, dst.dim)
        if x.shape != want_x or y.shape != want_y:
            msg = f"x must be {want_x} and y {want_y}; got {x.shape}, {y.shape}"
            raise ValueError(msg)
        i0, i1 = int(self.items_first[d0]), int(self.items_first[d1])
        n_items = i1 - i0
        if n_items == 0:
            return
        item_block, item_k0 = self._item_block[i0:i1], self._item_k0[i0:i1]
        threads, min_blocks = self.threads, self.min_blocks
        if nvec > CAP_MAX_NVEC:  # more accumulators: a register cap would spill (3-6x slower)
            threads, min_blocks = THREADS, 1
        blocks = min((n_items * WARP + threads - 1) // threads, 65535 * 8)
        _kernel(nvec, inter, self.packed, self.lean, threads, min_blocks)(
            (blocks,),
            (threads,),
            (
                item_block,
                item_k0,
                np.int64(n_items),
                *self._args,
                cp.ascontiguousarray(x, dtype=cp.float64),
                y,
                np.int64(width),
                np.int64(dst.dim),
                np.int64(src.block_offset[s0]),
                np.int32(s0),
                np.int32(s1),
                np.int32(1 if accumulate else 0),
            ),
        )

    def matvec(self, X: DeviceArray) -> DeviceArray:
        """``X O^T`` for a device block of row vectors ``(nvec, D)`` (transposes if interleaved)."""
        X = cp.ascontiguousarray(cp.atleast_2d(cp.asarray(X, dtype=cp.float64)))
        if self.layout == "interleaved":
            Xt = cp.ascontiguousarray(X.T)
            Yt = cp.empty((self.tables.dst.dim, X.shape[0]))
            self.apply(Xt, Yt)
            return cp.ascontiguousarray(Yt.T)
        Y = cp.empty((X.shape[0], self.tables.dst.dim))
        self.apply(X, Y)
        return Y


def hamiltonian_cuda(
    basis: MBasis, inter: MSchemeInteraction, layout: str = "rows"
) -> CudaOperator:
    """H on ``basis`` on the GPU; ``setup_seconds`` covers the tables and the upload."""
    t0 = time.perf_counter()
    op = CudaOperator(hamiltonian_tables(basis, inter), layout)
    op.setup_seconds = time.perf_counter() - t0
    return op


def jplus_cuda(basis: MBasis, layout: str = "rows") -> CudaOperator:
    """J+ from ``basis`` (2M) to the basis with 2M + 2, on the GPU."""
    tables, _ = jplus_tables(basis)
    return CudaOperator(tables, layout)


def j_expectation_cuda(vectors: DeviceArray, jplus: CudaOperator, m2: int) -> NDArray[np.float64]:
    """``J`` of each row of device ``vectors`` from ``<J^2> = |J+ psi|^2 + M (M + 1)``."""
    M = m2 / 2.0
    jp = jplus.matvec(vectors)
    jsq = cp.einsum("ij,ij->i", jp, jp) / cp.einsum("ij,ij->i", vectors, vectors) + M * (M + 1)
    return j_from_jsq(np.asarray(cp.asnumpy(jsq), dtype=np.float64))
