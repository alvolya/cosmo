# GPU shell-model CI code: design and plan

Nov 19, 2026 · @Kyle · revised Nov 20, 2026 after a review of the CoSMo source and measured CoSMo baselines

A Python-orchestrated M-scheme CI code with a factorized on-the-fly matvec, one CUDA C kernel through CuPy and one Numba CPU kernel. The first production target is psdpf with Nħω truncation; the jj model spaces follow once psdpf runs in production. Strength functions are aimed at dimensions of a few times 10^9 on a single multi-GPU node, and eigenstates at about 1 to 2 x 10^9.

## Decisions and scope

The stack is fixed: Python orchestration, CuPy with a CUDA C `RawKernel` on GPU, Numba on CPU, mpi4py for more than one GPU. JAX is dropped.

| Topic | Decision |
| --- | --- |
| Basis | M-scheme, proton and neutron Slater determinants as 64-bit strings, partition-based truncation |
| Model spaces | sd for validation; psdpf with Nħω truncation as the first production target (CoSMo's `spsdpf` with FSU or WBP); jj spaces (jj44, jj45, jj55, jj56) after psdpf |
| Hamiltonian | Never stored at scale; factorized on-the-fly matvec. Stored CSR only in the reference path |
| Eigensolver | Block thick-restart Lanczos, FP64 |
| Observables | Eigenstates, one-body densities, spectroscopic factors, Lanczos strength functions |
| GPU targets | NVIDIA only: V100/V100S, A100, L40S, H200 on MSU ICER |
| CPU targets | AMD EPYC and Intel Xeon nodes, laptops |
| Interop | Read KSHELL `.snt`, OXBASH/NuShell `.sps` + `.int` (CoSMo's interaction library and BIGSTICK both use these), and CoSMo XML systems |

Deferred: Intel GPUs, AMD GPUs, non-Hermitian (continuum) solver, time propagation, J-scheme.

Not a goal: translating CoSMo source. The code is written from the published algorithms; CoSMo and KSHELL are used as numerical oracles.

## Layers

Only flat typed arrays cross the kernel boundary, so a backend is one function plus six vector operations.

| Layer | Owns | Runs in |
| --- | --- | --- |
| Model space | Orbits, interaction files, truncation rules | Python |
| Basis | Partitions, sub-blocks, block offsets | Numba |
| Operator | Jump tables and coupling tables for H, J², transition operators, a† and a | Numba |
| Kernel | `y = O x` for a block of vectors | Numba (CPU) or CUDA C (GPU) |
| Vector backend | `dot`, `axpy`, `scal`, `copy`, `norm`, allocation | NumPy or CuPy |
| Solver | Lanczos, restart, J filtering, strength functions | Python, backend-agnostic |
| Distribution | Block ownership per rank, source-block exchange | mpi4py |
| I/O | `.snt`, CoSMo XML, HDF5 checkpoints and results | Python |

The solver never imports CuPy or NumPy directly. It receives a backend object, which is what keeps a later Intel or Kokkos kernel a drop-in.

The operator layer is the one abstraction to get right. H, J², Eλ and Mλ operators, Gamow-Teller, and single-nucleon addition or removal are all instances of the same table format, differing only in which basis they map from and to.

## Data layouts

The many-body vector is a list of dense 2D blocks indexed by (proton state, neutron state), stored end to end in one flat array.

**Single-species level**

- An m-scheme orbit is (n, l, j, m). Each species has at most 64 of them, so a Slater determinant is one `uint64`. psdpf needs 38 bits per species, or 40 with the 0s shell included as in CoSMo's `spsdpf`.
- A partition is the occupation tuple over j-orbits. It fixes parity and the ħω count.
- A sub-block is (partition, M). Its determinants are stored contiguously and sorted, so a state is addressed as (sub-block id, local index), with the local index an `int32`.

**Many-body level**

- A block is a pair (proton sub-block, neutron sub-block) with Mp + Mn = M, the right total parity, and a partition pair the truncation allows.
- Amplitudes of a block are row-major with shape (dp, dn). The global index is `block_offset[B] + ip * dn + in`, an `int64`.
- Truncation only removes blocks. Inside a block the product structure is exact, which is what the factorization needs.

**Operator tables**

The proton-neutron part of H is a sum of products of one-body jumps:

```latex
H_{pn} = \sum_{q,r} V_{qr}\, P_q \otimes N_r, \qquad P_q = a^\dagger_{k'} a_k \ (\text{protons}), \quad N_r = a^\dagger_{l'} a_l \ (\text{neutrons})
```

- One-body jump table, per species and per pair (destination sub-block, source sub-block): CSR over the destination index. Each entry is (source index `int32`, operator label `int16` with the sign in the top bit).
- Coupling table `V[q, r]`: dense FP64, one per transfer class (Δm, Δparity). Proton class and neutron class must be conjugate.
- Same-species part (pp, nn, single-particle energies): CSR sparse matrix per pair of sub-blocks, values FP64, already contracted with the interaction. It acts on one index with the other as spectator.
- Task list, CSR over destination block: each task names a source block, a proton jump table, a neutron jump table and a coupling table.

**The pn gather, per output amplitude**

```latex
w_{B'}[i',j'] \mathrel{+}= \sum_{(q,i,s)\,\in\,J_p(i')} \ \sum_{(r,j,t)\,\in\,J_n(j')} s\,t\,V_{qr}\; v_B[i,j]
```

Each output amplitude is written by exactly one thread, so there are no races and no atomics. The cost is the number of nonzeros of the row, which is the floor for any on-the-fly scheme.

Table sizes scale with the single-species dimensions, not with D. That is why they fit on the device when H itself would need terabytes.

A rectangular operator (a†, charge exchange) uses the same tables with different source and destination bases. For a† the jump table of one species holds one-body creation entries and the other species is the identity.

## Kernel ABI

A kernel computes `y = O x` for a block of vectors from read-only flat arrays, and does nothing else.

```python
def apply(
    # destination and source bases (the same object for H and J^2)
    dst_block_offset,   # int64[nb_dst + 1]
    dst_block_psub,     # int32[nb_dst]   proton sub-block id
    dst_block_nsub,     # int32[nb_dst]   neutron sub-block id
    src_block_offset,   # int64[nb_src + 1]
    src_block_nsub_dim, # int32[nb_src]   dn of each source block
    # task list, CSR over destination block
    task_ptr,           # int64[nb_dst + 1]
    task_src_block,     # int32[nt]
    task_pjump,         # int32[nt]   proton jump-table id, -1 = identity
    task_njump,         # int32[nt]   neutron jump-table id, -1 = identity
    task_coupling,      # int32[nt]   coupling-table id, or same-species CSR id
    task_kind,          # int8[nt]    0 = pn, 1 = pp spectator n, 2 = nn spectator p
    # jump tables, all concatenated
    jump_row_ptr,       # int64[...]  per table, CSR over destination local index
    jump_table_start,   # int64[n_tables + 1]
    jump_src,           # int32[n_jumps]
    jump_op,            # int16[n_jumps]  label, sign in the top bit
    # coupling tables, concatenated dense
    coup_start,         # int64[n_coup + 1]
    coup_ncol,          # int32[n_coup]
    coup_data,          # float64[...]
    # same-species CSR matrices, concatenated
    ss_start, ss_row_ptr, ss_col, ss_val,
    # vectors
    x,                  # float64[nvec, D_chunk], C-contiguous: source amplitudes of blocks src_first..src_last
    y,                  # float64[nvec, D_dst],   C-contiguous
    dst_first, dst_last,# destination block range this call owns
    src_first, src_last,# source block range present in x; tasks whose source is outside are skipped
    accumulate,         # bool: False overwrites y on the owned range, True adds to it
): ...
```

**Contract**

1. With `accumulate=False`, `y` is overwritten on the owned block range; with `accumulate=True` the contributions are added. `y` is untouched outside the owned range, and every other array is read-only.
2. Gather only: within one call, each output amplitude is written by exactly one thread, for all `nvec` vectors. No atomics, no scatter.
3. No allocation inside the kernel and no Python objects across the boundary.
4. Deterministic: the summation order per output amplitude is fixed by the task and jump order, and by the order of the source chunks when a matvec is split into several calls. The chunk order is fixed by the distribution layer, so CPU and GPU agree to roundoff and reruns agree bitwise on one backend.
5. Index widths are fixed: `int64` for anything that scales with D or with the total jump count, `int32` for local indices and ids, `int16` for operator labels.
6. The destination range, the source range and `accumulate` are the only distribution hooks. A single-device matvec is one call with the full source range and `accumulate=False`. A distributed matvec is one call per source chunk: the first overwrites and the rest accumulate.

The source range and accumulate flag are in the ABI from stage 3, even though only stage 5 needs them, because retrofitting them changes every kernel and every test.

**Consequences worth knowing**

- Symmetry of H is not exploited. Computing each row in full costs up to 2x the flops of a symmetric scatter and buys race-freedom. On a GPU that trade is clearly right; on CPU it can be revisited behind the same ABI.
- The multi-vector axis is inside the kernel, so the index arithmetic and table reads are paid once per output amplitude and amortized over the block of vectors.
- A kernel may reorganize the pn sum internally, for example into gather, dense multiply by `V`, gather. That is an optimization experiment for stage 4, invisible to the solver.
- GPU baseline thread mapping: one warp per destination row i' of a block, with lanes across j'. The proton jump list Jp(i') is then read once per warp, and reads of `v_B[i, :]` are coalesced. A plain one-thread-per-amplitude mapping re-reads proton jumps for every j' and is load-imbalanced.
- Task kinds (pn, pp, nn) are launched as separate kernels, or the task list is sorted by kind, so warps do not diverge on `task_kind`.

**Vector backend**

`zeros(nvec, D)`, `dot`, `axpy`, `scal`, `copy`, `norm`, plus `to_host` and `from_host` for checkpointing. NumPy and CuPy both satisfy this without wrappers beyond dispatch.

## Solvers

One Lanczos implementation serves both eigenstates and strength functions; they differ in the pivot, the restart policy and the reorthogonalization policy.

**Eigenstates: block thick-restart Lanczos**

- Block size b of 4 to 16, basis capped at m blocks, restart keeps the wanted Ritz vectors plus a few guards.
- Full reorthogonalization against the retained basis. Convergence on the residual norm, not on energy differences.
- The basis lives wherever there is room: device first, then host RAM, then NVMe. Reorthogonalization streams through it in chunks.
- FP64 throughout. FP32 storage of the retained basis is a later memory option, to be tested against the residual criterion before it is trusted.

**Angular momentum**

- Run at M equal to the lowest J wanted. States of higher J appear in the same run.
- J² is an operator in the same table format, so the expectation value per Ritz vector costs one matvec.
- To target one J, filter the pivot block with a polynomial in J² that vanishes on the unwanted J values up to the maximum J of the space.
- Failure mode: roundoff reintroduces other J components over many iterations. The check is the J² expectation value at each restart, with re-filtering when it drifts.

**Strength functions**

For an operator O and initial state ψ0, take the pivot Oψ0 normalized, and run N plain Lanczos steps with no restart:

```latex
S(E) = \langle \psi_0 | O^\dagger O | \psi_0 \rangle \sum_{k=1}^{N} |y_{1k}|^2\, \delta(E - \theta_k)
```

- θ and y are the eigenvalues and eigenvectors of the N by N tridiagonal matrix. The first 2N moments of the true distribution are reproduced exactly in exact arithmetic.
- The same tridiagonal matrix gives the resolvent matrix element as a continued fraction, which is the form a self-energy calculation wants.
- Only three vectors are needed, so strength functions reach larger D than eigenstates on the same hardware.
- Loss of orthogonality produces duplicate Ritz values. Their summed weight stays correct, so a smoothed S(E) and low moments are unaffected; individual peak strengths at high N are not to be read off without reorthogonalization.
- The pivot has good J only if ψ0 has J = 0 or the pivot is projected. Otherwise the strength is a sum over final J and must be labelled as such.
- For addition and removal operators the pivot lives in the A+1 or A-1 basis, which is why operators are rectangular.

**Center of mass**

In psdpf, add the Lawson term to H. It separates spurious states cleanly only for a complete Nħω truncation. With any other truncation, report the center-of-mass expectation value per state instead of assuming separation. E1 strength is the case where this matters most.

## Memory budget and multi-GPU distribution

Vectors set the ceiling: one FP64 vector costs 8 GB per 10^9 basis states, and the tables are small by comparison.

The table gives the largest D with three working vectors on the device, using 75% of device memory and leaving the rest for tables. Hardware is from the [ICER cluster resources page](https://docs.icer.msu.edu/Cluster_Resources/).

| Node type | Device memory | Largest D, 3 vectors on device | Host RAM |
| --- | --- | --- | --- |
| 1x V100 or V100S | 32 GB | 1.0 x 10^9 | 178 to 367 GB |
| 1x A100 PCIe | 40 GB | 1.2 x 10^9 | 256 GB |
| 1x L40S | 48 GB | 1.5 x 10^9 | 768 GB |
| 1x A100 SXM | 80 GB | 2.5 x 10^9 | 512 GB |
| 1x H200 | 141 GB | 4.4 x 10^9 | 768 GB |
| 4x H200 | 564 GB | 1.8 x 10^10 | 768 GB |
| 8x H200 | 1128 GB | 3.5 x 10^10 | 768 GB |
| CPU only, EPYC 9654 | none | 2.3 TB holds 28 vectors at 10^10 | 768 GB to 2.3 TB |

Three vectors is the strength-function case. The multi-GPU rows assume the chunked exchange below, so that each GPU holds only its own share of x and y plus one receive buffer. Without chunking, each GPU would need most of x and the per-GPU limit would not grow with the GPU count.

Eigenstates need the retained Krylov basis as well. With block size 8 and about 100 retained vectors, D = 10^10 needs about 8 TB, which exceeds both device memory (1.1 TB on 8x H200) and host RAM (768 GB). The practical single-node eigenstate limit is therefore about 1 to 2 x 10^9 with the basis on device. Above that the basis streams from NVMe on every reorthogonalization, or the run goes multi-node.

On the 8x H200 nodes the host has less memory than the devices combined. Spilling to host does not work there; the basis must stay small and on device, or go to NVMe.

The first psdpf targets are much smaller than these ceilings. The 1ħω negative-parity spaces of the sd-shell nuclei benchmarked below run from 2.7 x 10^4 (20Ne) to 3.8 x 10^6 (28Si) and fit on one GPU. Multi-GPU becomes necessary at 2ħω and above, or for the heavier jj spaces.

**Distribution**

- One MPI rank per GPU. Each rank owns a contiguous range of destination blocks, ordered by proton sub-block and balanced by nonzero count, not by dimension.
- Which terms leave a rank's blocks:

| Term | Proton sub-block | Neutron sub-block | Source local? |
| --- | --- | --- | --- |
| Diagonal and single-particle energies | same | same | yes |
| nn two-body | same | changes | yes |
| pp two-body | changes | same | often not |
| pn | changes | changes | often not |

- So a rank may need a large fraction of x. It never needs it all at once. The matvec runs as one kernel call per source chunk:

```
rank r:
  y_r  = kernel(x_r,         src = own blocks,    accumulate=False)
  for each other rank s, in a fixed order:
      receive the blocks of x_s that rank r's tasks reference
      y_r += kernel(chunk_s, src = those blocks,  accumulate=True)
      discard chunk_s
```

- Memory per rank: about D/n for each of x and y, plus one receive buffer. Communication per matvec: up to (n-1)/n of x per rank, about 6 GB at D = 10^9 on 4 GPUs. Over NVLink that is tens of milliseconds against a matvec of seconds; over InfiniBand it is about 10x slower. Single-node multi-GPU should be cheap, and multi-node is the real risk.
- The remote fraction, meaning the bytes of x each rank needs from the others, is computed from the task list in stage 2 for real psdpf cases and splits of 2, 4 and 8 ranks, before any GPU code exists. The ħω structure of psdpf may make it much smaller than the worst case.
- Within one node the exchange is device-to-device. Across nodes it goes over InfiniBand through MPI.

## Staged plan

Each stage is the regression oracle for the next, so no stage starts before the previous gate passes.

| Stage | Deliverable | Gate |
| --- | --- | --- |
| 0. Baseline | CoSMo built with `-O3 -march=x86-64-v3 -DNDEBUG` in `bin-O3/`. H construction time, stored-H matvec time and nonzero count (`bench-psdpf/csb_matvec_bench.cpp`), and Davidson time, each recorded separately for spsdpf FSU9 1ħω cases. KSHELL (OpenMP, on-the-fly) time per matvec on the same cases: it stores no H, so it is the stronger CPU baseline | Partly done Nov 20, 2026 (see the speedup section). Still to add: one pf-shell case, 48Cr, and the KSHELL timings for 24Mg/28Si sd, 48Cr and the psdpf 1ħω cases |
| 0b. Smoke test | A trivial gather `RawKernel` under CuPy on a V100 and an H200 node, with the ICER CUDA modules | Compiles and matches NumPy on both |
| 1. CPU reference | Basis with partitions and truncations, stored CSR, block thick-restart Lanczos, J², `.snt` and CoSMo XML readers | Energies within 1e-6 MeV of CoSMo (its printed precision, with float32 inputs) and within 1e-8 MeV of KSHELL on identical `.snt` input; dimensions equal to both; J exact. Progress Nov 22, 2026 on branch `gpu-rewrite` of github.com/beykyle/cosmo (`cosmo_gpu/`): Nħω and occupation truncation, partition-first basis, mixed parity and the `.snt` reader are done. CoSMo: dimensions match through 28Si 1ħω (3,834,452); energies and J match for sd USDB (including odd-odd 22Na and 26Al) and psdpf FSU9 (33Al 0ħω/1ħω, 20Ne and 22Na 1ħω), once inputs are rounded to float32 as CoSMo's `std::stof` parser does. KSHELL (`tests/test_kshell.py`): dimensions, energies to 1e-8 MeV and J match for sd USDB (20Ne, 22Na, 24Mg, 25Mg, 26Al, 28Si), pf KB3G (44Ti, 46Ti) and psdpf FSU9 (33Al 0ħω⁺/1ħω⁻, 20Ne and 22Na 1ħω⁻); KSHELL's own partitions give the same truncated dimensions as ours and CoSMo's. Every coupled matrix element of this code's isospin-to-pn conversion matches KSHELL's `nushell2snt.py` for USDB, KB3G and FSU9. 48Cr moved to stage 3: a stored H at D ≈ 2 x 10^6 needs about 100 GB. Still to do: CoSMo XML reader |
| 2. Operator layer | Jump tables, rectangular operators, one-body densities, spectroscopic factors, strength functions; remote-fraction count from the task list | Non-energy-weighted sum rules to 1e-10; agreement with CoSMo `XSHLSF` and `XSHLEMB`; B(E2), B(M1) and orbit occupations also against KSHELL `transit.exe` and its `<p Nj>`/`<n Nj>` output, with the same effective charges, g-factors and oscillator length (record KSHELL's `hw_type` convention for b); remote fraction reported for 24Mg and 28Si 1ħω at 2, 4 and 8 ranks |
| 3. On-the-fly CPU kernel | Numba kernel behind the ABI, including the source range and accumulate flag | Matches the stage 1 matvec to 1e-13 relative on random vectors, both as one call and split into source chunks; 48Cr kb3g (D ≈ 2 x 10^6, too big for the stored oracle) energies match KSHELL to 1e-8 MeV; a solver comparison with KSHELL's thick-restart block Lanczos at the same `n_block`, `n_restart_vec` and `max_lanc_vec` (matvec counts and convergence history on 28Si and 48Cr); table setup time and matvec time compared with the stage 0 CoSMo and KSHELL numbers |
| 4. GPU kernel | CUDA C kernel behind the same ABI; solver vector operations on device | Matches stage 3 to 1e-13; throughput measured on A100 and H200 for 24Mg and 28Si 1ħω and reported as speedup over both CoSMo and KSHELL on a CPU node of the same generation |
| 5. Multi-GPU | Block ownership, chunked source exchange, checkpointing | A psdpf case with D near 10^9 (2ħω) on one 4x H200 node, with energies checked against an MPI KSHELL run (`-DMPI`, OpenMPI 4.1.5 module) that also gives the node-count baseline; exchange cost measured against the matvec |
| 6. psdpf production | Lawson term, center-of-mass expectation values, psdpf runs | A published psdpf spectrum with FSU or WBP reproduced, plus CoSMo agreement on 33Al and 24Mg 1ħω and KSHELL agreement on 24Mg and 28Si 1ħω (per-A FSU9 `.snt` from `tests/data/kshell/make_snt.py`) |
| 7. jj spaces | Same-species tables generated on the fly if they do not fit; jj production runs | A published jj44 or jj55 spectrum reproduced; jun45 or jj44 energies match KSHELL on the same `.snt` (KSHELL ships both) |
| Deferred | Non-Hermitian solver, time propagation, Intel and AMD GPUs |  |

**Validation set, fixed in stage 1 and reused at every gate**

- 20Ne, 24Mg, 28Si with USDB: full spectra by dense diagonalization where the dimension allows.
- One odd-A and one odd-odd sd-shell nucleus, to exercise half-integer M and non-zero ground-state J.
- 48Cr in the pf shell (M-scheme dimension about 2 x 10^6): the first case that stresses the on-the-fly path, gated against KSHELL in stage 3.
- spsdpf with FSU9 at 1ħω, negative parity: 20Ne (D = 26,600), 22Na (218,924), 24Mg (978,948) and 28Si (3,834,452). These are the stage 0 benchmark cases and the first psdpf regression cases. 33Al at 0ħω and 1ħω (119 and 23,842 states) has known CoSMo energies.
- Random-vector matvec comparisons, which catch phase and index errors that converged energies can hide.
- KSHELL on identical `.snt` input (`tests/test_kshell.py`): sd USDB 20Ne, 22Na, 24Mg, 25Mg, 26Al, 28Si; pf KB3G 44Ti, 46Ti; psdpf FSU9 33Al 0ħω⁺ and 1ħω⁻, 20Ne and 22Na 1ħω⁻. Dimensions, energies to 1e-8 MeV and J. The `.snt` files are converted from the repo's OXBASH files by KSHELL's `nushell2snt.py`, and every coupled matrix element is also compared with this code's own isospin-to-pn conversion.
- BIGSTICK ships a 25Mg USDB reference (`examples/sd/output/mg25ref.res`, ground state −94.40128 MeV) that is used directly. See "Prior art and code-code validation".

- CoSMo quirks found while validating: it parses interaction files with `std::stof`, which is single precision. Its `CoreQuanta` (`include/sm/SingleParticleStates.cxx`) sorts an empty range, so the Nħω reference configuration is the first N m-orbits in file order. That is harmless for spaces listed in ascending quanta, such as spsdpf, but differs from the true minimum otherwise. Some of its Davidson states stop unconverged at the iteration cap (24Mg USDB 3⁺₁, eigenvector error 0.21).

**Order of risk**

Stages 1 and 2 are routine. Stage 3 is where phase conventions and jump-table bugs surface. Stages 4 and 5 carry the performance risk.

## Speedup against CoSMo: measured baseline and estimates

The CoSMo numbers below are measured. The numbers for the new code are estimates until stages 3 and 4 replace them.

**What CoSMo does today (public repository, checked Nov 20, 2026)**

- Neither diagonalizer uses an on-the-fly matvec. Both build H once from the operator form (`PPMatrix`) and store it. `DEXHVJT` stores row-sparse rows (`PPMatrix2SparseMatrix`, `include/xsm/DEXHV.cxx:62`). `DEXHVJTcsb` stores compressed sparse blocks: the upper triangle only, at 12 bytes per stored nonzero (`PPMatrix2CSBBlk`, `include/xsm/DEXHVcsb.cxx:66`). Below `-l` (default 1000 states), both build a fully dense matrix and diagonalize it with Eigen.
- The choice depends on dimension only, not on the model space.
- On-the-fly evaluation is used in two places: to build the stored H, and in observables that never store a matrix (`OverlapNH` in `XSHLEMB` and `XSHLSF`). The search-and-atomic multiply `PPMatrix::TimesVectorSymmetric` exists but no diagonalizer calls it.
- Building H costs a generic operator application per basis state, a search of the full basis for each resulting ket, and a sorted container. That is the expensive step, and it is paid once per run.
- The basis is not factorized into proton and neutron parts.
- Davidson (`include/kdav/fsdavidson.cxx`) runs one matvec for each unconverged vector in every iteration, up to the block size. With 10 states, 400 iterations cost about 4,150 matvecs.
- `build.sh` and the Makefile pass `-std=c++17 -fopenmp` with no optimization level. OpenMP only: no MPI, no GPU.

**Measured: default build against optimized build**

Setup: spsdpf with FSU9, 1ħω, negative parity. `DEXHVJTcsb -n 10 -l 1` with the default 400-iteration cap. 16 OpenMP threads on an AMD EPYC 9654 node (amd24, one job per node type, nodes shared). The optimized build is `-O3 -march=x86-64-v3 -DNDEBUG` with g++ 14.3, in `~/cosmo/bin-O3`. Matvec times come from `bench-psdpf/csb_matvec_bench.cpp`, which builds H exactly as `DEXHVJTcsb` does and times 20 stored multiplies.

| | 22Na | 24Mg | 28Si (optimized only) |
| --- | --- | --- | --- |
| Dimension | 218,924 | 978,948 | 3,834,452 |
| Nonzeros of H (full matrix) | 7.6 x 10^7 (349 per row) | 4.5 x 10^8 (456 per row) | 2.2 x 10^9 (576 per row) |
| Stored H, CSB | 0.46 GB | 2.7 GB | 13.3 GB (15 GB peak) |
| Build H, default / optimized | 335 s / 159 s (2.1x) | 1961 s / 501 s (3.9x) | 2984 s |
| Davidson, 400 iterations, default / optimized | 268 s / 49 s (5.5x) | 1392 s / 305 s (4.6x) | about 1,500 s, projected from 4,150 matvecs |
| One stored matvec, default / optimized | 19.9 ms / 11.8 ms (1.7x) | 103 ms / 73 ms (1.4x) | 365 ms |
| Whole run, default / optimized | 643 s / 220 s (2.9x) | 3559 s / 857 s (4.2x) | about 75 min, projected |

What the numbers say:

- **The build flags are worth 3 to 4x end to end.** Most of the gain is outside the stored multiply: the stored multiply itself gains only 1.4 to 1.7x, while the Davidson vector operations, which take about 70% of Davidson time in the default build, become negligible. Adding `-O3` to `build.sh` and the Makefile is worth doing for current CoSMo users regardless of this project.
- **In the optimized build, Davidson time is almost all matvecs.** About 4,150 matvecs at 11.8 ms gives 49 s for 22Na, and at 73 ms gives 305 s for 24Mg, which matches the measured Davidson times.
- **Building H costs more than all 400 Davidson iterations.** It takes 1.2 to 2.4 µs per nonzero on 16 threads, against 0.16 ns per nonzero for a stored multiply, so one build costs as much as about 7,000 to 15,000 multiplies.
- **Costs per nonzero are flat with D.** From 22Na to 28Si (a factor of 17.5 in D and 29 in nonzeros), building H costs 1.2 to 2.4 µs per nonzero and a stored matvec costs 0.155 to 0.165 ns per nonzero. CoSMo's early progress estimate of 12 hours for the 28Si build was wrong; it took 50 minutes. Both costs scale with the nonzeros of H, which grow faster than D (349 to 576 per row).
- **Memory, not time, is CoSMo's limit.** 28Si needs 13.3 GB of stored H at D = 3.8 x 10^6. At the same 12 bytes per stored nonzero and about 600 or more per row, D = 10^8 needs about 400 GB and D = 10^9 about 4 TB or more.
- **Convergence caveat.** At 400 iterations only the lowest one to three states reach eigenvector errors below 1e-5, and the highest requested states stay at 1e-2 to 5e-2. Both builds do the same 400 iterations, so the comparison is fair, but production runs need more iterations. That shifts the total toward matvecs.

**Estimated: new code against optimized CoSMo**

| Comparison | Estimate | Basis |
| --- | --- | --- |
| Setup: jump tables against CoSMo's H build | 10 to 100x or more | CoSMo pays about 1 to 2 µs per nonzero of H. Jump tables scale with the single-species dimensions, not with the nonzeros of H |
| CPU matvec, same cores | 0.5 to 2x | The stored CSB multiply is already efficient, at 0.16 ns per nonzero on 16 cores. The factorized kernel replaces matrix streaming with table reads and index arithmetic |
| One H200 against CoSMo on 16 to 32 cores, per matvec | 10 to 30x | CoSMo's stored multiply streams about 40 GB/s of matrix. An H200 has 4.8 TB/s, and an irregular gather typically reaches 10 to 30% of that, while reading fewer bytes per nonzero than a stored matrix |
| 24Mg 1ħω end to end, one H200 against optimized CoSMo (857 s) | 10 to 40x | Setup of 10 to 60 s, plus 4,150 matvecs at 2 to 9 ms each |
| Same, against default CoSMo (3559 s) | 40 to 150x | As above |
| Memory | No stored H | 24Mg stores 2.7 GB of H. The new code holds vectors and tables only, which is what reaches 10^9 states |

What would change these estimates:

- Stage 3 measures the CPU kernel against the 0.16 ns per nonzero baseline. If the kernel is slower than about 0.3 ns per nonzero on the same cores, the CPU path is a memory win but a speed loss.
- Many small blocks in heavily truncated psdpf spaces lower GPU efficiency. Stage 4 measures this on 24Mg and 28Si 1ħω.
- No speedup is claimed from the eigensolver. Davidson and block Lanczos need comparable numbers of matvecs.

**Reproducing**

`bench-psdpf/bench.sb` runs `DEXHVJTcsb` and the matvec benchmark for given cases and builds, and `bench-psdpf/matvec_only.sb` runs the matvec benchmark alone. The system XMLs and bases are in `bench-psdpf/`, and the per-run logs are in `bench-psdpf/runs/`.

## Prior art and code-code validation

Checked Nov 20, 2026 from the public repositories, read only, nothing built. Neither KSHELL nor BIGSTICK has a public GPU implementation: a search for OpenACC, OpenMP target offload, CUDA and HIP across all branches, the full history and the active forks found nothing. The GPU work that does exist in nuclear CI is MFDn (OpenACC, papers listed in the sources) and Shimizu's Monte Carlo shell model on H100 (a 2025 talk), not KSHELL.

| | KSHELL | BIGSTICK |
| --- | --- | --- |
| Public source | No official repository. Official tarballs (v4, 2023-07) on the KSHELL site; maintained fork GaffaSnobb/kshell (v4, last commit 2025-10) | cwjsdsu/BigstickPublick, version 8.0.0 (master tip 2026-07) |
| GPU | None in code | None in code; the Perlmutter build target is CPU only |
| License | Unclear: GPLv3 in the CPC paper, "academic use only" on the site, no license file | MIT (LLNS and SDSU) |
| Basis | Blocks of (proton partition, neutron partition, Mp); each block a dense proton x neutron matrix | Proton and neutron sectors by (Jz, parity, W); each proton sector lists its conjugate neutron sectors |
| pn matvec | One-body jump lists per species, then `lwf(p,n) += V(ik,jl) rwf(p',n')` (`bridge_partitions.F90`) | Double loop over proton x neutron one-body jumps, element `hmatpn(cpnpair + dpnpair)` (`bapplyhlibOMP.f90`) |
| Multi-vector | Up to 64 vectors per pass (`bp_block.F90`) | Block-Lanczos variants (`bapplyhlibOMP_block.f90`) |
| Write conflicts | Threads over partition tasks | Owner-computes: each thread owns a contiguous output range and skips writes outside it, no atomics |
| MPI distribution | 2-D rank grid; the input vector moves around a ring (`shift_mpi_*`), partial results reduced | Vector split into fragments of contiguous sectors (default 500 M states); needs nproc ≥ nfragments² |
| Precision | FP64 | FP32 Lanczos vectors and pn matrix elements by default |
| Interaction input | `.snt` (pn formalism); NuShellX converter `nushell2snt.py` | OXBASH/NuShell `.sps` + `.int` (iso or pn), MFDn TBME files; no `.snt` |
| Nħω truncation | Yes, through `gen_partition.py` into the `.ptn` file; Lawson term via `beta_cm` | Yes, W truncation with a maximum excitation |
| Observables | Lanczos strength (E1, E2, M1, GT and others), transitions, spectroscopic factors, OBTD/TBTD | Lanczos strength, one-body densities, occupations; spectroscopic factors need the companion code Rhodium |
| psdpf interactions shipped | None (closest: ysox for psd, sdpf-m, sdpf-mu) | None |
| Reference outputs | Sparse (one 20Ne USDA level list in an older fork) | `examples/sd/output/` (25Mg USDB, densities, strengths) |

**What this means for the design**

- Both codes use the factorization in this plan: proton x neutron blocks, one-body jump lists per species, and a pn term built as an outer product. That supports the core design; the new parts are the GPU kernel, the Python solver layer and the observables interface.
- BIGSTICK's owner-computes thread ranges are the CPU form of contract item 2 (gather, no atomics). Its "opbundles" (operator type, sector pair, jump range, split by measured cost) are the same idea as the task list here, and its cost-based splitting is worth copying for GPU load balance.
- KSHELL's ring shift of the input vector is the closest published precedent for the chunked exchange in the distribution section. BIGSTICK's fragment scheme is the alternative, and its nproc ≥ nfragments² requirement is a cost to avoid. Stage 5 should compare the ring order with the per-rank chunk order described above.
- KSHELL's 64-vector blocking confirms that putting the multi-vector axis inside the kernel pays off.

**Using them as oracles**

- BIGSTICK reads OXBASH-format `.int` files, the same format as CoSMo's interaction library (`interactions/spsdpf/fsu9/*.int`). It is likely, but not yet verified, that FSU and WBP run in BIGSTICK with only a new `.sps` file giving the W (ħω) value of each orbit.
- KSHELL needs `.snt` files. `tests/data/kshell/make_snt.py` converts the repo's OXBASH files with KSHELL's own `nushell2snt.py`. FSU9 sums six files with different mass scaling while a `.snt` has one mass factor, so it gets one file per mass number. `nushell2snt.py` reads exactly the declared number of TBMEs while CoSMo reads them all and warns: `fp/kb3g.int` declares 194 and lists 195, and the converter would silently drop p1/2² J=1 T=0 (−0.674 MeV). The script passes it a corrected copy; KSHELL's shipped `kb3g.snt` includes the element.
- KSHELL's log prints energies to 5–6 decimals, but the `.wav` header stores them as float64 (`neig`, `mtotal`, eigenvalues, 2J as a raw stream), which is what the 1e-8 MeV comparison uses. Built outside the repo with `tools/kshell/build_kshell.sh` (pinned commit 41ea03a, conda gfortran 14.3, OpenMP).
- The Nħω conventions must be checked before comparing truncated spaces: KSHELL counts relative to the lowest configuration, BIGSTICK uses W relative to the orbit assignments, and CoSMo uses its own `hwrejection`. A matching dimension for the same nucleus and truncation is the first check.
- BIGSTICK defaults to FP32 vectors and matrix elements, so comparisons with it are good to about 1e-5 relative unless it is built for FP64. CoSMo and KSHELL are FP64 and support the 1e-8 MeV gates.
- Both build with gfortran or Intel Fortran plus LAPACK (BIGSTICK's MPI build needs `mpi_f08` bindings). KSHELL is built on ICER (Nov 22, 2026); BIGSTICK is not, and its FP32 default makes it a weaker oracle than KSHELL.

## Risks and open questions

The largest technical risk is the cost of exchanging source blocks between GPUs; the largest non-technical one is licensing.

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| Source-block exchange dominates the multi-GPU matvec | Poor scaling past one GPU | Order blocks by proton sub-block; chunked exchange with an accumulating kernel; count the remote fraction in stage 2 and measure the exchange in stage 5 before multi-node work |
| Each rank needs most of x | Per-GPU memory does not grow with the GPU count | Chunked exchange: each rank holds its own x and y plus one receive buffer |
| Same-species two-body tables grow too large in jj spaces | Tables no longer fit on device | Generate pp and nn jumps on the fly from bit operations inside the kernel, behind the same ABI |
| Table construction in Python is slow for 10^5 or more partitions | Setup time rivals solve time | Numba for all construction; cache tables in HDF5 keyed by model space and truncation |
| V100 nodes need an older CUDA than the newer cards | Two CuPy builds to maintain | Stage 0b settles it; unverified until then |
| Mixed CPU instruction sets across ICER nodes | Illegal-instruction crashes from cached native code | Numba compiles per node; key any on-disk cache by CPU model |
| J contamination over long Lanczos runs | Wrong J labels, slow convergence | J² check at every restart, re-filter on drift |
| Center-of-mass contamination in truncated psdpf | Spurious states in spectra and E1 strength | Lawson term plus reported center-of-mass expectation values |
| CoSMo is GPL v3 | Code derived from its source would carry the GPL | Implement from published algorithms; decide the project license before any porting. Not legal advice |
| KSHELL license is unclear (GPLv3 per the CPC paper, "academic use only" per its site, no license file in any repository) | Reusing its code is legally unclear | Use KSHELL only as a numerical oracle and as a published algorithm description. BIGSTICK is MIT and can be read freely, with attribution if anything is adapted |
| Interaction files carry their own terms (the FSU files state distribution terms in their headers) | Redistributing converted files may be restricted | Convert interaction files locally for validation; do not ship them without checking their terms |

**Open questions**

- [x] First production target: psdpf with Nħω truncation. jj spaces come after (stage 7).
- [ ] Which jj space and nucleus comes first after psdpf? It fixes the largest D and the table sizes for stage 7.
- [ ] Which psdpf truncation: 0+1 ħω only, or higher? 1ħω fits on one GPU for the sd-shell nuclei measured; 2ħω is what makes stage 5 necessary.
- [x] Is there a CoSMo version that runs Davidson with an on-the-fly matvec? No: neither the public repository nor any other known version does, so the stored path is the baseline.
- [x] Do KSHELL or BIGSTICK have a public GPU path? No (checked Nov 20, 2026; see "Prior art and code-code validation"). The reach argument for this code stands.
- [ ] MFDn has published OpenACC GPU work (arXiv:2109.00485, arXiv:2110.10765), but it is a no-core code. Can it run valence-space psdpf with an empirical interaction such as FSU, and is it available to us? If both, it is the one GPU comparison point.
- [ ] Is CoSMo XML output compatibility needed, or only reading its inputs?
- [ ] What license should the new code carry?
- [ ] How much H200 time is realistically available? It decides whether stage 5 targets 4 or 8 GPUs.

## Sources

- [CoSMo repository](https://github.com/alvolya/cosmo): source read for the matvec and build analysis.
- [ICER cluster resources](https://docs.icer.msu.edu/Cluster_Resources/): node types and memory.
- [Thick-restart block Lanczos method for large-scale shell-model calculations](https://arxiv.org/pdf/1902.02064): KSHELL on-the-fly matvec and partition structure.
- [BIGSTICK: a flexible configuration-interaction shell-model code](https://arxiv.org/pdf/1801.08432): factorized on-the-fly algorithm.
- [CuPy v14 release notes](https://www.preferred.jp/en/news/pr20260219): current CuPy status.
- [KSHELL site](https://sites.google.com/alumni.tsukuba.ac.jp/kshell-nuclear/) and [GaffaSnobb/kshell](https://github.com/GaffaSnobb/kshell): KSHELL v4 source, read for the prior-art survey.
- [BigstickPublick](https://github.com/cwjsdsu/BigstickPublick): BIGSTICK 8.0.0 source, read for the prior-art survey.
- [MFDn LOBPCG with OpenACC](https://arxiv.org/abs/2109.00485) and [OpenMP versus OpenACC offload in MFDn](https://arxiv.org/abs/2110.10765): existing GPU work in nuclear CI.
