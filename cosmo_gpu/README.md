# cosmo_gpu

An M-scheme shell-model CI code written in Python, with a factorized on-the-fly matvec: Numba on
CPU and CuPy on GPU. It lives inside a fork of [CoSMo](https://github.com/alvolya/cosmo). It reads
CoSMo's interaction catalogue (`../interactions`) and KSHELL `.snt` files, and uses CoSMo and
KSHELL as numerical oracles. The design and staged plan are in
[`docs/design-plan.md`](docs/design-plan.md).

## Status

Stages 1 (CPU reference), 0b (GPU smoke test), 2 (operator layer), 3 (on-the-fly CPU kernel) and 4 (GPU kernel) are complete; stage 5 (multi-GPU) is next. The code has:

- valence spaces and OXBASH isospin interactions (multi-file, scaled, mass-dependent) from
  CoSMo's `int.xml`, KSHELL `.snt` proton-neutron interactions (`--snt`), and CoSMo system
  XML files (nucleus, space, Hamiltonian and truncations of a CoSMo job);
- Nħω and orbit-occupation truncations with CoSMo's semantics, and mixed-parity spaces;
- M-scheme basis as proton × neutron blocks of (partition, M) sub-blocks, built partition by
  partition so a truncated space is never enumerated in full;
- stored-CSR Hamiltonian and J+ built from bitstrings (the reference path, not for production);
- block thick-restart Lanczos behind a vector-backend interface (NumPy and CuPy backends);
- operator tables for H and J+ (jump, coupling, same-species and task tables) and the
  remote-fraction count for multi-GPU planning;
- an on-the-fly matvec (Numba) that applies H and J+ from those tables without storing H, with
  the kernel ABI's source/destination ranges and accumulate flag (`--matvec otf`, automatic above
  2 x 10^5 states); 48Cr in the pf shell (D = 1,963,461) matches KSHELL to 1e-8 MeV;
- the same matvec as one CUDA kernel with the Lanczos solve on the GPU (`--matvec gpu`, needs the
  `gpu` extra): 28Si spsdpf 1ħω (D = 3.8 x 10^6, 10 states) in 18 s on one H200, against 607 s
  for our CPU kernel and 661 s for KSHELL on 16 cores;
- E2, M1, Gamow-Teller, occupations, one-body densities and spectroscopic factors with reduced
  matrix elements, and Lanczos strength functions.

Agreement with CoSMo, to its printed precision (`tests/test_cosmo_*.py`):
- sd USDB: 18O, 19O, 20Ne, 22Na, 24Mg, 25Mg, 26Al and 28Si;
- psdpf (`spsdpf`) FSU9: 33Al at 0ħω and 1ħω, 20Ne at 1ħω, and 22Na at 1ħω (slow test);
- basis dimensions up to 28Si at 1ħω (3,834,452).

Agreement with KSHELL on the same `.snt` input, to 1e-8 MeV (`tests/test_kshell.py`): sd USDB,
pf KB3G (44Ti, 46Ti) and psdpf FSU9 cases; B(E2), B(M1) and occupations against `transit.exe`.
See the KSHELL section below.

Operators: sum rules to 1e-10 (`tests/test_operators.py`, `tests/test_strength.py`); E2, M1,
occupations and spectroscopic factors against CoSMo `XSHLEMB`, `XSHLAO` and `XSHLSF`
(`tests/test_cosmo_ops.py`).

## Install and run

The package is a [uv](https://docs.astral.sh/uv/) project. Python ≥ 3.12 is required, and
`uv.lock` pins every dependency. From `cosmo_gpu/`:

```bash
uv sync                                     # creates .venv with the package and dev tools
uv run cosmo-gpu levels 24Mg sd usdb -n 10  # or: uv run python -m cosmo_gpu levels ...
uv run cosmo-gpu levels 33Al spsdpf fsu9 --parity -1 --hw 1:1 -n 10
uv run cosmo-gpu levels 48Cr fp kb3g --matvec otf -n 10   # H never stored (use a Slurm job)
uv run cosmo-gpu levels 48Cr fp kb3g --matvec gpu -n 10   # on a GPU node, with the gpu extra
uv sync --extra gpu                         # add CuPy (also: mpi, io)
uv run cosmo-gpu levels 24Mg --snt tests/data/kshell/snt/usdb_from_int.snt  # a KSHELL file
uv run cosmo-gpu levels tests/data/cosmo_fsu9/33Al_1hw-.xml   # rerun a CoSMo system file
uv run cosmo-gpu transitions 24Mg sd usdb --op E2 --charges 1.5,0.5 -n 6
uv run cosmo-gpu sf 24Mg 25Mg sd usdb -n 3                    # adding one neutron
uv run cosmo-gpu strength 20Ne sd usdb --op M1 --steps 100
```

Interaction files are found through `$SMINTPATH` (as set by CoSMo's `cosmovar.sh`), otherwise
in `$COSMO/interactions`, otherwise in this repository's `interactions/`.

## Checks

Run all of these from `cosmo_gpu/`. Every commit must pass them, and CI runs the same set:

```bash
uv run ruff check .          # lint (rules in pyproject.toml)
uv run ruff format --check . # formatting
uv run mypy                  # strict type checking
uv run pytest                # fast tests; `uv run pytest -m slow` for 28Si (use a Slurm job)
```

`uvx pre-commit install` at the repository root runs ruff and mypy on staged `cosmo_gpu` files.

Conventions: angular momenta are stored doubled (`j2`, `m2`), so they are integers. Arrays carry
their dtype in the type (`NDArray[np.uint64]`, `NDArray[np.int64]`, ...), because index widths
are part of the kernel ABI. Numba kernels live in `*_kernels.py`, and typed wrappers carry the
public API.

## KSHELL references

KSHELL (fork GaffaSnobb/kshell, pinned commit) is the second oracle: double precision, the same
`.snt` input, and the same thick-restart block Lanczos family. Its license is unclear, so it is
built outside this repository and only its outputs are committed.

```bash
tools/kshell/build_kshell.sh                                 # clone + build into ~/kshell
uv run python tests/data/kshell/make_snt.py                  # repo .int -> snt/ (KSHELL's nushell2snt)
sbatch tests/data/kshell/kshell.sb                           # run kshell.exe -> kshell_reference.json
```

`KSHELL_DIR` overrides `~/kshell`. `run_kshell.py` writes partitions with KSHELL's
`gen_partition`, so truncated dimensions are KSHELL's own, and reads eigenvalues at full
precision from the `.wav` header (the log prints only 5–6 decimals).

## GPU smoke test

`tools/gpu/gpu_smoke.sb` runs `python -m cosmo_gpu.gpu_smoke` on one GPU per job (see the script
header for the five ICER GPU types). It needs a separate environment with the `gpu` extra,
created once on the login node: `UV_PROJECT_ENVIRONMENT=.venv-gpu uv sync --locked --extra gpu`.
Results are in `tools/gpu/smoke_results.jsonl`.
