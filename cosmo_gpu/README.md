# cosmo_gpu

An M-scheme shell-model CI code written in Python, with a factorized on-the-fly matvec: Numba on
CPU and CuPy on GPU. It lives inside a fork of [CoSMo](https://github.com/alvolya/cosmo). It reads
CoSMo's interaction catalogue (`../interactions`) and uses CoSMo as its numerical oracle. The
design and staged plan are in [`docs/design-plan.md`](docs/design-plan.md).

## Status

Stage 1 (CPU reference) in progress:

- valence spaces and OXBASH isospin interactions (multi-file, scaled, mass-dependent) from
  CoSMo's `int.xml`;
- Nħω and orbit-occupation truncations with CoSMo's semantics, and mixed-parity spaces;
- M-scheme basis as proton × neutron blocks of (partition, M) sub-blocks, built partition by
  partition so a truncated space is never enumerated in full;
- stored-CSR Hamiltonian and J+ built from bitstrings (the reference path, not for production);
- block thick-restart Lanczos behind a vector-backend interface.

Agreement with CoSMo, to its printed precision (`tests/test_cosmo_*.py`):
- sd USDB: 18O, 19O, 20Ne, 22Na, 24Mg, 25Mg, 26Al and 28Si;
- psdpf (`spsdpf`) FSU9: 33Al at 0ħω and 1ħω, 20Ne at 1ħω, and 22Na at 1ħω (slow test);
- basis dimensions up to 28Si at 1ħω (3,834,452).

## Install and run

The package is a [uv](https://docs.astral.sh/uv/) project. Python ≥ 3.12 is required, and
`uv.lock` pins every dependency. From `cosmo_gpu/`:

```bash
uv sync                                     # creates .venv with the package and dev tools
uv run cosmo-gpu levels 24Mg sd usdb -n 10  # or: uv run python -m cosmo_gpu levels ...
uv run cosmo-gpu levels 33Al spsdpf fsu9 --parity -1 --hw 1:1 -n 10
uv sync --extra gpu                         # add CuPy (also: mpi, io)
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
