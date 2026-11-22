# cosmo_gpu

An M-scheme shell-model CI code written in Python, with a factorized on-the-fly matvec: Numba on
CPU and CuPy on GPU. It lives inside a fork of [CoSMo](https://github.com/alvolya/cosmo). It reads
CoSMo's interaction catalogue (`../interactions`) and uses CoSMo as its numerical oracle. The
design and staged plan are in [`docs/design-plan.md`](docs/design-plan.md).

## Status

Stage 1 (CPU reference), first slice:

- valence spaces and OXBASH isospin interactions from CoSMo's `int.xml`, including mass scaling;
- M-scheme basis as proton × neutron blocks of (partition, M) sub-blocks;
- stored-CSR Hamiltonian and J+ built from bitstrings (the reference path, not for production);
- block thick-restart Lanczos behind a vector-backend interface.

USDB energies and spins for 18O, 19O, 20Ne, 24Mg, 25Mg and 28Si agree with CoSMo to its printed
precision (`tests/test_cosmo_usdb.py`).

## Install and run

```bash
pip install -e "cosmo_gpu[dev]"
python -m cosmo_gpu levels 24Mg sd usdb -n 10
```

Interaction files are found through `$SMINTPATH` (as set by CoSMo's `cosmovar.sh`), otherwise
in `$COSMO/interactions`, otherwise in this repository's `interactions/`.

## Checks

Run all of these from `cosmo_gpu/`. Every commit must pass them, and CI runs the same set:

```bash
ruff check .          # lint (rules in pyproject.toml)
ruff format --check . # formatting
mypy                  # strict type checking
pytest                # fast tests; `pytest -m slow` for the 28Si check (use a Slurm job)
```

`pre-commit install` at the repository root runs ruff and mypy on staged `cosmo_gpu` files.

Conventions: angular momenta are stored doubled (`j2`, `m2`), so they are integers. Arrays carry
their dtype in the type (`NDArray[np.uint64]`, `NDArray[np.int64]`, ...), because index widths
are part of the kernel ABI. Numba kernels live in `*_kernels.py`, and typed wrappers carry the
public API.
