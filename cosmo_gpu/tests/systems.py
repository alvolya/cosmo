"""Test helper: shell-model systems built once per test session."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

import scipy.sparse as sp

from cosmo_gpu.basis import MBasis, build_mbasis, build_species_basis
from cosmo_gpu.interaction import (
    CoupledInteraction,
    MSchemeInteraction,
    read_interaction,
    to_mscheme,
)
from cosmo_gpu.modelspace import ModelSpace, load_hamiltonian_spec, load_model_space, parse_nucleus
from cosmo_gpu.reference_csr import hamiltonian_csr


@dataclass(frozen=True)
class System:
    space: ModelSpace
    coupled: CoupledInteraction
    inter: MSchemeInteraction
    basis: MBasis
    H: sp.csr_matrix


@cache
def build_system(
    nucleus: str, space_name: str = "sd", ham: str = "usdb", cosmo_float32: bool = False
) -> System:
    space = load_model_space(space_name)
    Z, N = parse_nucleus(nucleus)
    zv, nv = Z - space.core_z, N - space.core_n
    spec = load_hamiltonian_spec(ham, space_name)
    coupled = read_interaction(spec, space, zv + nv, cosmo_float32)
    inter = to_mscheme(coupled)
    P, Nb = build_species_basis(space, zv), build_species_basis(space, nv)
    basis = build_mbasis(space, P, Nb, (zv + nv) % 2)
    return System(space, coupled, inter, basis, hamiltonian_csr(basis, inter))
