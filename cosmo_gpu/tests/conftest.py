"""Shared pytest fixtures."""

from __future__ import annotations

import pytest

from cosmo_gpu.modelspace import ModelSpace, load_model_space


@pytest.fixture(scope="session")
def sd_space() -> ModelSpace:
    return load_model_space("sd")
