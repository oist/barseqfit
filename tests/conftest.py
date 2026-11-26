"""Shared pytest fixtures."""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def data_dir() -> Path:
    """Repo-root tests/data, shared by the Python and R tests."""
    return REPO_ROOT / "tests" / "data"
