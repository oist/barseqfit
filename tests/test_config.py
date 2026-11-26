"""The committed run config loads and is no looser than the fixture defaults."""

from __future__ import annotations

from pathlib import Path

from validation.compare import Tolerance
from validation.config import load_run_config

REPO = Path(__file__).resolve().parent.parent


def test_fixture_defaults_are_no_stricter_than_the_run_config():
    """A smaller tolerance is stricter, so each fixture value is >= the Run config value."""
    config = load_run_config(REPO / "validation.toml")
    assert config["fit"]["chunk_size"] == 2000
    assert config["gate"]["published_pair_rows"] == 3518696
    fixture = Tolerance.fixture()
    run = Tolerance.from_run_config(config)
    assert fixture.abs_tol >= run.abs_tol
    assert fixture.rel_tol >= run.rel_tol
    assert fixture.noise_floor >= run.noise_floor
