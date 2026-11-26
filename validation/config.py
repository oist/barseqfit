"""Loader for the committed run config ``validation.toml``.

Only the validation scripts call this; the library functions take their
parameters as arguments.
"""

from __future__ import annotations

import tomllib
from pathlib import Path


def load_run_config(path: Path | str) -> dict:
    """Return the parsed run config at ``path`` as a nested dict."""
    with Path(path).open("rb") as handle:
        return tomllib.load(handle)
