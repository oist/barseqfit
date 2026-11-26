"""Resolver for the machine environment: store path and scratch directory.

Only the validation scripts call this; the library functions take their
parameters as arguments. The values come from the shell and a repo-root
``.env``: a variable set in the shell, even to an empty string, wins and the
``.env`` is not consulted for it; an empty value means unset. Nothing is read
at import time and ``os.environ`` is never modified.
"""

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

# Whitespace is space, tab, CR, FF and VT, as the R reader's \s.
_WS = " \t\r\f\v"
_ENV_LINE = (
    r"[ \t\r\f\v]*(?:export[ \t\r\f\v]+)?"
    r"([A-Za-z_][A-Za-z0-9_]*)[ \t\r\f\v]*=(.*)"
)
_STORE_VAR = "PAPKOU2023RUGGED_DATA"
_SCRATCH_VAR = "BARSEQFIT_SCRATCH"


@dataclass(frozen=True)
class MachineEnv:
    """Resolved store path (None when unset) and scratch directory.

    ``scratch_source`` is ``"shell"``, ``"env_file"`` or ``"default"``.
    """

    store: Path | None
    scratch: Path
    scratch_source: str


def read_env_file(path: Path) -> dict[str, str]:
    """Parse a ``.env`` file; ``{}`` when it does not exist.

    Lines are ``[export ]KEY=value``. Blank lines, ``#`` lines and lines
    without a valid key are skipped; the last assignment of a key wins.
    Whitespace around the key and value is trimmed, a trailing CR included.
    One matching pair of surrounding quotes is removed. Nothing else is
    interpreted: no escapes, inline comments or ``${VAR}`` expansion.
    """
    if not path.is_file():
        return {}
    values: dict[str, str] = {}
    for line in path.read_bytes().decode("utf-8").split("\n"):
        if line.lstrip(_WS).startswith("#"):
            continue
        match = re.fullmatch(_ENV_LINE, line)
        if match is None:
            continue
        value = match.group(2).strip(_WS)
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        values[match.group(1)] = value
    return values


def _lookup(
    key: str, environ: Mapping[str, str], file_values: Mapping[str, str]
) -> tuple[str, str]:
    """Return (value, source): a shell key wins even when empty."""
    if key in environ:
        return environ[key], "shell"
    return file_values.get(key, ""), "env_file"


def _as_path(value: str, repo_root: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else repo_root / path


def resolve_machine_env(
    repo_root: Path, environ: Mapping[str, str] | None = None
) -> MachineEnv:
    """Resolve the store and scratch from the shell and ``repo_root/.env``.

    ``environ`` defaults to ``os.environ``, read at call time. Relative paths
    resolve against ``repo_root``; the scratch default is
    ``repo_root/.scratch/papkou2023``. The store is not validated.
    """
    if environ is None:
        environ = os.environ
    file_values = read_env_file(repo_root / ".env")
    store, _ = _lookup(_STORE_VAR, environ, file_values)
    scratch, source = _lookup(_SCRATCH_VAR, environ, file_values)
    if not scratch:
        source = "default"
    return MachineEnv(
        store=_as_path(store, repo_root) if store else None,
        scratch=(
            _as_path(scratch, repo_root)
            if scratch
            else repo_root / ".scratch" / "papkou2023"
        ),
        scratch_source=source,
    )
