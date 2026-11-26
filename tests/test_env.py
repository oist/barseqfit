"""The machine-environment resolver: .env grammar, shell-wins rule, import purity."""

from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

from validation import env as _env
from validation.env import MachineEnv, read_env_file, resolve_machine_env

STORE = "PAPKOU2023RUGGED_DATA"
SCRATCH = "BARSEQFIT_SCRATCH"

# Joined with "\n" and a final "\n", so only the last line ends in "\r\n".
GRAMMAR_LINES = [
    "# comment line",
    "   # indented comment",
    "",
    "PLAIN=value",
    "export EXPORTED=exp",
    "  SPACED  =  spaced value  ",
    'DQ="double quoted"',
    "SQ='single quoted'",
    'QEMPTY=""',
    'QPAD="  padded  "',
    "MISMATCH=\"mismatched'",
    "HASH=a # not a comment",
    "EMPTY=",
    "DUP=first",
    "DUP=second",
    "INTERP=${HOME}/x",
    "no equals sign here",
    "1BAD=invalid key",
    "CRLF=crlf\r",
]
GRAMMAR_PARSED = {
    "PLAIN": "value",
    "EXPORTED": "exp",
    "SPACED": "spaced value",
    "DQ": "double quoted",
    "SQ": "single quoted",
    "QEMPTY": "",
    "QPAD": "  padded  ",
    "MISMATCH": "\"mismatched'",
    "HASH": "a # not a comment",
    "EMPTY": "",
    "DUP": "second",
    "INTERP": "${HOME}/x",
    "CRLF": "crlf",
}


def _write_env(root: Path, text: str) -> None:
    (root / ".env").write_text(text)


def test_read_env_file_grammar(tmp_path):
    path = tmp_path / ".env"
    path.write_bytes(("\n".join(GRAMMAR_LINES) + "\n").encode())
    assert read_env_file(path) == GRAMMAR_PARSED


def test_read_env_file_missing_file_is_empty(tmp_path):
    assert read_env_file(tmp_path / ".env") == {}


def test_read_env_file_empty_file_and_no_final_newline(tmp_path):
    path = tmp_path / ".env"
    path.write_bytes(b"")
    assert read_env_file(path) == {}
    path.write_bytes(b"A=1\nB=2")
    assert read_env_file(path) == {"A": "1", "B": "2"}


def test_read_env_file_trims_space_tab_cr_ff_and_vt(tmp_path):
    path = tmp_path / ".env"
    path.write_bytes(b"A=value\r \n\x0bB=\t x\x0c\n")
    assert read_env_file(path) == {"A": "value", "B": "x"}


def test_shell_value_wins_over_env_file(tmp_path):
    _write_env(tmp_path, f"{STORE}=/file/store\n")
    menv = resolve_machine_env(tmp_path, {STORE: "/shell/store"})
    assert menv.store == tmp_path / "/shell/store"


def test_store_from_env_file_when_shell_lacks_it(tmp_path):
    _write_env(tmp_path, f"{STORE}=/file/store\n")
    assert resolve_machine_env(tmp_path, {}).store == tmp_path / "/file/store"


def test_store_expands_the_home_directory_and_resolves_relative_paths(tmp_path):
    menv = resolve_machine_env(tmp_path, {STORE: "~/store"})
    assert menv.store == Path("~/store").expanduser()
    menv = resolve_machine_env(tmp_path, {STORE: "data/store"})
    assert menv.store == tmp_path / "data/store"


def test_empty_shell_value_wins_and_means_unset(tmp_path):
    _write_env(tmp_path, f"{STORE}=/file/store\n")
    assert resolve_machine_env(tmp_path, {STORE: ""}).store is None


def test_scratch_from_env_file_when_shell_lacks_it(tmp_path):
    _write_env(tmp_path, f"{SCRATCH}={tmp_path / 'file_scratch'}\n")
    menv = resolve_machine_env(tmp_path, {})
    assert menv == MachineEnv(None, tmp_path / "file_scratch", "env_file")


def test_empty_shell_scratch_gives_the_default(tmp_path):
    _write_env(tmp_path, f"{SCRATCH}={tmp_path / 'file_scratch'}\n")
    menv = resolve_machine_env(tmp_path, {SCRATCH: ""})
    assert menv.scratch == tmp_path / ".scratch" / "papkou2023"
    assert menv.scratch_source == "default"


def test_empty_env_file_scratch_gives_the_default(tmp_path):
    _write_env(tmp_path, f"{SCRATCH}=\n")
    menv = resolve_machine_env(tmp_path, {})
    assert menv.scratch == tmp_path / ".scratch" / "papkou2023"
    assert menv.scratch_source == "default"


def test_scratch_defaults_when_nothing_is_set(tmp_path):
    menv = resolve_machine_env(tmp_path, {})
    assert menv == MachineEnv(None, tmp_path / ".scratch" / "papkou2023", "default")


def test_relative_env_file_scratch_resolves_against_the_repo_root(tmp_path):
    _write_env(tmp_path, f"{SCRATCH}=rel/dir\n")
    menv = resolve_machine_env(tmp_path, {})
    assert menv.scratch == tmp_path / "rel/dir"
    assert menv.scratch_source == "env_file"


def test_shell_scratch_expands_the_home_directory(tmp_path):
    menv = resolve_machine_env(tmp_path, {SCRATCH: "~/x"})
    assert menv.scratch == Path("~/x").expanduser()
    assert menv.scratch_source == "shell"


def test_store_is_none_when_unset(tmp_path):
    assert resolve_machine_env(tmp_path, {}).store is None


def test_resolver_reads_but_does_not_change_os_environ(tmp_path, monkeypatch):
    _write_env(tmp_path, f"{STORE}=/file/store\n{SCRATCH}=from_file\n")
    monkeypatch.setenv(STORE, "/shell/store")
    monkeypatch.delenv(SCRATCH, raising=False)
    before = dict(os.environ)
    menv = resolve_machine_env(tmp_path)
    assert dict(os.environ) == before
    # "/shell/store" has no drive on Windows, so the resolver joins it to the
    # root there; tmp_path / "/x" does the same and is just "/x" on POSIX.
    expected = MachineEnv(tmp_path / "/shell/store", tmp_path / "from_file", "env_file")
    assert menv == expected


def test_resolver_does_not_change_the_environ_it_is_given(tmp_path):
    _write_env(tmp_path, f"{SCRATCH}=from_file\n")
    environ = {STORE: "/shell/store"}
    resolve_machine_env(tmp_path, environ)
    assert environ == {STORE: "/shell/store"}


def test_env_module_runs_nothing_at_import_time():
    """The module body holds a docstring, stdlib imports, constants, defs and classes."""
    tree = ast.parse(Path(_env.__file__).read_text())
    for node in tree.body:
        if isinstance(node, ast.Expr | ast.Assign | ast.AnnAssign):
            assert isinstance(node.value, ast.Constant), ast.dump(node)
        elif isinstance(node, ast.Import):
            assert all(
                a.name.split(".")[0] in sys.stdlib_module_names for a in node.names
            )
        elif isinstance(node, ast.ImportFrom):
            assert node.module.split(".")[0] in sys.stdlib_module_names
        else:
            assert isinstance(node, ast.FunctionDef | ast.ClassDef), ast.dump(node)


def test_importing_barseqfit_does_not_import_validation():
    code = "import sys, barseqfit; assert 'validation' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)
