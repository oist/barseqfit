"""Adapter layer shared by the full-data WT and pair gates.

The gate rule itself (``GateRule``, ``compare_chain``) lives in
``validation.gate``; this module holds what is script-side: reading
``validation.toml`` ( call ``load_gate_config`` from ``main``, nothing
is read at import time), the Export checks, the fit timing / fallback capture,
the Markdown / CSV / timing writers. The
store and scratch directory come from ``validation.env.resolve_machine_env``.
"""

from __future__ import annotations

import hashlib
import json
import logging
import sys
import time
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

import polars as pl

REPO = Path(__file__).resolve().parent.parent
# The non-installed `validation` package sits at the repo root. Every gate
# script imports this module first, so this is the one place that puts it on the path.
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from validation.config import load_run_config
from validation.gate import ChainVerdict, GateRule, summarise

TIME_UNIT = "generations"


def load_gate_config() -> dict:
    """Read the committed Run config; scripts call this from main()."""
    return load_run_config(REPO / "validation.toml")


def refuse(
    prog: str, message: str, hint: str = "Run `just export-papkou` first."
) -> None:
    print(f"{prog}: {message}\n{hint}", file=sys.stderr)
    sys.exit(2)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_export(
    prog: str,
    scratch: Path,
    store: Path | None,
    needed: list[str],
    r_table: str,
    r_timing: str,
    r_hint: str,
) -> None:
    """Refuse a missing or stale export; check R freshness when its table exists.

    ``store`` is the resolved store path; None skips the same-store check.
    """
    missing = [n for n in needed if not (scratch / n).exists()]
    manifest_path = scratch / "export_manifest.json"
    if missing or not manifest_path.exists():
        refuse(
            prog,
            f"export incomplete in {scratch} (missing {missing or manifest_path.name})",
        )
    manifest = json.loads(manifest_path.read_text())
    if store is not None and store.resolve() != Path(manifest["store"]).resolve():
        refuse(
            prog, "export was made from a different store than PAPKOU2023RUGGED_DATA"
        )
    for entry in manifest["store_files"]:
        path = Path(entry["path"])
        if not path.exists():
            refuse(prog, f"store file {path} no longer exists")
        stamp = datetime.fromisoformat(entry["mtime"]).timestamp()
        if abs(path.stat().st_mtime - stamp) > 1.0 or md5(path) != entry["md5"]:
            refuse(prog, f"store file {path} changed since the export (stale manifest)")
    if (scratch / r_table).exists():
        check_r_fresh(prog, scratch, manifest, r_timing, r_hint)


def check_r_fresh(
    prog: str, scratch: Path, manifest: dict, r_timing: str, hint: str
) -> None:
    """Refuse an R table not fitted from the current export."""
    path = scratch / r_timing
    if not path.exists():
        refuse(prog, f"{path.name} missing next to the R table", hint)
    timing = json.loads(path.read_text())
    want_md5 = [e["md5"] for e in manifest["store_files"]]
    if timing.get("export_timestamp") != manifest["timestamp"]:
        refuse(
            prog, "R table is stale: export timestamp differs from the manifest", hint
        )
    if timing.get("store_md5") != want_md5:
        refuse(prog, "R table is stale: store md5s differ from the manifest", hint)


class _Count(logging.Handler):
    """Capture the scalar-oracle fallback count the vectorised engines log."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.fallbacks: int | None = None

    def emit(self, record: logging.LogRecord) -> None:
        if "re-fitted by the scalar oracle" in record.msg:
            self.fallbacks = int(record.args[0])


def fit_and_count(
    fit_callable: Callable[[], pl.DataFrame],
) -> tuple[pl.DataFrame, float, int | None]:
    """Run a fit, returning (table, seconds, scalar-oracle fallback count)."""
    counter = _Count()
    log = logging.getLogger("barseqfit._comparison")
    log.addHandler(counter)
    old_level = log.level
    log.setLevel(logging.DEBUG)
    try:
        t0 = time.perf_counter()
        table = fit_callable()
        seconds = time.perf_counter() - t0
    finally:
        log.setLevel(old_level)
        log.removeHandler(counter)
    return table, seconds, counter.fallbacks


def render_verdict_md(
    *,
    title: str,
    header_lines: list[str],
    sv_stats: dict,
    lines_sets: list[str],
    before_only_lines: list[str] | None,
    status_ls: list[str],
    status_disagree: list[str],
    warn_lines: list[str],
    pseudo_lines: list[str],
    chain: ChainVerdict,
    failures: list[str],
    scalar_ok: bool,
    rule: GateRule,
    extra_sections: list[tuple[str, list[str]]] | None = None,
) -> tuple[list[str], bool]:
    """Assemble the verdict Markdown; returns (lines, passed).

    ``header_lines`` are the script-specific bullets after the shared tolerance
    bullet. ``before_only_lines`` is None for gates without that section.
    ``extra_sections`` are (heading, lines) pairs placed before the verdict.
    """
    tol = rule.tolerance
    md = [
        title,
        "",
        f"- Tolerance: abs {tol.abs_tol:g}, plus rel {tol.rel_tol:g} when |expected| < 1 (fixture rule; relative bound floored at {tol.noise_floor:g} for m); p with both values below {tol.p_underflow:g} (IEEE underflow) counts as equal.",
        *header_lines,
        "",
        "## Scalar vs vectorised (seeded sample, seed 0)",
        f"- Entries differing by more than {rule.scalar_tol:g}: {sv_stats['bad']}; status mismatches: {sv_stats['status_bad']}; iteration mismatches: {sv_stats['iter_bad']}; max abs difference {sv_stats['max']:.3e}",
        "",
        "## Set differences",
        *lines_sets,
        *chain.r_lines,
    ]
    if before_only_lines is not None:
        md += ["", "## Before-only genotypes", *before_only_lines]
    md += ["", "## Python vs published", *summarise(chain.res_pub, "pub")]
    if chain.res_r is not None:
        md += ["", "## Python vs R", *summarise(chain.res_r, "r")]
    md += ["", "## Declared rule outcome", *chain.rule_lines]
    if chain.tie_lines:
        md += ["", "## Threshold ties", *chain.tie_lines]
    md += [
        "",
        "## Status counts (Python)",
        *status_ls,
        "",
        "## Status disagreements",
        *status_disagree,
        "",
        "## Warning messages (Python)",
        *warn_lines,
        "",
        "## Pseudocount check",
        *pseudo_lines,
    ]
    for heading, lines in extra_sections or []:
        md += ["", f"## {heading}", *lines]
    passed = not failures and scalar_ok
    md += [
        "",
        f"## Verdict: {'PASS' if passed else 'FAIL'}",
        f"- Asserted failures: {len(failures)}",
        *[f"  - {f}" for f in failures],
        f"- Scalar/vectorised agreement: {'ok' if scalar_ok else 'FAILED'}",
    ]
    return md, passed


def write_verdict(
    scratch: Path,
    stem: str,
    md_lines: list[str],
    verdict_csv: pl.DataFrame,
    timing: dict,
    timing_name: str,
) -> None:
    """Write ``{stem}.md``, ``{stem}.csv`` and the timing JSON; print the Markdown."""
    (scratch / f"{stem}.md").write_text("\n".join(md_lines) + "\n")
    verdict_csv.write_csv(scratch / f"{stem}.csv")
    (scratch / timing_name).write_text(json.dumps(timing, indent=2) + "\n")
    print("\n".join(md_lines))
