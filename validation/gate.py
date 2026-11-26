"""The gate rule as a pure module (``GateRule``, ``compare_chain``).

This module takes a ``GateRule`` value built by the caller, performs no I/O,
loads no config at import time, logs nothing and prints nothing. It lives in the
non-installed ``validation`` package at the repository root, next to
``compare.py``, so the installed ``barseqfit`` package holds only fitting code.
The scripts under ``scripts/`` are the only readers of the Run config; they
load the Export, fit, build the ``GateRule`` and render.

``keys`` names the row-identifying columns: ``["sequence_id"]`` for WT,
``["from", "to"]`` for pairs. Callers pass already-filtered frames (the WT gate
drops the reference row before calling).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import polars as pl

from validation.compare import Tolerance

COLUMNS = ("m", "se", "p")


@dataclass(frozen=True)
class GateRule:
    """Tolerance, tie caps and scalar-vs-vectorised tolerance as one value."""

    tolerance: Tolerance
    tie_fraction: float
    tie_se_abs: float
    tie_p_abs: float
    scalar_tol: float
    max_tie_rows: int = 50

    @classmethod
    def from_run_config(cls, config: Mapping) -> GateRule:
        gate = config["gate"]
        return cls(
            tolerance=Tolerance.from_run_config(config),
            tie_fraction=gate["tie_fraction"],
            tie_se_abs=gate["tie_se_abs"],
            tie_p_abs=gate["tie_p_abs"],
            scalar_tol=gate["scalar_tol"],
        )


@dataclass(frozen=True)
class ChainVerdict:
    """Residual frames, report lines, the per-row verdict table and failures."""

    res_pub: pl.DataFrame
    res_r: pl.DataFrame | None
    verdict_csv: pl.DataFrame
    r_lines: list[str]
    rule_lines: list[str]
    tie_lines: list[str]
    failures: list[str]

    @property
    def passed(self) -> bool:
        return not self.failures


def keyset(frame: pl.DataFrame, keys: list[str]) -> set:
    """Set of row keys: scalars for one key column, tuples for several."""
    if len(keys) == 1:
        return set(frame[keys[0]].to_list())
    return set(zip(*(frame[k].to_list() for k in keys), strict=True))


def residuals(
    fitted: pl.DataFrame,
    other: pl.DataFrame,
    tag: str,
    keys: list[str],
    tolerance: Tolerance,
) -> pl.DataFrame:
    """Per-row abs residual and out-of-tolerance flag per column vs ``other``."""
    joined = fitted.select(*keys, *COLUMNS).join(
        other.select(*keys, *COLUMNS),
        on=keys,
        how="inner",
        suffix="_o",
    )
    out = joined.select(*keys)
    for c in COLUMNS:
        a = joined[c].fill_null(np.nan).to_numpy()
        e = joined[f"{c}_o"].fill_null(np.nan).to_numpy()
        res = tolerance.compare(a, e, c)
        if c == "p":
            out = out.with_columns(pl.Series(f"p_uf_{tag}", res.underflow))
        out = out.with_columns(
            pl.Series(f"{c}_abs_{tag}", res.abs_residual),
            pl.Series(f"{c}_rel_{tag}", res.rel_residual),
            pl.Series(f"{c}_out_{tag}", res.outside),
        )
    return out


def summarise(res: pl.DataFrame, tag: str) -> list[str]:
    lines = ["| column | max abs | max rel | rows outside |", "|---|---|---|---|"]
    for c in COLUMNS:
        absr = res[f"{c}_abs_{tag}"].drop_nans().drop_nulls()
        rel = res[f"{c}_rel_{tag}"].to_numpy()
        rel = rel[np.isfinite(rel)]
        lines.append(
            f"| {c} | {absr.max() if len(absr) else float('nan'):.3e} "
            f"| {rel.max() if len(rel) else float('nan'):.3e} "
            f"| {int(res[f'{c}_out_{tag}'].sum())} |"
        )
    return lines


def subset_residual_lines(
    fitted: pl.DataFrame,
    expected: pl.DataFrame,
    keys: list[str],
    subset: pl.DataFrame,
    tolerance: Tolerance,
) -> list[str]:
    """Residuals vs ``expected`` on the rows of ``fitted`` whose keys are in ``subset``.

    ``subset`` holds the key columns. Reporting only; nothing is asserted.
    """
    wanted = subset.select(keys).unique()
    res = residuals(
        fitted.join(wanted, on=keys, how="semi"),
        expected.join(wanted, on=keys, how="semi"),
        "sub",
        keys,
        tolerance,
    )
    return [
        f"- Rows compared with the published table: {res.height}",
        *summarise(res, "sub"),
    ]


def set_diff_lines(name_a: str, a: set, name_b: str, b: set) -> list[str]:
    return [
        f"- {name_a} not in {name_b}: {len(a - b)}"
        + (f" (e.g. {sorted(a - b)[:3]})" if a - b else ""),
        f"- {name_b} not in {name_a}: {len(b - a)}"
        + (f" (e.g. {sorted(b - a)[:3]})" if b - a else ""),
    ]


def warning_counts(fitted: pl.DataFrame) -> pl.DataFrame:
    return (
        fitted.select(pl.col("diagnostics").struct.field("warning").alias("warning"))
        .filter(pl.col("warning").is_not_null())
        .group_by("warning")
        .len()
        .sort("len", descending=True)
    )


def scalar_vs_vectorised(sv: pl.DataFrame, scalar_tol: float) -> dict:
    """Compare a scalar fit joined to the vectorised one (suffix ``_v``)."""
    status_bad = int((sv["status"] != sv["status_v"]).sum())
    iter_bad = int(
        (
            sv["diagnostics"].struct.field("iterations").fill_null(-1)
            != sv["diagnostics_v"].struct.field("iterations").fill_null(-1)
        ).sum()
    )
    bad = 0
    max_diff = 0.0
    for c in COLUMNS:
        a = sv[c].fill_null(np.nan).to_numpy()
        e = sv[f"{c}_v"].fill_null(np.nan).to_numpy()
        both_nan = np.isnan(a) & np.isnan(e)
        with np.errstate(invalid="ignore"):
            diff = np.where(both_nan, 0.0, np.abs(a - e))
        bad += int((~(diff <= scalar_tol)).sum())
        finite = diff[np.isfinite(diff)]
        max_diff = max(max_diff, float(finite.max()) if len(finite) else 0.0)
    return {
        "bad": bad,
        "status_bad": status_bad,
        "iter_bad": iter_bad,
        "max": max_diff,
        "ok": bad == 0 and status_bad == 0 and iter_bad == 0,
    }


def status_lines(vec: pl.DataFrame, r_table: pl.DataFrame | None, keys: list[str]):
    counts = vec.group_by("status").len().sort("status")
    lines = [f"- {r['status']}: {r['len']}" for r in counts.iter_rows(named=True)]
    if r_table is not None and "status" in r_table.columns:
        dis = vec.select(*keys, "status").join(
            r_table.select(*keys, pl.col("status").alias("status_r")), on=keys
        )
        dis = dis.filter(pl.col("status") != pl.col("status_r"))
        disagree = [f"- Python vs R status disagreements: {dis.height}"]
        disagree += [
            f"  - {' / '.join(str(r[k]) for k in keys)}: python {r['status']}, R {r['status_r']}"
            for r in dis.head(10).iter_rows(named=True)
        ]
    else:
        n_failed = vec.filter(pl.col("status") == "failed").height
        disagree = [
            f"- Python vs published: published table carries no status; {n_failed} Python rows failed."
        ]
    return lines, disagree


def compare_chain(
    vec: pl.DataFrame,
    expected: pl.DataFrame,
    r_table: pl.DataFrame | None,
    keys: list[str],
    rule: GateRule,
) -> ChainVerdict:
    """Apply the declared tolerance / tie rule.

    Returns a ``ChainVerdict``: residual frames, report lines, the per-row
    verdict table and the list of asserted failures.
    """
    tolerance = rule.tolerance
    failures: list[str] = []
    res_pub = residuals(vec, expected, "pub", keys, tolerance)
    verdict_csv = vec.select(*keys, "status", "m", "se", "p").join(
        res_pub, on=keys, how="left"
    )
    r_lines: list[str] = []
    rule_lines: list[str] = []
    tie_lines: list[str] = []
    res_r = None
    n_pub_out = int(res_pub["m_out_pub"].sum())
    if n_pub_out:
        failures.append(f"m outside tolerance vs published: {n_pub_out}")
    if r_table is None:
        n_bad = int(sum(res_pub[f"{c}_out_pub"].sum() for c in ("se", "p")))
        rule_lines.append("- R table absent: se, p asserted strictly vs published.")
        if n_bad:
            failures.append(f"se/p outside tolerance vs published: {n_bad}")
        return ChainVerdict(
            res_pub, None, verdict_csv, r_lines, rule_lines, tie_lines, failures
        )

    fitted_ids, pub_ids, r_ids = (
        keyset(vec, keys),
        keyset(expected, keys),
        keyset(r_table, keys),
    )
    r_lines += set_diff_lines("R", r_ids, "fitted", fitted_ids)
    r_lines += set_diff_lines("R", r_ids, "published", pub_ids)
    if r_ids != pub_ids:
        failures.append(
            f"R ids != published ids ({len(r_ids - pub_ids)} only R, "
            f"{len(pub_ids - r_ids)} only published)"
        )
    res_r = residuals(vec, r_table, "r", keys, tolerance)
    res_rp = residuals(r_table, expected, "rp", keys, tolerance)
    verdict_csv = verdict_csv.join(res_r, on=keys, how="left")
    if "status" in r_table.columns:
        verdict_csv = verdict_csv.join(
            r_table.select(*keys, pl.col("status").alias("status_r")),
            on=keys,
            how="left",
        )
    iters = vec.select(
        *keys,
        pl.col("diagnostics").struct.field("iterations").fill_null(-1).alias("it_py"),
    ).join(
        r_table.select(
            *keys,
            pl.col("diagnostics")
            .str.extract(r"iterations=(\d+)", 1)
            .cast(pl.Int64)
            .fill_null(-1)
            .alias("it_r"),
        ),
        on=keys,
    )
    tab = (
        res_r.join(res_pub, on=keys)
        .join(res_rp, on=keys, how="left")
        .join(iters, on=keys)
        .with_columns(
            (pl.col("it_py") != pl.col("it_r")).alias("tie"),
            (pl.col("se_out_rp") | pl.col("p_out_rp")).fill_null(False).alias("rp_out"),
            (pl.col("se_out_pub") | pl.col("p_out_pub")).alias("pub_out"),
            (pl.col("se_out_r") | pl.col("p_out_r")).alias("r_out"),
        )
    )
    verdict_csv = verdict_csv.join(
        tab.select(*keys, "it_py", "it_r", "tie"), on=keys, how="left"
    )
    n_rows = tab.height
    ties = tab.filter(pl.col("tie")).sort("se_abs_r", descending=True)
    n_ties = ties.height
    tie_cap = int(rule.tie_fraction * n_rows)
    tie_se_max = float(ties["se_abs_r"].max()) if n_ties else 0.0
    tie_p_max = float(ties["p_abs_r"].max()) if n_ties else 0.0

    # m vs R: every row.
    n_m_r = int(tab["m_out_r"].sum())
    if n_m_r:
        failures.append(f"m outside tolerance vs R: {n_m_r}")
    # se, p vs R on matched rows.
    matched = tab.filter(~pl.col("tie"))
    n_match_bad = int(matched["r_out"].sum())
    if n_match_bad:
        failures.append(f"se/p outside tolerance vs R on matched rows: {n_match_bad}")
    if n_ties > tie_cap:
        failures.append(f"threshold ties {n_ties} exceed cap {tie_cap}")
    if tie_se_max > rule.tie_se_abs:
        failures.append(f"tie se residual {tie_se_max:.3e} exceeds {rule.tie_se_abs:g}")
    if tie_p_max > rule.tie_p_abs:
        failures.append(f"tie p residual {tie_p_max:.3e} exceeds {rule.tie_p_abs:g}")
    # se, p vs published on rows that are neither ties nor R-vs-published outliers.
    gated = tab.filter(~pl.col("tie") & ~pl.col("rp_out"))
    n_gated_bad = int(gated["pub_out"].sum())
    if n_gated_bad:
        failures.append(
            f"se/p outside tolerance vs published on gated rows: {n_gated_bad}"
        )
    n_pub_se_out = int(tab["pub_out"].sum())
    n_rp_out = int(tab["rp_out"].sum())
    if n_rp_out > tie_cap:
        failures.append(f"R-vs-published outliers {n_rp_out} exceed cap {tie_cap}")
    if n_pub_se_out > n_rp_out + n_ties:
        failures.append(
            f"Python-vs-published outliers {n_pub_se_out} exceed R-vs-published {n_rp_out} + ties {n_ties}"
        )
    rest = tab.filter(pl.col("pub_out") & (pl.col("tie") | pl.col("rp_out")))

    rule_lines += [
        f"- rows compared: {n_rows}",
        f"- m vs published / vs R: outside tolerance {n_pub_out} / {n_m_r} (asserted 0)",
        f"- matched rows (equal iterations): {matched.height}; se/p outside tolerance vs R: {n_match_bad} (asserted 0)",
        f"- threshold ties (iterations differ): {n_ties} (cap {tie_cap} = {rule.tie_fraction:.2%} of rows); max se residual {tie_se_max:.3e} (cap {rule.tie_se_abs:g}); max p residual {tie_p_max:.3e} (cap {rule.tie_p_abs:g})",
        f"- R vs published: rows with se/p outside tolerance: {n_rp_out} (cap {tie_cap}); max se residual {float(tab['se_abs_rp'].max()):.3e}",
        f"- Python vs published: rows with se/p outside tolerance: {n_pub_se_out} (asserted <= {n_rp_out} + {n_ties})",
        f"- se/p vs published asserted on {gated.height} rows (not ties, not R-vs-published outliers): outside tolerance {n_gated_bad} (asserted 0)",
        f"- p underflow (both sides < {tolerance.p_underflow:g}, treated as equal): {int(tab['p_uf_r'].sum())} rows vs R, {int(tab['p_uf_pub'].sum())} vs published",
        f"- Python-vs-published outliers excused by a tie or an R outlier (reported): {rest.height}",
    ]
    header = " / ".join(keys)
    tie_lines = [
        f"| {header} | R iter | Py iter | se residual vs R |",
        "|---|---|---|---|",
    ] + [
        f"| {' / '.join(str(r[k]) for k in keys)} | {r['it_r']} | {r['it_py']} | {r['se_abs_r']:.3e} |"
        for r in ties.head(rule.max_tie_rows).iter_rows(named=True)
    ]
    if n_ties > rule.max_tie_rows:
        tie_lines.append("")
    tie_lines.append(
        f"Showing {min(n_ties, rule.max_tie_rows)} of {n_ties} threshold ties (largest se residual first)."
    )
    return ChainVerdict(
        res_pub, res_r, verdict_csv, r_lines, rule_lines, tie_lines, failures
    )
