"""Full-data pair gate: Python (vectorised) vs scalar sample vs published vs R.

Run via ``just full-pair-gate`` after ``just export-papkou`` and
``just full-pair-r``. Writes ``verdict_full_pair.md``, ``verdict_full_pair.csv``
and ``timing_pair_python.json`` to ``$BARSEQFIT_SCRATCH`` (default
``.scratch/papkou2023/``). Inputs are Parquet; parameters come from
``validation.toml``, are read in ``main`` and recorded in
``timing_pair_python.json``.

The tolerance, tie and R-outlier rules and their caps are those of the WT gate
implemented once in ``validation.gate``. Differences:
rows are ``(from, to)`` pairs, there is no reference row, and there is no
before-only list. The published orientation is copied as data, never derived.

The fitted, published and (when present) R ``(from, to)`` row sets must be
equal (order-insensitive); the R table must have been fitted from the current
export (``timing_pair_r.json`` carries the export manifest timestamp and store
md5s; a mismatch exits 2).

Exit code is non-zero only for an asserted failure. Status disagreements,
warning counts and the zero before-count list are findings.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import polars as pl
from _gate_common import (
    REPO,
    TIME_UNIT,
    check_export,
    fit_and_count,
    load_gate_config,
    render_verdict_md,
    write_verdict,
)

import barseqfit as bsf
from validation.env import resolve_machine_env
from validation.gate import (
    GateRule,
    compare_chain,
    keyset,
    scalar_vs_vectorised,
    set_diff_lines,
    status_lines,
    subset_residual_lines,
    warning_counts,
)

PROG = "full_pair_gate"
KEYS = ["from", "to"]


def fit(
    d: dict, engine: str, pairs: pl.DataFrame, chunk_size: int, **kwargs
) -> pl.DataFrame:
    return bsf.fit_pair(
        d["counts"],
        d["samples"],
        d["sequences"],
        pairs=pairs,
        elapsed_time=d["elapsed"],
        time_unit=TIME_UNIT,
        pseudocount="all_zeros",
        dispersion="fixed",
        engine=engine,
        chunk_size=chunk_size,
        **kwargs,
    )


def main() -> int:
    config = load_gate_config()
    rule = GateRule.from_run_config(config)
    chunk_size = config["fit"]["chunk_size"]
    sample_size = config["gate"]["scalar_sample_size"]
    published_rows = config["gate"]["published_pair_rows"]
    menv = resolve_machine_env(REPO)
    scratch = menv.scratch
    check_export(
        PROG,
        scratch,
        menv.store,
        [
            "counts.parquet",
            "samples.parquet",
            "sequences.parquet",
            "pairs.parquet",
            "expected_pair.parquet",
        ],
        "fit_pair_full_r.parquet",
        "timing_pair_r.json",
        "Re-run `just full-pair-r` (or `just validate-full-pair`).",
    )

    manifest = json.loads((scratch / "export_manifest.json").read_text())
    dropped = manifest.get("row_counts", {}).get("edge_rows_dropped")
    failures: list[str] = []
    # The manifest records the published table the export read; a subset smoke
    # run trims the parquet files but leaves the manifest, so it is checked here.
    published = manifest.get("row_counts", {}).get("expected_pair")
    if published != published_rows:
        failures.append(
            f"export_manifest.json row_counts.expected_pair is {published!r}, "
            f"expected {published_rows} (validation.toml published_pair_rows)"
        )
    if dropped != 0:
        failures.append(
            f"export_manifest.json row_counts.edge_rows_dropped is {dropped!r}, expected 0"
        )

    samples = pl.read_parquet(scratch / "samples.parquet")
    d = {
        "counts": pl.read_parquet(scratch / "counts.parquet"),
        "samples": samples,
        "sequences": pl.read_parquet(scratch / "sequences.parquet"),
        "elapsed": float(samples["time"].max() - samples["time"].min()),
    }
    pairs = pl.read_parquet(scratch / "pairs.parquet").select(KEYS)
    expected = pl.read_parquet(scratch / "expected_pair.parquet")
    r_path = scratch / "fit_pair_full_r.parquet"
    r_table = pl.read_parquet(r_path) if r_path.exists() else None

    for name, table in (("pairs.parquet", pairs), ("expected_pair.parquet", expected)):
        if table.height != published_rows:
            failures.append(
                f"{name} has {table.height} rows, expected {published_rows} "
                "(validation.toml published_pair_rows)"
            )

    vec, t_vec, fallbacks = fit_and_count(
        lambda: fit(d, "vectorised", pairs, chunk_size)
    )
    print(f"vectorised fit: {t_vec:.1f}s, {vec.height} pairs", flush=True)

    # Scalar vs vectorised on a seeded sample of pairs.
    rng = np.random.default_rng(0)
    n_sample = min(sample_size, vec.height)
    idx = np.sort(rng.choice(vec.height, n_sample, replace=False))
    sample_pairs = pairs[idx.tolist()]
    scal, t_scalar, _ = fit_and_count(
        lambda: fit(d, "scalar", sample_pairs, chunk_size)
    )
    sv = scal.join(vec, on=KEYS, suffix="_v")
    sv_stats = scalar_vs_vectorised(sv, rule.scalar_tol)
    scalar_ok = sv_stats["ok"] and sv.height == n_sample

    # Exact (from, to) row-set match across Python, published and R.
    fitted_ids, pub_ids = keyset(vec, KEYS), keyset(expected, KEYS)
    lines_sets = set_diff_lines("published", pub_ids, "fitted", fitted_ids)
    if pub_ids != fitted_ids:
        failures.append(
            f"published pairs != fitted pairs ({len(pub_ids - fitted_ids)} only published, "
            f"{len(fitted_ids - pub_ids)} only fitted)"
        )
    if vec.height != len(fitted_ids) or expected.height != len(pub_ids):
        failures.append("duplicate (from, to) rows in the fitted or published table")

    chain = compare_chain(vec, expected, r_table, KEYS, rule)
    failures += chain.failures
    if r_table is not None and r_table.height != len(keyset(r_table, KEYS)):
        failures.append("duplicate (from, to) rows in the R table")

    status_ls, status_disagree = status_lines(vec, r_table, KEYS)
    warn = warning_counts(vec)
    warn_lines = [
        f"- {r['len']} x {r['warning']}" for r in warn.iter_rows(named=True)
    ] or ["- none"]

    # Pseudocount check on pairs touching a genotype with a zero before-count.
    before_t = float(samples["time"].min())
    before_samples = samples.filter(pl.col("time") == before_t)["sample_id"]
    zero_ids = (
        d["counts"]
        .filter(
            pl.col("sample_id").is_in(before_samples.implode()) & (pl.col("count") == 0)
        )["sequence_id"]
        .unique()
        .to_list()
    )
    zero_pairs = pairs.filter(
        pl.col("from").is_in(zero_ids) | pl.col("to").is_in(zero_ids)
    )
    pseudo_header = f"Genotypes with any zero before-count: {len(zero_ids)}; pairs touching one: {zero_pairs.height}"
    pseudo_lines: list[str] = [pseudo_header]
    pseudo_json: dict = {"n_zero_before": len(zero_ids), "n_pairs": zero_pairs.height}
    zero_lines = subset_residual_lines(
        vec, expected, KEYS, zero_pairs.select("from", "to"), rule.tolerance
    )

    header_lines = [
        f"- Declared rule: m asserted on every pair vs published and vs R. se, p vs R asserted on rows with equal Python and R iteration counts; rows with different counts are threshold ties, asserted at <= {rule.tie_fraction:.2%} of rows with se and p residuals <= {rule.tie_se_abs:g} (se) and {rule.tie_p_abs:g} (p). se, p vs published asserted on rows that are neither ties nor R-vs-published outliers (outliers also capped at {rule.tie_fraction:.2%} of rows); the rest are reported. Fitted, published and R (from, to) row sets asserted equal.",
        f"- Pairs fitted: {vec.height} (published rows: {expected.height})",
        f"- Vectorised fit: {t_vec:.1f} s; scalar fit on {n_sample} sampled pairs: {t_scalar:.1f} s",
        f"- Scalar-oracle fallbacks inside the vectorised fit: {fallbacks}",
        f"- R table present: {r_table is not None}",
    ]
    md, passed = render_verdict_md(
        title="# Full-data pair gate verdict",
        header_lines=header_lines,
        sv_stats=sv_stats,
        lines_sets=lines_sets,
        before_only_lines=None,
        chain=chain,
        status_ls=status_ls,
        status_disagree=status_disagree,
        warn_lines=warn_lines,
        pseudo_lines=pseudo_lines,
        failures=failures,
        scalar_ok=scalar_ok,
        rule=rule,
        extra_sections=[
            (
                "Zero before-count subset (default pseudocount rule)",
                zero_lines,
            )
        ],
    )
    write_verdict(
        scratch,
        "verdict_full_pair",
        md,
        chain.verdict_csv,
        {
            "vectorised_s": t_vec,
            "scalar_sample_s": t_scalar,
            "pairs": vec.height,
            "scalar_sample_size": n_sample,
            "fallbacks": fallbacks,
            "pseudocount": pseudo_json,
            "config": config,
        },
        "timing_pair_python.json",
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
