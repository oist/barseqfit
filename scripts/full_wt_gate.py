"""Full-data WT gate: Python (vectorised) vs scalar sample vs published vs R.

Run via ``just full-wt-gate`` after ``just export-papkou``. Writes
``verdict_full_wt.md``, ``verdict_full_wt.csv`` and ``timing_python.json`` to
``$BARSEQFIT_SCRATCH`` (default ``.scratch/papkou2023/``). Inputs are Parquet;
parameters come from ``validation.toml``, are read in ``main`` and
recorded in ``timing_python.json``. The rule itself is ``validation.gate``.

Declared rule: ``m`` is held to the fixture rule on every
non-reference row, vs published and vs R. ``se`` and ``p`` vs R are held to the
fixture rule on rows where the Python and R IRLS iteration counts are equal;
rows with different counts are "threshold ties" (the 1e-8 deviance convergence
test falls on the other side of its threshold), bounded in number
(``tie_fraction``) and size (``tie_se_abs``, ``tie_p_abs``), all from
``validation.toml`` ``[gate]``.
``se`` and ``p`` vs published are held to the fixture rule on rows that are
neither ties nor rows where R itself is outside tolerance vs published; those
R-vs-published outliers are capped at the same ``tie_fraction`` of rows.

The published, fitted and (when present) R id sets must be equal, and the R
table must have been fitted from the current export (``timing_r.json`` carries
the export manifest timestamp and store md5s; a mismatch exits 2).

Exit code is non-zero only for an asserted failure. Status disagreements, the
before-only list, warning counts and the zero before-count list are findings.
"""

from __future__ import annotations

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
from validation.identity import IdentityRule, sequences_identity

PROG = "full_wt_gate"
REFERENCE_ID = "GCCGATCTC"
NOT_REFERENCE = pl.col("sequence_id") != REFERENCE_ID
KEYS = ["sequence_id"]
BEFORE_ONLY_SHOWN = 10


def fit(
    d: dict,
    engine: str,
    chunk_size: int,
    ids: list[str] | None = None,
    **kwargs,
) -> pl.DataFrame:
    sequences = d["sequences"]
    if ids is not None:
        keep = set(ids)
        sequences = sequences.filter(pl.col("sequence_id").is_in(list(keep)))
    counts = d["counts"].join(sequences.select("sequence_id"), on="sequence_id")
    return bsf.fit_wt(
        counts,
        d["samples"],
        sequences,
        reference_id=REFERENCE_ID,
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
    identity_rule = IdentityRule.from_run_config(config)
    chunk_size = config["fit"]["chunk_size"]
    sample_size = config["gate"]["scalar_sample_size"]
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
            "expected_wt.parquet",
        ],
        "fit_wt_full_r.parquet",
        "timing_r.json",
        "Re-run `just full-wt-r` (or `just validate-full-wt`).",
    )

    samples = pl.read_parquet(scratch / "samples.parquet")
    d = {
        "counts": pl.read_parquet(scratch / "counts.parquet"),
        "samples": samples,
        "sequences": pl.read_parquet(scratch / "sequences.parquet"),
        "elapsed": float(samples["time"].max() - samples["time"].min()),
    }
    expected = pl.read_parquet(scratch / "expected_wt.parquet")
    r_path = scratch / "fit_wt_full_r.parquet"
    r_table = pl.read_parquet(r_path) if r_path.exists() else None

    vec, t_vec, fallbacks = fit_and_count(lambda: fit(d, "vectorised", chunk_size))
    print(f"vectorised fit: {t_vec:.1f}s, {vec.height} rows", flush=True)

    # Scalar vs vectorised on a seeded sample.
    rng = np.random.default_rng(0)
    others = [s for s in vec["sequence_id"].to_list() if s != REFERENCE_ID]
    sample_ids = [
        others[i] for i in rng.choice(len(others), sample_size, replace=False)
    ]
    scal, t_scalar, _ = fit_and_count(
        lambda: fit(d, "scalar", chunk_size, [REFERENCE_ID, *sample_ids])
    )
    vec_sample = vec.filter(pl.col("sequence_id").is_in(scal["sequence_id"].implode()))
    sv = scal.join(vec_sample, on="sequence_id", suffix="_v").filter(NOT_REFERENCE)
    sv_stats = scalar_vs_vectorised(sv, rule.scalar_tol)
    scalar_ok = sv_stats["ok"]

    # Three-way comparison.
    fitted_ids = keyset(vec, KEYS)
    pub_ids = keyset(expected, KEYS)
    lines_sets = set_diff_lines("published", pub_ids, "fitted", fitted_ids)
    failures: list[str] = []
    if pub_ids != fitted_ids:
        failures.append(
            f"published ids != fitted ids ({len(pub_ids - fitted_ids)} only published, "
            f"{len(fitted_ids - pub_ids)} only fitted)"
        )
    before_only_path = scratch / "before_only_ids.parquet"
    if before_only_path.exists():
        before_only = pl.read_parquet(before_only_path)["sequence_id"].to_list()
        before_only_lines = [
            f"- before-only genotypes (in the before matrix, absent from after; not fitted): {len(before_only)}",
            f"- first {min(BEFORE_ONLY_SHOWN, len(before_only))}: "
            + ", ".join(sorted(before_only)[:BEFORE_ONLY_SHOWN]),
        ]
        stray = set(before_only) & fitted_ids
        if stray:
            failures.append(
                f"before-only ids present in the fitted table: {len(stray)}"
            )
    else:
        before_only_lines = ["- before_only_ids.parquet absent"]

    # The reference row must be present in every table before it is dropped.
    ref_missing = [
        name
        for name, ids in (
            ("fitted", keyset(vec, KEYS)),
            ("published", keyset(expected, KEYS)),
            *((("R", keyset(r_table, KEYS)),) if r_table is not None else ()),
        )
        if REFERENCE_ID not in ids
    ]
    for name in ref_missing:
        failures.append(f"reference row {REFERENCE_ID} missing from {name} table")

    # The reference row is not compared: drop it before the chain, then restore
    # it (with nulls) in the verdict CSV so the file keeps one row per genotype.
    chain = compare_chain(
        vec.filter(NOT_REFERENCE),
        expected.filter(NOT_REFERENCE),
        None if r_table is None else r_table.filter(NOT_REFERENCE),
        KEYS,
        rule,
    )
    failures += chain.failures
    verdict_csv = vec.select(*KEYS, "status", "m", "se", "p").join(
        chain.verdict_csv.drop(["status", "m", "se", "p"]), on=KEYS, how="left"
    )
    if r_table is not None and "status" in r_table.columns:
        # R carries a status for the reference row too; keep it in the CSV.
        verdict_csv = (
            verdict_csv.join(
                r_table.select(*KEYS, pl.col("status").alias("_status_r_all")),
                on=KEYS,
                how="left",
            )
            .with_columns(pl.col("_status_r_all").alias("status_r"))
            .drop("_status_r_all")
        )

    status_ls, status_disagree = status_lines(vec, r_table, KEYS)
    warn = warning_counts(vec)
    warn_lines = [
        f"- {r['len']} x {r['warning']}" for r in warn.iter_rows(named=True)
    ] or ["- none"]

    # Pseudocount check.
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
    zero_ids = [s for s in zero_ids if s != REFERENCE_ID]
    pseudo_header = f"Genotypes with any zero before-count: {len(zero_ids)}"
    pseudo_lines: list[str] = [pseudo_header]
    pseudo_json: dict = {"n_zero_before": len(zero_ids)}
    zero_lines = subset_residual_lines(
        vec,
        expected,
        KEYS,
        pl.DataFrame({"sequence_id": zero_ids}, schema={"sequence_id": pl.String}),
        rule.tolerance,
    )

    # fit_sequences identity: m == r_i - r_ref and se == sqrt(se_pooled_i^2 + se_pooled_ref^2).
    seq_fit, t_seq, _ = fit_and_count(
        lambda: bsf.fit_sequences(
            d["counts"],
            d["samples"],
            d["sequences"],
            elapsed_time=d["elapsed"],
            time_unit=TIME_UNIT,
            pseudocount="all_zeros",
            chunk_size=chunk_size,
        )
    )
    print(f"fit_sequences: {t_seq:.1f}s, {seq_fit.height} rows", flush=True)
    identity = sequences_identity(seq_fit, vec, REFERENCE_ID, identity_rule)
    failures += [f"fit_sequences identity: {f}" for f in identity.failures]
    identity_lines = [
        f"- fit_sequences fit: {t_seq:.1f} s, {seq_fit.height} rows",
        *identity.lines(identity_rule),
    ]

    header_lines = [
        f"- Declared rule: m asserted on every non-reference row vs published and vs R. se, p vs R asserted on rows with equal Python and R iteration counts; rows with different counts are threshold ties, asserted at <= {rule.tie_fraction:.2%} of rows with se and p residuals <= {rule.tie_se_abs:g} (se) and {rule.tie_p_abs:g} (p). se, p vs published asserted on rows that are neither ties nor R-vs-published outliers (outliers also capped at {rule.tie_fraction:.2%} of rows); the rest are reported. Published, fitted and R id sets asserted equal.",
        f"- Genotypes fitted: {vec.height} (published rows: {expected.height})",
        f"- Vectorised fit: {t_vec:.1f} s; scalar fit on {sample_size} sampled genotypes: {t_scalar:.1f} s",
        f"- Scalar-oracle fallbacks inside the vectorised fit: {fallbacks}",
        f"- R table present: {r_table is not None}",
    ]
    md, passed = render_verdict_md(
        title="# Full-data WT gate verdict",
        header_lines=header_lines,
        sv_stats=sv_stats,
        lines_sets=lines_sets,
        before_only_lines=before_only_lines,
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
            ),
            ("fit_sequences identity with fit_wt", identity_lines),
        ],
    )
    write_verdict(
        scratch,
        "verdict_full_wt",
        md,
        verdict_csv,
        {
            "vectorised_s": t_vec,
            "scalar_sample_s": t_scalar,
            "genotypes": vec.height,
            "scalar_sample_size": sample_size,
            "fallbacks": fallbacks,
            "pseudocount": pseudo_json,
            "fit_sequences_s": t_seq,
            "identity": {
                "checked": identity.checked,
                "skipped": identity.skipped,
                "max_m": identity.max_m,
                "max_se": identity.max_se,
                "outside_m": len(identity.outside_m),
                "outside_se": len(identity.outside_se),
            },
            "config": config,
        },
        "timing_python.json",
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
