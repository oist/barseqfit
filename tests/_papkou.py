"""Helpers shared by the Papkou2023 fixture tests."""

import polars as pl

from validation.compare import assert_close  # noqa: F401

REFERENCE_ID = "GCCGATCTC"
TIME_UNIT = "generations"


def load(data_dir):
    def read(name, **kwargs):
        return pl.read_csv(data_dir / name, **kwargs)

    samples = read("samples.csv")
    elapsed = samples["time"].max() - samples["time"].min()
    return {
        "counts": read("counts.csv"),
        "samples": samples,
        "sequences": read("sequences.csv"),
        "pairs": read("pairs.csv"),
        "elapsed": elapsed,
    }
