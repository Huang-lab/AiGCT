"""Removing variants shared between two sets.

Two distinct operations were previously tangled in one function:

`remove_overlapping` enforces a clean split between the positive and negative
sets of a task. A de novo variant seen in both a case cohort and the shared
control set carries no signal, so it is dropped from both sides. This is what
produces the final negative-set sizes quoted in the Methods (1,770 for ASD,
1,798 for CHD, 1,788 for DDD).

`exclude_clinvar_variants` produces the ClinVar-excluded developmental
disorder benchmarks (Figure 4B, D, F). Here one group is the ClinVar variant
set and the aim is a single concatenated table of everything not present in
it, used to quantify how much of a predictor's ranking depends on variants it
may have been trained on.

Both write an append-only log of the sizes before and after, so a run can be
audited after the fact.
"""
import os
from contextlib import redirect_stdout

import pandas as pd

from .columns import VARIANT_KEY
from .config import output_path

_LOG_NAME = "overlap.log"


def _load(paths, base_dir):
    """Read and concatenate variant tables, de-duplicated on the variant key."""
    frames = [
        pd.read_csv(os.path.join(base_dir, p)).drop_duplicates(subset=VARIANT_KEY)
        for p in paths
    ]
    return pd.concat(frames, ignore_index=True)


def _difference(df, other):
    """Rows of `df` whose variant key does not appear in `other`."""
    df = df.copy()
    other = other[VARIANT_KEY].copy()
    # Chromosome is read as int for some inputs and str for others.
    df["#chr"] = df["#chr"].astype(str)
    other["#chr"] = other["#chr"].astype(str)
    merged = df.merge(other, on=VARIANT_KEY, how="left", indicator=True)
    return merged[merged["_merge"] == "left_only"].drop(columns="_merge")


def remove_overlapping(
    group_a,
    group_b,
    base_dir=None,
    log_file=None,
    mode="both",
    suffix="_no_overlap",
):
    """Drop variants shared between two groups of files.

    Each file is written back out as `<name><suffix>.csv` containing only the
    variants absent from the opposite group. `mode` selects which side is
    filtered: "left" (group A only), "right" (group B only) or "both".
    """
    base_dir = base_dir or output_path("processed")
    log_file = log_file or output_path("logs", _LOG_NAME)

    written = []
    with open(log_file, "a", encoding="utf-8") as f:
        with redirect_stdout(f):
            df_a_all = _load(group_a, base_dir)
            df_b_all = _load(group_b, base_dir)
            print(f"Initial: total A={len(df_a_all)}, total B={len(df_b_all)}")

            plan = []
            if mode in ("left", "both"):
                plan.append(("A", group_a, df_b_all))
            if mode in ("right", "both"):
                plan.append(("B", group_b, df_a_all))

            for label, group, opposite in plan:
                print(f"\nProcessing Group {label} files:")
                for name in group:
                    src = os.path.join(base_dir, name)
                    df = pd.read_csv(src).drop_duplicates(subset=VARIANT_KEY)
                    kept = _difference(df, opposite)
                    out = os.path.join(
                        base_dir, os.path.basename(name).replace(".csv", suffix + ".csv")
                    )
                    kept.to_csv(out, index=False)
                    written.append(out)
                    print(f"{name}: kept {len(kept)} of {len(df)} -> {out}")

            print("\nDone. Non-overlapping files written to:", base_dir)
    return written


def exclude_clinvar_variants(
    dataset_files,
    clinvar_files,
    output_name,
    base_dir=None,
    log_file=None,
):
    """Concatenate `dataset_files` with every ClinVar variant removed.

    Used to build the ClinVar-excluded ASD, CHD and DDD benchmarks. Returns
    the path of the single concatenated table that was written.
    """
    base_dir = base_dir or output_path("processed")
    log_file = log_file or output_path("logs", _LOG_NAME)
    out = os.path.join(base_dir, output_name)

    with open(log_file, "a", encoding="utf-8") as f:
        with redirect_stdout(f):
            clinvar = _load(clinvar_files, base_dir)
            print(f"\nExcluding {len(clinvar)} ClinVar variants from {output_name}:")

            kept_frames = []
            for name in dataset_files:
                df = pd.read_csv(os.path.join(base_dir, name)).drop_duplicates(
                    subset=VARIANT_KEY
                )
                kept = _difference(df, clinvar)
                kept_frames.append(kept)
                print(f"{name}: kept {len(kept)} of {len(df)}")

            combined = pd.concat(kept_frames, ignore_index=True)
            combined.to_csv(out, index=False)
            print(f"-> {out} ({len(combined)} variants)")
    return out
