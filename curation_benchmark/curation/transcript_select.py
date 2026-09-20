"""Transcript hierarchy for resolving duplicate variants.

dbNSFP reports one record per transcript, so a variant that falls in several
transcripts appears several times with identical coordinates. This module
collapses each such set to a single representative transcript, following the
hierarchy described in the Methods and Supplementary Figure S6:

  1. the canonical transcript as annotated by dbNSFP (VEP_canonical == 'YES'),
     where there is exactly one;
  2. failing that, the Ensembl canonical transcript;
  3. failing that, the CCDS-annotated transcript with the longest CDS, ties
     broken by overall transcript length;
  4. failing that, the transcript with the longest overall length (UTRs + CDS).

A tiebreak whose column is missing for every candidate is skipped rather than
applied, so it cannot empty the group; see `_keep_max`.

Reference tables (Ensembl BioMart export and CCDS CDS lengths) are loaded
lazily on first use so that importing this module is cheap and does not fail
when the tables are absent.
"""
from functools import lru_cache

import pandas as pd

from .columns import MULTI_TRANSCRIPT_COLUMNS, VARIANT_KEY
from .config import load_config, output_path, repo_path

_LENGTH_COLUMN = "Transcript length (including UTRs and CDS)"


@lru_cache(maxsize=1)
def transcript_reference():
    """Ensembl transcripts joined to their CCDS CDS length.

    Rows present only in the CCDS table (no matching Ensembl transcript) are
    dropped; Ensembl transcripts with no CCDS entry are kept with a null
    'count' so that the transcript-length tiebreak can still use them.
    """
    cfg = load_config()["repo"]
    biomart = pd.read_csv(repo_path(cfg["biomart_export"]), sep="\t")
    ccds = pd.read_csv(repo_path(cfg["ccds_length_table"]))

    # CCDS IDs arrive as "CCDS12345.1|<description>"; reduce to the bare ID.
    ccds["ID"] = ccds["ID"].str.split("|").str[0].str.split(".").str[0]

    merged = pd.merge(
        biomart, ccds, left_on="CCDS ID", right_on="ID", how="outer", indicator=True
    )
    merged = merged.drop(merged[merged["_merge"] == "right_only"].index)
    merged = merged.drop(columns=["ID", "_merge"])
    merged["Ensembl Canonical"] = merged["Ensembl Canonical"].astype(str)
    return merged


def _canonical_index(value):
    """Position of the canonical transcript within a ';'-delimited field."""
    try:
        return int(value.split(";").index("YES"))
    except (AttributeError, ValueError):
        return -1


def _element_at(value, index):
    """Element `index` of a ';'-delimited field, or the field itself if -1."""
    if index == -1:
        return value
    return value.split(";")[index]


def _keep_max(frame, column):
    """Rows where `column` is at its maximum, or `frame` unchanged if all NaN.

    `frame[column].max()` is NaN when every value in the column is missing, and
    `NaN == NaN` is False, so comparing against it drops every row and deletes
    the variant from the benchmark with no message and no duplicates entry. An
    absent CDS length is not evidence about which transcript to keep, so the
    tiebreak falls through to the next one instead of emptying the group.
    """
    if frame[column].isna().all():
        return frame
    return frame[frame[column] == frame[column].max()]


def _apply_hierarchy(group):
    """Pick one row from `group` using the CCDS/length hierarchy."""
    if group[group["_merge"] != "left_only"].empty:
        print("no reference match for variant:\n", group)
        return group

    # Step 1: Ensembl canonical.
    canonical = group[group["Ensembl Canonical"] == "1.0"]
    if canonical.shape[0] == 1:
        return canonical
    if not canonical.empty:
        canonical = _keep_max(canonical, "count")
        canonical = _keep_max(canonical, _LENGTH_COLUMN)
        if canonical.shape[0] > 1:
            print("tie among canonical transcripts:\n", canonical, "\n")
        return canonical

    # Step 2: CCDS-annotated, longest CDS then longest transcript.
    ccds = group[~pd.isna(group["CCDS ID"])]
    if ccds.shape[0] == 1:
        return ccds
    if not ccds.empty:
        ccds = _keep_max(ccds, "count")
        ccds = _keep_max(ccds, _LENGTH_COLUMN)
        if ccds.shape[0] > 1:
            print("tie among CCDS transcripts:\n", ccds, "\n")
        return ccds

    # Step 3: longest overall transcript.
    longest = _keep_max(group, _LENGTH_COLUMN)
    if longest.shape[0] > 1:
        print("tie on transcript length:\n", longest, "\n")
    return longest


def _select_by_reference(value):
    """Explode the per-transcript fields and apply the hierarchy."""
    for col in MULTI_TRANSCRIPT_COLUMNS:
        value[col] = value[col].str.split(";")
    value = value.explode(MULTI_TRANSCRIPT_COLUMNS).reset_index(drop=True)
    value = pd.merge(
        value,
        transcript_reference(),
        left_on="Ensembl_transcriptid",
        right_on="Transcript stable ID",
        how="left",
        indicator=True,
    )
    value = value.drop_duplicates(subset=VARIANT_KEY + ["Ensembl_transcriptid"])
    return _apply_hierarchy(value).reset_index(drop=True)


def _resolve_group(group):
    """Resolve one variant's rows, preferring dbNSFP's own canonical flag.

    A variant with exactly one dbNSFP canonical transcript is settled there.
    With none, or with more than one, the choice is delegated to the
    Ensembl/CCDS/length hierarchy in `_select_by_reference`.
    """
    if group.shape[0] == 1:
        if group["VEP_canonical"].isin(["YES", "."]).any():
            return group
        return _select_by_reference(group)
    canonical = group[group["VEP_canonical"] == "YES"]
    if canonical.shape[0] == 1:
        return canonical
    return _select_by_reference(group)


def choose_canonical(file, dataset_name):
    """Reduce `file` to one row per variant.

    Variants that still have more than one row afterwards are written to the
    configured duplicates directory as `<dataset_name>.csv` for inspection.
    """
    file = file.copy()
    file["num"] = file["VEP_canonical"].apply(_canonical_index)
    for col in MULTI_TRANSCRIPT_COLUMNS:
        file[col] = file.apply(lambda row: _element_at(row[col], row["num"]), axis=1)

    grouped = file.groupby(VARIANT_KEY, as_index=False)
    result = grouped.apply(_resolve_group).reset_index(drop=True)
    result = result.drop_duplicates(subset=VARIANT_KEY + ["Ensembl_transcriptid"])

    unresolved = result[result.duplicated(subset=VARIANT_KEY, keep=False)]
    if not unresolved.empty:
        unresolved.to_csv(output_path("duplicates", dataset_name + ".csv"), index=False)

    result = result.drop(
        columns=[
            "Transcript stable ID",
            "Transcript stable ID version",
            _LENGTH_COLUMN,
            "Ensembl Canonical",
            "CCDS ID",
            "count",
            "_merge",
            "VEP_canonical",
            "num",
        ],
        axis=1,
        errors="ignore",
    )
    print("initial len:", len(file), "\nfinal len:", len(result), "\n")
    return result
