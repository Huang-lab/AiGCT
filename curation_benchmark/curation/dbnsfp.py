"""Annotate variant coordinate lists with dbNSFP v5.0a records.

Each cohort contributes a whitespace-delimited list of variants (chrom, pos,
ref, alt) under data/annotation. This module joins those lists against the
per-chromosome dbNSFP release files, keeps the annotation and VEP score
columns, and hands the result to the transcript hierarchy.

Coordinates are matched against dbNSFP's primary (hg38) columns or its hg19
columns depending on the assembly the source study reported. dbNSFP carries
hg19 positions but only hg38 alleles, so an hg19 list is joined on hg19
chromosome and position and on the hg38 ref/alt; a locus whose reference base
changed between builds therefore fails to match. hg18 lists are not supported.

Every variant that finds no dbNSFP record is written, in the input format, to
the configured `unmatched` directory as `<dataset>.txt`, and the per-dataset
match count is appended to `match_rates.csv` there, so that the loss at this
stage is on record rather than in a terminal scrollback.
"""
import os

import pandas as pd
import tqdm

from .columns import EXTRACT_COLUMNS, VEP_RENAME
from .config import dbnsfp_chromosome_file, output_path
from .transcript_select import choose_canonical

ANNOTATION_COLUMNS = ["Chr", "Pos", "Ref", "Alt"]

# dbNSFP columns to join against, per source assembly.
_JOIN_COLUMNS = {
    "hg38": ["#chr", "pos(1-based)", "ref", "alt"],
    "hg19": ["hg19_chr", "hg19_pos(1-based)", "ref", "alt"],
}

_CHUNK_SIZE = 1_000_000
_MATCH_LOG = "match_rates.csv"


def read_annotation(path):
    """Read a whitespace-delimited, headerless variant list.

    Split on runs of whitespace rather than on a single space: two of the
    committed coordinate lists carry a trailing space on every line, which a
    single-space split turns into a fifth, empty field. pandas then shifts the
    columns left, the alternate allele reads as missing, and nothing matches.
    """
    frame = pd.read_csv(
        path, sep=r"\s+", header=None, names=ANNOTATION_COLUMNS, dtype=str,
        engine="python",
    )
    if frame[ANNOTATION_COLUMNS].isna().any().any():
        raise ValueError(f"{path}: malformed line(s); expected 'chr pos ref alt'")
    return frame


def extract(annotation, assembly):
    """Pull dbNSFP records for every variant in `annotation`.

    Reads one chromosome file at a time in chunks, since the full release is
    far too large to hold in memory. Returns the matched records and the
    annotation rows that matched nothing.
    """
    if assembly not in _JOIN_COLUMNS:
        raise ValueError(
            f"assembly must be one of {sorted(_JOIN_COLUMNS)}, got {assembly!r}"
        )
    right_on = _JOIN_COLUMNS[assembly]

    print("input variant size = " + str(len(annotation)))
    pieces = []

    for chrom, group in annotation.groupby("Chr"):
        print("processing chr" + str(chrom))
        chunks = pd.read_csv(
            dbnsfp_chromosome_file(chrom),
            sep="\t",
            usecols=EXTRACT_COLUMNS,
            dtype=str,
            chunksize=_CHUNK_SIZE,
        )
        for data in tqdm.tqdm(chunks):
            both = pd.merge(
                group, data, left_on=ANNOTATION_COLUMNS, right_on=right_on, how="inner"
            )
            if len(both):
                pieces.append(both)

    if pieces:
        matched = pd.concat(pieces, ignore_index=True)
    else:
        matched = pd.DataFrame(columns=ANNOTATION_COLUMNS + EXTRACT_COLUMNS)

    hit = matched[ANNOTATION_COLUMNS].drop_duplicates()
    unmatched = annotation.merge(hit, on=ANNOTATION_COLUMNS, how="left", indicator=True)
    unmatched = unmatched.loc[unmatched["_merge"] == "left_only", ANNOTATION_COLUMNS]

    extracted = matched.drop(columns=ANNOTATION_COLUMNS).rename(columns=VEP_RENAME)
    return extracted, unmatched


def _record_unmatched(dataset_name, assembly, annotation_file, annotation, unmatched):
    """Write the unmatched coordinates and append the match count to the log."""
    unmatched_path = output_path("unmatched", dataset_name + ".txt")
    if os.path.abspath(unmatched_path) == os.path.abspath(annotation_file):
        raise ValueError(
            f"unmatched output {unmatched_path} would overwrite the input "
            f"coordinate list; refusing to run"
        )
    unmatched.to_csv(unmatched_path, sep=" ", header=False, index=False)

    n_in, n_unmatched = len(annotation), len(unmatched)
    row = pd.DataFrame([{
        "dataset": dataset_name,
        "assembly": assembly,
        "n_input": n_in,
        "n_matched": n_in - n_unmatched,
        "n_unmatched": n_unmatched,
        "match_rate": round((n_in - n_unmatched) / n_in, 4) if n_in else float("nan"),
    }])
    log_path = output_path("unmatched", _MATCH_LOG)
    if os.path.exists(log_path):
        log = pd.read_csv(log_path)
        log = pd.concat([log[log["dataset"] != dataset_name], row], ignore_index=True)
    else:
        log = row
    log.to_csv(log_path, index=False)
    print(f"matched {n_in - n_unmatched} of {n_in} in dbNSFP "
          f"({row['match_rate'].iloc[0]:.1%}); unmatched -> {unmatched_path}")


def annotate_dataset(annotation_file, out_path, assembly, dataset_name):
    """Annotate one cohort file and write the transcript-resolved table.

    Returns the resulting DataFrame, or None if no variant was found in
    dbNSFP.
    """
    if os.path.abspath(out_path) == os.path.abspath(annotation_file):
        raise ValueError(
            f"output {out_path} would overwrite the input coordinate list; "
            f"refusing to run"
        )

    annotation = read_annotation(annotation_file)
    raw, unmatched = extract(annotation, assembly)
    _record_unmatched(dataset_name, assembly, annotation_file, annotation, unmatched)
    if raw.empty:
        print("not in dbNSFP data set")
        return None

    resolved = choose_canonical(raw, dataset_name)
    print("output variant size = " + str(len(resolved)))
    resolved.to_csv(out_path, index=False)
    return resolved
