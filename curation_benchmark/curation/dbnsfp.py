"""Annotate variant coordinate lists with dbNSFP v5.0a records.

Each cohort contributes a whitespace-delimited list of variants (chrom, pos,
ref, alt) under data/annotation. This module joins those lists against the
per-chromosome dbNSFP release files, keeps the annotation and VEP score
columns, and hands the result to the transcript hierarchy.

Coordinates are matched against dbNSFP's primary (hg38) columns or its hg19
columns depending on the assembly the source study reported.
"""
import pandas as pd
import tqdm

from .columns import EXTRACT_COLUMNS, VEP_RENAME
from .config import dbnsfp_chromosome_file
from .transcript_select import choose_canonical

ANNOTATION_COLUMNS = ["Chr", "Pos", "Ref", "Alt"]

# dbNSFP columns to join against, per source assembly.
_JOIN_COLUMNS = {
    "hg38": ["#chr", "pos(1-based)", "ref", "alt"],
    "hg19": ["hg19_chr", "hg19_pos(1-based)", "ref", "alt"],
}

_CHUNK_SIZE = 1_000_000


def read_annotation(path):
    """Read a whitespace-delimited, headerless variant list."""
    return pd.read_csv(
        path, sep=" ", header=None, names=ANNOTATION_COLUMNS, dtype=str
    )


def extract(annotation, assembly):
    """Pull dbNSFP records for every variant in `annotation`.

    Reads one chromosome file at a time in chunks, since the full release is
    far too large to hold in memory.
    """
    if assembly not in _JOIN_COLUMNS:
        raise ValueError(
            f"assembly must be one of {sorted(_JOIN_COLUMNS)}, got {assembly!r}"
        )
    right_on = _JOIN_COLUMNS[assembly]

    print("input variant size = " + str(len(annotation)))
    extracted = pd.DataFrame(columns=EXTRACT_COLUMNS)

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
            merged = pd.merge(
                group,
                data,
                left_on=ANNOTATION_COLUMNS,
                right_on=right_on,
                how="left",
                indicator=True,
            )
            if "both" in merged["_merge"].values:
                both = merged[merged["_merge"] == "both"]
                both = both.drop(columns=ANNOTATION_COLUMNS + ["_merge"])
                extracted = pd.concat([extracted, both], ignore_index=True)

    return extracted.rename(columns=VEP_RENAME)


def annotate_dataset(annotation_file, out_path, assembly, dataset_name):
    """Annotate one cohort file and write the transcript-resolved table.

    Returns the resulting DataFrame, or None if no variant was found in
    dbNSFP.
    """
    annotation = read_annotation(annotation_file)
    raw = extract(annotation, assembly)
    if raw.empty:
        print("not in dbNSFP data set")
        return None

    resolved = choose_canonical(raw, dataset_name)
    print("output variant size = " + str(len(resolved)))
    resolved.to_csv(out_path, index=False)
    return resolved
