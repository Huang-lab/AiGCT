"""Build the ClinVar benchmark from the annotations embedded in dbNSFP.

The ClinVar task does not start from a coordinate list. Instead the whole
dbNSFP v5.0a release is scanned and every record carrying a ClinVar
classification is retained, using dbNSFP's own `clinvar_clnsig` and
`clinvar_review` fields. dbNSFP therefore fixes the ClinVar release used; see
the dbNSFP readme for the exact date.

Retained variants are those classified pathogenic / likely pathogenic or
benign / likely benign with a review status of one to four stars. Records
whose alternate amino acid is a stop ('X') or missing ('.') are dropped so
that the benchmark covers missense variants only.
"""
import pandas as pd
import tqdm

from .columns import CLINVAR_COLUMNS, COLUMN_INFO, VARIANT_KEY, VEP_LIST, VEP_RENAME
from .config import dbnsfp_chromosome_file, load_config
from .transcript_select import choose_canonical

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["M", "X", "Y"]

_SCAN_COLUMNS = COLUMN_INFO + CLINVAR_COLUMNS + VEP_LIST
_CHUNK_SIZE = 10_000


def scan_dbnsfp():
    """Return every dbNSFP record carrying a retained ClinVar classification."""
    cfg = load_config()["clinvar"]
    clnsig_labels = cfg["clnsig_labels"]
    review_stars = cfg["review_stars"]

    collected = []
    for chrom in CHROMOSOMES:
        print("processing chr" + chrom)
        chunks = pd.read_csv(
            dbnsfp_chromosome_file(chrom),
            sep="\t",
            usecols=_SCAN_COLUMNS,
            dtype=str,
            chunksize=_CHUNK_SIZE,
        )
        for chunk in tqdm.tqdm(chunks):
            chunk = chunk[chunk["clinvar_clnsig"].isin(clnsig_labels)]
            chunk = chunk[chunk["clinvar_review"].isin(review_stars)]
            if chunk.empty:
                continue
            chunk = chunk.copy()
            chunk["BINARY_LABEL"] = chunk["clinvar_clnsig"].map(clnsig_labels)
            chunk["clinvar_review"] = chunk["clinvar_review"].map(review_stars)
            collected.append(chunk)

    if not collected:
        return pd.DataFrame(columns=_SCAN_COLUMNS + ["BINARY_LABEL"])
    return pd.concat(collected, ignore_index=True).rename(columns=VEP_RENAME)


def drop_non_missense(df):
    """Remove records whose alternate amino acid is a stop or is missing."""
    return df[~df["aaalt"].isin(["X", "."])]


def build(out_path, dataset_name="clinvar"):
    """Scan, transcript-resolve and filter the full ClinVar benchmark.

    Note the ordering: `aaalt` is a ';'-delimited per-transcript field in the
    raw release, so the stop/missing filter is only meaningful once a single
    transcript has been chosen.
    """
    raw = scan_dbnsfp()
    resolved = choose_canonical(raw, dataset_name)
    resolved = resolved.drop_duplicates(subset=VARIANT_KEY)
    resolved = drop_non_missense(resolved)
    resolved.to_csv(out_path, index=False)
    print("ClinVar benchmark variants:", len(resolved))
    return resolved
