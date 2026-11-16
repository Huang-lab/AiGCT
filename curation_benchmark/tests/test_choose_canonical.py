"""choose_canonical: the whole transcript-resolution step, end to end.

`test_transcript_select.py` pins the hierarchy (`_apply_hierarchy`) on a frame
that is already merged. Nothing there calls `choose_canonical` itself, which is
the function every `annotate` and `clinvar` run ends with, hours after the
dbNSFP scan has finished. On pandas 3 it raised before returning anything:

* `DataFrameGroupBy.apply` no longer hands the grouping columns to the function,
  so the `drop_duplicates(subset=VARIANT_KEY + ...)` that follows raised
  `KeyError: ['#chr', 'pos(1-based)', 'ref', 'alt']`;
* `DataFrame.drop(columns=..., axis=1)` raises "Cannot specify both 'axis' and
  'index'/'columns'".

These tests run it on a small raw table with a stand-in for the BioMart/CCDS
reference, so they need neither dbNSFP nor the bundled reference tables, and
they hold on every pandas version the pipeline claims to support.

    pytest curation_benchmark/tests
"""
import io
import os
import sys
from contextlib import redirect_stdout

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation import transcript_select  # noqa: E402
from curation.columns import (  # noqa: E402
    EXTRACT_COLUMNS,
    MULTI_TRANSCRIPT_COLUMNS,
    VARIANT_KEY,
)

LENGTH = transcript_select._LENGTH_COLUMN


def _raw(chrom, pos, ref, alt, transcripts, flags):
    """One dbNSFP row. Per-transcript fields are ';'-joined, as dbNSFP writes them;
    `genename` carries a distinct marker per transcript so the element that was
    picked can be read back."""
    row = {c: "." for c in EXTRACT_COLUMNS}
    for col in MULTI_TRANSCRIPT_COLUMNS:
        row[col] = ";".join(["."] * len(transcripts))
    row.update({
        "#chr": chrom, "pos(1-based)": pos, "ref": ref, "alt": alt,
        "Ensembl_transcriptid": ";".join(transcripts),
        "genename": ";".join(f"G_{t}" for t in transcripts),
        "VEP_canonical": ";".join(flags),
    })
    return row


def _reference(*rows):
    """(transcript, ensembl_canonical, cds_length, transcript_length, ccds_id)."""
    return pd.DataFrame([{
        "Transcript stable ID": t,
        "Transcript stable ID version": t + ".1",
        LENGTH: length,
        "Ensembl Canonical": canonical,
        "CCDS ID": ccds,
        "count": cds,
    } for t, canonical, cds, length, ccds in rows])


def _run(raw_rows, reference, tmp_path, monkeypatch):
    monkeypatch.setattr(transcript_select, "transcript_reference", lambda: reference)
    monkeypatch.setattr(
        transcript_select, "output_path",
        lambda kind, name=None: str(tmp_path / (name or kind)),
    )
    raw = pd.DataFrame(raw_rows, columns=EXTRACT_COLUMNS)
    with redirect_stdout(io.StringIO()):
        return transcript_select.choose_canonical(raw, "ds")


# Five variants, deliberately not in key order, one per way of being resolved.
RAW = [
    # dbNSFP's own canonical flag settles it (single transcript)
    _raw("2", "100", "A", "G", ["ENST_A"], ["YES"]),
    # canonical flag on the second of two transcripts: that element is taken
    _raw("1", "200", "C", "T", ["ENST_B1", "ENST_B2"], [".", "YES"]),
    # no dbNSFP flag: the Ensembl canonical transcript wins
    _raw("1", "300", "G", "A", ["ENST_C1", "ENST_C2"], [".", "."]),
    # no dbNSFP flag and neither transcript is in the reference: kept whole
    _raw("1", "400", "T", "C", ["ENST_X1", "ENST_X2"], [".", "."]),
    # two Ensembl-canonical transcripts, neither with a CDS length: the longer wins
    _raw("1", "500", "A", "T", ["ENST_D1", "ENST_D2"], [".", "."]),
]
REFERENCE = _reference(
    ("ENST_C1", "nan", np.nan, 5000, np.nan),
    ("ENST_C2", "1.0", np.nan, 1000, np.nan),
    ("ENST_D1", "1.0", np.nan, 3000, np.nan),
    ("ENST_D2", "1.0", np.nan, 2000, np.nan),
)


def test_every_variant_is_resolved_and_the_helper_columns_are_dropped(tmp_path, monkeypatch):
    result = _run(RAW, REFERENCE, tmp_path, monkeypatch)

    chosen = {tuple(k): t for *k, t in zip(
        result["#chr"], result["pos(1-based)"], result["ref"], result["alt"],
        result["Ensembl_transcriptid"])}
    assert chosen[("2", "100", "A", "G")] == "ENST_A"
    assert chosen[("1", "200", "C", "T")] == "ENST_B2"
    assert chosen[("1", "300", "G", "A")] == "ENST_C2"
    assert chosen[("1", "500", "A", "T")] == "ENST_D1"

    helper = {"Transcript stable ID", "Transcript stable ID version", LENGTH,
              "Ensembl Canonical", "CCDS ID", "count", "_merge", "VEP_canonical", "num"}
    assert helper.isdisjoint(result.columns)
    assert set(VARIANT_KEY) <= set(result.columns)


def test_the_element_at_the_canonical_position_is_taken_for_every_field(tmp_path, monkeypatch):
    result = _run(RAW, REFERENCE, tmp_path, monkeypatch)
    row = result[(result["#chr"] == "1") & (result["pos(1-based)"] == "200")].iloc[0]

    assert row["Ensembl_transcriptid"] == "ENST_B2"
    assert row["genename"] == "G_ENST_B2"


def test_rows_come_back_in_sorted_variant_order(tmp_path, monkeypatch):
    """balance.py samples with a fixed seed, which depends on row order."""
    result = _run(RAW, REFERENCE, tmp_path, monkeypatch)
    keys = list(zip(result["#chr"], result["pos(1-based)"]))

    assert keys == sorted(keys)


def test_an_unresolvable_variant_is_kept_and_reported_in_duplicates(tmp_path, monkeypatch):
    result = _run(RAW, REFERENCE, tmp_path, monkeypatch)

    unresolved = result[result["pos(1-based)"] == "400"]
    assert sorted(unresolved["Ensembl_transcriptid"]) == ["ENST_X1", "ENST_X2"]
    reported = pd.read_csv(tmp_path / "ds.csv", dtype=str)
    assert sorted(reported["Ensembl_transcriptid"]) == ["ENST_X1", "ENST_X2"]
    assert len(result) == 6          # five variants, one of them left as two rows


def test_a_run_with_nothing_unresolved_writes_no_duplicates_file(tmp_path, monkeypatch):
    result = _run(RAW[:3], REFERENCE, tmp_path, monkeypatch)

    assert len(result) == 3
    assert not (tmp_path / "ds.csv").exists()
