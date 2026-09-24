"""Transcript hierarchy: a tiebreak must never delete the variant.

`_apply_hierarchy` narrows a group of candidate transcripts by taking the rows
at the maximum of a column. When that column is missing for every candidate,
`max()` is NaN and `NaN == NaN` is False, so a naive comparison drops every row
and the variant leaves the pipeline with no message, no exception and no entry
in the duplicates file. These tests pin the group as non-empty on every path.

They need neither dbNSFP nor the BioMart/CCDS reference tables: the hierarchy
operates on an already-merged frame, so the frame is built directly.

    pytest curation_benchmark/tests
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation.transcript_select import _LENGTH_COLUMN, _apply_hierarchy  # noqa: E402


def _row(transcript, canonical, cds_length, transcript_length, ccds_id):
    return {
        "_merge": "both",
        "Ensembl_transcriptid": transcript,
        "Ensembl Canonical": canonical,
        "count": cds_length,
        _LENGTH_COLUMN: transcript_length,
        "CCDS ID": ccds_id,
    }


def _group(*rows):
    return pd.DataFrame(list(rows))


# ── Ensembl canonical tier ────────────────────────────────────────────────────

def test_two_canonical_transcripts_without_a_cds_length_keep_one():
    """Both Ensembl-canonical, neither in the CCDS length table.

    The CCDS tiebreak has nothing to say here, so the transcript-length
    tiebreak has to settle it. Before the fix this returned an empty frame.
    """
    group = _group(
        _row("ENST1", "1.0", np.nan, 3000, np.nan),
        _row("ENST2", "1.0", np.nan, 2000, np.nan),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1, "variant deleted by an all-NaN CDS-length tiebreak"
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST1"


def test_two_canonical_transcripts_with_no_usable_tiebreak_at_all_are_kept():
    """Neither column can decide. Keep both and let the caller report them."""
    group = _group(
        _row("ENST1", "1.0", np.nan, np.nan, np.nan),
        _row("ENST2", "1.0", np.nan, np.nan, np.nan),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 2, "an undecidable group must survive to duplicates/"


def test_a_present_cds_length_still_wins_over_a_missing_one():
    """The fix must not weaken the tiebreak when the column is usable."""
    group = _group(
        _row("ENST1", "1.0", 900, 2000, "CCDS1"),
        _row("ENST2", "1.0", np.nan, 3000, np.nan),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST1"


def test_the_longest_cds_wins_among_canonical_transcripts():
    group = _group(
        _row("ENST1", "1.0", 900, 2000, "CCDS1"),
        _row("ENST2", "1.0", 1200, 3000, "CCDS2"),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST2"


# ── CCDS tier ─────────────────────────────────────────────────────────────────

def test_two_ccds_transcripts_missing_from_the_length_table_keep_one():
    """CCDS IDs that BioMart carries but the NCBI CCDS release has retired."""
    group = _group(
        _row("ENST1", "nan", np.nan, 3000, "CCDS1"),
        _row("ENST2", "nan", np.nan, 2000, "CCDS2"),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1, "variant deleted by an all-NaN CDS-length tiebreak"
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST1"


def test_the_longest_cds_wins_among_ccds_transcripts():
    group = _group(
        _row("ENST1", "nan", 900, 3000, "CCDS1"),
        _row("ENST2", "nan", 1200, 2000, "CCDS2"),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST2"


# ── Length tier and the no-reference-match path ───────────────────────────────

def test_the_longest_transcript_wins_when_no_row_is_canonical_or_ccds():
    group = _group(
        _row("ENST1", "nan", np.nan, 2000, np.nan),
        _row("ENST2", "nan", np.nan, 3000, np.nan),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 1
    assert result["Ensembl_transcriptid"].iloc[0] == "ENST2"


def test_a_group_with_no_transcript_length_at_all_is_kept():
    group = _group(
        _row("ENST1", "nan", np.nan, np.nan, np.nan),
        _row("ENST2", "nan", np.nan, np.nan, np.nan),
    )
    result = _apply_hierarchy(group)

    assert len(result) == 2, "an undecidable group must survive to duplicates/"


def test_a_group_matching_no_reference_row_is_returned_whole():
    """The designed-for failure: loud, and reported through duplicates/."""
    group = _group(
        _row("ENST1", "nan", np.nan, np.nan, np.nan),
        _row("ENST2", "nan", np.nan, np.nan, np.nan),
    )
    group["_merge"] = "left_only"
    result = _apply_hierarchy(group)

    assert len(result) == 2


# ── The invariant the whole module rests on ───────────────────────────────────

@pytest.mark.parametrize(
    "cds_lengths,transcript_lengths,canonical,ccds",
    [
        ((np.nan, np.nan), (np.nan, np.nan), "1.0", np.nan),
        ((np.nan, np.nan), (3000, 2000), "1.0", np.nan),
        ((np.nan, np.nan), (np.nan, np.nan), "nan", "CCDS"),
        ((np.nan, np.nan), (3000, 2000), "nan", "CCDS"),
        ((np.nan, np.nan), (np.nan, np.nan), "nan", np.nan),
        ((900, 1200), (2000, 3000), "1.0", "CCDS"),
        ((900, np.nan), (2000, np.nan), "nan", "CCDS"),
    ],
)
def test_the_hierarchy_never_returns_an_empty_group(
    cds_lengths, transcript_lengths, canonical, ccds
):
    """No combination of missing reference data may delete a variant."""
    group = _group(
        _row("ENST1", canonical, cds_lengths[0], transcript_lengths[0], ccds),
        _row("ENST2", canonical, cds_lengths[1], transcript_lengths[1], ccds),
    )
    assert len(_apply_hierarchy(group)) >= 1
