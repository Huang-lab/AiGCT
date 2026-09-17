"""Gene-level balancing: the balance itself, and the coverage it reports.

`balance_by_gene` prints how many sampled variants lack a score from the
priority VEP. That line is the only readout of how well priority-first sampling
worked, and it is easy to make it always say "none". These tests cover the
sampling invariants and pin the reported count to the real one.

    pytest curation_benchmark/tests
"""
import io
import os
import sys
from contextlib import redirect_stdout

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation.balance import balance_by_gene  # noqa: E402

PRIORITY = "AlphaMissense_rankscore"
GENE = "Ensembl_geneid"


def _variants(gene, n_pos, n_neg, pos_scored=0, neg_scored=0):
    """One gene's records, `*_scored` of each class carrying a priority score."""
    rows = []
    for label, total, scored in ((1, n_pos, pos_scored), (0, n_neg, neg_scored)):
        for i in range(total):
            rows.append({
                GENE: gene,
                "BINARY_LABEL": label,
                "clinvar_review": 2,
                PRIORITY: 0.5 if i < scored else np.nan,
            })
    return rows


def _run(df, tmp_path):
    out = tmp_path / "balanced.csv"
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        result = balance_by_gene(df, output=str(out))
    return result, buffer.getvalue()


# ── The reported coverage ─────────────────────────────────────────────────────

def test_the_missing_priority_score_count_is_the_real_one(tmp_path):
    """Four of the eight sampled variants have no priority score."""
    df = pd.DataFrame(_variants("G1", 4, 4, pos_scored=2, neg_scored=2))
    result, output = _run(df, tmp_path)

    assert len(result) == 8
    assert result[PRIORITY].isna().sum() == 4
    assert f"without a {PRIORITY} score: 4" in output, output


def test_a_fully_scored_set_reports_none_missing(tmp_path):
    df = pd.DataFrame(_variants("G1", 3, 3, pos_scored=3, neg_scored=3))
    result, output = _run(df, tmp_path)

    assert result[PRIORITY].isna().sum() == 0
    assert f"without a {PRIORITY} score: 0" in output, output


def test_a_wholly_unscored_set_reports_every_variant_missing(tmp_path):
    """The case the previous count could not distinguish from "none missing"."""
    df = pd.DataFrame(_variants("G1", 3, 3))
    result, output = _run(df, tmp_path)

    assert len(result) == 6
    assert f"without a {PRIORITY} score: 6" in output, output


# ── The sampling itself ───────────────────────────────────────────────────────

def test_each_gene_is_balanced_to_its_smaller_class(tmp_path):
    df = pd.DataFrame(_variants("G1", 10, 3))
    result, _ = _run(df, tmp_path)

    assert (result["BINARY_LABEL"] == 1).sum() == 3
    assert (result["BINARY_LABEL"] == 0).sum() == 3


def test_a_gene_missing_either_class_is_excluded(tmp_path):
    df = pd.DataFrame(_variants("G1", 4, 4) + _variants("G2", 5, 0))
    result, _ = _run(df, tmp_path)

    assert set(result[GENE]) == {"G1"}


def test_scored_variants_are_preferred_over_unscored_ones(tmp_path):
    """Two of five positives carry a score; both must survive a draw of two."""
    df = pd.DataFrame(_variants("G1", 5, 2, pos_scored=2, neg_scored=2))
    result, _ = _run(df, tmp_path)

    positives = result[result["BINARY_LABEL"] == 1]
    assert len(positives) == 2
    assert positives[PRIORITY].notna().all()


def test_only_one_to_four_star_records_are_eligible(tmp_path):
    df = pd.DataFrame(_variants("G1", 4, 4))
    df.loc[df.index[:2], "clinvar_review"] = 0
    result, _ = _run(df, tmp_path)

    assert (result["clinvar_review"] == 0).sum() == 0
    assert (result["BINARY_LABEL"] == 1).sum() == 2


def test_the_sample_is_reproducible(tmp_path):
    df = pd.DataFrame(_variants("G1", 12, 12, pos_scored=4, neg_scored=4))
    first, _ = _run(df.copy(), tmp_path)
    second, _ = _run(df.copy(), tmp_path)

    pd.testing.assert_frame_equal(first, second)
