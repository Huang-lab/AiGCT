"""Gene-level balancing: the sampling invariants.

`balance_by_gene` draws, within each gene, an equal number of pathogenic and
benign variants uniformly at random under a fixed seed. These tests pin the
invariants that the released benchmark depends on, plus one regression guard:
dbNSFP encodes a missing score as the string ".", and an earlier version of
this module let that value steer the draw without meaning to. The draw must be
independent of any score column.

    pytest curation_benchmark/tests
"""
import io
import os
import sys
from contextlib import redirect_stdout

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation.balance import balance_by_gene  # noqa: E402

GENE = "Ensembl_geneid"
SCORE = "AlphaMissense_rankscore"
GOLDEN = [
    "G1_1_13", "G1_1_4", "G1_1_12", "G1_1_9", "G1_1_7", "G1_1_6", "G1_1_0", "G1_1_11",
    "G1_0_4", "G1_0_5", "G1_0_6", "G1_0_7", "G1_0_1", "G1_0_0", "G1_0_3", "G1_0_2",
]


def _variants(gene, n_pos, n_neg, scored=None):
    """One gene's records. `scored` marks how many of each class carry a score;
    the rest carry dbNSFP's "." placeholder, as the real curation output does."""
    rows = []
    for label, total in ((1, n_pos), (0, n_neg)):
        for i in range(total):
            rows.append({
                "vid": f"{gene}_{label}_{i}",
                GENE: gene,
                "BINARY_LABEL": label,
                "clinvar_review": 2,
                SCORE: "0.5" if scored is None or i < scored else ".",
            })
    return rows


def _run(df, tmp_path):
    out = tmp_path / "balanced.csv"
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        result = balance_by_gene(df, output=str(out))
    return result, buffer.getvalue()


def test_each_gene_is_balanced_to_its_smaller_class(tmp_path):
    df = pd.DataFrame(_variants("G1", 10, 3))
    result, _ = _run(df, tmp_path)

    assert (result["BINARY_LABEL"] == 1).sum() == 3
    assert (result["BINARY_LABEL"] == 0).sum() == 3


def test_a_gene_missing_either_class_is_excluded(tmp_path):
    df = pd.DataFrame(_variants("G1", 4, 4) + _variants("G2", 5, 0))
    result, output = _run(df, tmp_path)

    assert set(result[GENE]) == {"G1"}
    assert "excluded (one label only): 1" in output, output


def test_only_one_to_four_star_records_are_eligible(tmp_path):
    df = pd.DataFrame(_variants("G1", 4, 4))
    df.loc[df.index[:2], "clinvar_review"] = 0
    result, _ = _run(df, tmp_path)

    assert (result["clinvar_review"] == 0).sum() == 0
    assert (result["BINARY_LABEL"] == 1).sum() == 2


def test_the_sample_is_reproducible(tmp_path):
    df = pd.DataFrame(_variants("G1", 12, 12))
    first, _ = _run(df.copy(), tmp_path)
    second, _ = _run(df.copy(), tmp_path)

    pd.testing.assert_frame_equal(first, second)


def test_the_draw_ignores_which_variants_carry_a_score(tmp_path):
    """Scores on some variants and "." on others must not steer the draw.

    14 pathogenic against 8 benign draws the pathogenic class down to 8, so
    which 8 depends on how the draw is made. A priority-first draw takes the
    scored ones first and returns a different 8 than a blind draw does. (A
    check that scores every variant, or none, cannot tell the two apart.)"""
    frame = pd.DataFrame(_variants("G1", 14, 8, scored=5))
    mixed, _ = _run(frame.copy(), tmp_path)
    blind, _ = _run(frame.assign(**{SCORE: "."}), tmp_path)

    assert mixed["vid"].tolist() == blind["vid"].tolist()


def test_each_class_is_a_seeded_uniform_sample(tmp_path):
    """The draw is `DataFrame.sample(n, random_state=18)` on each class, nothing more."""
    frame = pd.DataFrame(_variants("G1", 14, 8, scored=5))
    result, _ = _run(frame, tmp_path)

    pos = frame[frame["BINARY_LABEL"] == 1].sample(8, random_state=18)
    neg = frame[frame["BINARY_LABEL"] == 0].sample(8, random_state=18)
    assert result["vid"].tolist() == pos["vid"].tolist() + neg["vid"].tolist()


def test_the_released_seed_still_selects_the_same_variants(tmp_path):
    """Pins random_state 18 and the per-class draw to literal ids.

    The balanced ClinVar set in the released database was drawn this way; if
    this fails, re-running the module no longer reproduces it. The ids were
    identical on pandas 2.3.3 and 3.0.6."""
    frame = pd.DataFrame(_variants("G1", 14, 8, scored=5))
    result, _ = _run(frame, tmp_path)

    assert result["vid"].tolist() == GOLDEN


def test_gene_stats_and_exclusions_are_written(tmp_path):
    df = pd.DataFrame(_variants("G1", 3, 5) + _variants("G2", 0, 2))
    stats = tmp_path / "stats.csv"
    excluded = tmp_path / "excluded.csv"
    with redirect_stdout(io.StringIO()):
        balance_by_gene(df, output=str(tmp_path / "b.csv"),
                        gene_stats_path=str(stats), excluded_genes_path=str(excluded))

    s = pd.read_csv(stats)
    assert s.to_dict("records") == [{"geneid": "G1", "num_variants": 6}]
    assert pd.read_csv(excluded)["excluded_geneid"].tolist() == ["G2"]
