"""dbNSFP annotation: reading coordinate lists, and recording what did not match.

These run against a three-row stand-in for a dbNSFP chromosome file, so they
need neither the real release nor the reference tables.

    pytest curation_benchmark/tests
"""
import gzip
import io
import os
import sys
from contextlib import redirect_stdout

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation import dbnsfp  # noqa: E402
from curation.columns import EXTRACT_COLUMNS  # noqa: E402


# ── read_annotation ──────────────────────────────────────────────────────────

def _write(path, text):
    with open(path, "w", newline="") as f:
        f.write(text)
    return str(path)


def test_reads_a_plain_list(tmp_path):
    p = _write(tmp_path / "a.txt", "1 100 A G\n2 200 C T\n")
    a = dbnsfp.read_annotation(p)
    assert a.values.tolist() == [["1", "100", "A", "G"], ["2", "200", "C", "T"]]


def test_tolerates_trailing_space_and_crlf(tmp_path):
    """The format of two committed coordinate lists. A single-space split turns
    the trailing space into a fifth field and shifts every column left."""
    p = _write(tmp_path / "a.txt", "1 100 A G \r\n2 200 C T \r\n")
    a = dbnsfp.read_annotation(p)
    assert a.values.tolist() == [["1", "100", "A", "G"], ["2", "200", "C", "T"]]
    assert a["Alt"].notna().all()


def test_decodes_iupac_ambiguity_codes(tmp_path):
    """Two committed lists write heterozygous alternates as ambiguity codes."""
    p = _write(tmp_path / "a.txt", "1 100 G R\n1 200 C T\n")
    a = dbnsfp.read_annotation(p)
    assert a["Alt"].tolist() == ["A", "T"]


def test_rejects_a_malformed_line(tmp_path):
    p = _write(tmp_path / "a.txt", "1 100 A G\n2 200 C\n")
    with pytest.raises(ValueError, match="malformed"):
        dbnsfp.read_annotation(p)


# ── extract: unmatched variants are returned, not dropped ────────────────────

def _fake_dbnsfp(tmp_path, monkeypatch):
    """A one-chromosome dbNSFP with two hg38 rows; one also carries hg19 coords."""
    rows = []
    for pos, hg19pos, ref, alt in (("100", "90", "A", "G"), ("200", ".", "C", "T")):
        r = {c: "." for c in EXTRACT_COLUMNS}
        r.update({"#chr": "1", "pos(1-based)": pos, "ref": ref, "alt": alt,
                  "hg19_chr": "1" if hg19pos != "." else ".", "hg19_pos(1-based)": hg19pos,
                  "Ensembl_transcriptid": "ENST1", "VEP_canonical": "YES"})
        rows.append(r)
    path = tmp_path / "chr1.gz"
    with gzip.open(path, "wt") as f:
        pd.DataFrame(rows, columns=EXTRACT_COLUMNS).to_csv(f, sep="\t", index=False)
    monkeypatch.setattr(dbnsfp, "dbnsfp_chromosome_file", lambda chrom: str(path))


def test_hg38_join_separates_matched_from_unmatched(tmp_path, monkeypatch):
    _fake_dbnsfp(tmp_path, monkeypatch)
    ann = pd.DataFrame([["1", "100", "A", "G"], ["1", "200", "C", "T"], ["1", "300", "G", "A"]],
                       columns=dbnsfp.ANNOTATION_COLUMNS)
    with redirect_stdout(io.StringIO()):
        matched, unmatched = dbnsfp.extract(ann, "hg38")
    assert len(matched) == 2
    assert unmatched.values.tolist() == [["1", "300", "G", "A"]]


def test_hg19_join_uses_hg19_position_and_hg38_alleles(tmp_path, monkeypatch):
    _fake_dbnsfp(tmp_path, monkeypatch)
    ann = pd.DataFrame([["1", "90", "A", "G"], ["1", "200", "C", "T"]],
                       columns=dbnsfp.ANNOTATION_COLUMNS)
    with redirect_stdout(io.StringIO()):
        matched, unmatched = dbnsfp.extract(ann, "hg19")
    assert len(matched) == 1                       # hg19 pos 90 -> row 1
    assert unmatched.values.tolist() == [["1", "200", "C", "T"]]   # 200 is an hg38 position


def test_unsupported_assembly_is_refused():
    with pytest.raises(ValueError, match="assembly"):
        dbnsfp.extract(pd.DataFrame(columns=dbnsfp.ANNOTATION_COLUMNS), "hg18")


# ── the record of what did not match ─────────────────────────────────────────

def test_unmatched_file_and_match_log_are_written(tmp_path, monkeypatch):
    monkeypatch.setattr(dbnsfp, "output_path",
                        lambda kind, name=None: str(tmp_path / kind / name))
    (tmp_path / "unmatched").mkdir()
    ann = pd.DataFrame([["1", "100", "A", "G"], ["1", "300", "G", "A"]],
                       columns=dbnsfp.ANNOTATION_COLUMNS)
    unmatched = ann.iloc[[1]]
    with redirect_stdout(io.StringIO()):
        dbnsfp._record_unmatched("ds", "hg38", str(tmp_path / "in.txt"), ann, unmatched)

    assert (tmp_path / "unmatched" / "ds.txt").read_text() == "1 300 G A\n"
    log = pd.read_csv(tmp_path / "unmatched" / "match_rates.csv")
    assert log.to_dict("records") == [{
        "dataset": "ds", "assembly": "hg38", "n_input": 2,
        "n_matched": 1, "n_unmatched": 1, "match_rate": 0.5,
    }]


def test_rerunning_a_dataset_replaces_its_log_row(tmp_path, monkeypatch):
    monkeypatch.setattr(dbnsfp, "output_path",
                        lambda kind, name=None: str(tmp_path / kind / name))
    (tmp_path / "unmatched").mkdir()
    ann = pd.DataFrame([["1", "100", "A", "G"]], columns=dbnsfp.ANNOTATION_COLUMNS)
    with redirect_stdout(io.StringIO()):
        dbnsfp._record_unmatched("ds", "hg38", str(tmp_path / "in.txt"), ann, ann.iloc[:0])
        dbnsfp._record_unmatched("ds", "hg38", str(tmp_path / "in.txt"), ann, ann)
    log = pd.read_csv(tmp_path / "unmatched" / "match_rates.csv")
    assert len(log) == 1 and log["n_unmatched"].iloc[0] == 1


def test_output_paths_may_not_be_the_input(tmp_path, monkeypatch):
    """The failure that destroyed MSK_passenger_annotation.txt."""
    p = _write(tmp_path / "in.txt", "1 100 A G\n")
    with pytest.raises(ValueError, match="overwrite"):
        dbnsfp.annotate_dataset(p, p, "hg38", "ds")

    monkeypatch.setattr(dbnsfp, "output_path", lambda kind, name=None: p)
    ann = dbnsfp.read_annotation(p)
    with pytest.raises(ValueError, match="overwrite"):
        dbnsfp._record_unmatched("ds", "hg38", p, ann, ann)
