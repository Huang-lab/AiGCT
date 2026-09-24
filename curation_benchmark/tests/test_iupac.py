"""IUPAC ambiguity codes in coordinate lists.

    pytest curation_benchmark/tests
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation.iupac import decode_allele, decode_annotation  # noqa: E402


@pytest.mark.parametrize("ref,alt,expected", [
    ("G", "R", "A"),   # R = A/G
    ("A", "R", "G"),
    ("C", "Y", "T"),   # Y = C/T
    ("T", "Y", "C"),
    ("G", "S", "C"),   # S = G/C
    ("A", "W", "T"),   # W = A/T
    ("G", "K", "T"),   # K = G/T
    ("C", "M", "A"),   # M = A/C
])
def test_a_code_resolves_against_the_reference_base(ref, alt, expected):
    assert decode_allele(ref, alt) == expected


@pytest.mark.parametrize("ref,alt", [("A", "G"), ("C", "T")])
def test_a_plain_base_is_left_alone(ref, alt):
    assert decode_allele(ref, alt) == alt


def test_a_code_the_reference_does_not_belong_to_is_left_alone():
    """R is A/G; a reference of C says nothing about which is the alternate."""
    assert decode_allele("C", "R") == "R"


def test_an_unrecognised_letter_is_left_alone():
    assert decode_allele("G", "N") == "N"
    assert decode_allele("G", "-") == "-"


def test_decode_annotation_reports_what_it_changed():
    frame = pd.DataFrame(
        [["1", "10", "G", "R"], ["1", "20", "C", "T"], ["1", "30", "C", "R"]],
        columns=["Chr", "Pos", "Ref", "Alt"],
    )
    out, changed, unresolved = decode_annotation(frame)

    assert out["Alt"].tolist() == ["A", "T", "R"]
    assert changed == 1        # only the first row
    assert unresolved == 1     # the C/R row stays ambiguous
    assert frame["Alt"].tolist() == ["R", "T", "R"]   # input untouched
