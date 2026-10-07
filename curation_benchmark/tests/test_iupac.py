"""IUPAC ambiguity codes in coordinate lists.

    pytest curation_benchmark/tests
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from curation.iupac import (  # noqa: E402
    BASES,
    IUPAC_PAIRS,
    decode_allele,
    decode_annotation,
)


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


@pytest.mark.parametrize("ref", list("ACGT"))
@pytest.mark.parametrize("code", ["B", "D", "H", "V", "N"])
def test_codes_that_cover_more_than_two_bases_are_never_decoded(ref, code):
    """No unique alternate exists, so guessing one would invent an allele."""
    assert decode_allele(ref, code) == code


def test_every_code_in_the_table_stands_for_exactly_two_bases():
    """decode_allele takes `(pair - {ref}).pop()`. That is only deterministic, and
    only correct, for a two-element set: for a larger one it would return an
    arbitrary member, and a different one from run to run under hash
    randomisation. Adding B, D, H, V or N to the table has to fail here first."""
    assert IUPAC_PAIRS
    for code, pair in IUPAC_PAIRS.items():
        assert len(pair) == 2 and pair <= BASES, code
