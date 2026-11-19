"""Resolve IUPAC ambiguity codes in coordinate lists.

Two of the source studies record the alternate allele of a heterozygous call
as an IUPAC ambiguity code -- the pair of bases at the site -- rather than as
the alternate base alone:

    22 19435300 G R      R = A/G, reference is G, so the alternate is A

dbNSFP stores one base, so an exact-match join drops these rows silently. The
decode is unambiguous whenever the reference base is one of the two the code
stands for: the alternate is the other one. Rows where it is not (a code the
reference does not belong to, or an unrecognised letter) are left untouched
and will simply fail to match, which is the honest outcome -- guessing would
invent an allele.

Both affected lists are drawn from the same sample collection as two other
studies in the same task, which record the same variants with ordinary single
bases, so the variants themselves reach the benchmark either way. Decoding
matters for the per-study named filters, which would otherwise be missing the
ambiguity-coded rows.
"""
import pandas as pd

# The two bases each ambiguity code stands for.
IUPAC_PAIRS = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
}

BASES = {"A", "C", "G", "T"}


def decode_allele(ref, alt):
    """The alternate base an ambiguity code stands for, or `alt` unchanged.

    Returns `alt` untouched when it is already a single base, when the code is
    unrecognised, or when `ref` is not one of the two bases the code covers --
    in which case the pair carries no information about which base is the
    alternate.
    """
    if alt in BASES:
        return alt
    pair = IUPAC_PAIRS.get(alt)
    if pair is None or ref not in pair:
        return alt
    return (pair - {ref}).pop()


def decode_annotation(annotation, ref_col="Ref", alt_col="Alt"):
    """Decode ambiguity codes in a coordinate list, and report what changed.

    Returns the frame with `alt_col` decoded, plus the number of rows decoded
    and the number left alone because the code could not be resolved.
    """
    frame = annotation.copy()
    original = frame[alt_col]
    frame[alt_col] = [
        decode_allele(r, a) for r, a in zip(frame[ref_col], original)
    ]
    changed = int((frame[alt_col] != original).sum())
    unresolved = int((~frame[alt_col].isin(BASES)).sum())
    return frame, changed, unresolved
