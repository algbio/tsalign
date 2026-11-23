"""End-to-end tests of the pure-Python wrapper in ``tsalign/__init__.py``.

These go through ``tsalign.align`` and ``tsalign.Aligner.align`` (not the compiled
extension directly), so they catch mismatches between the wrapper and the Rust
signature of ``_tsalign.Aligner.align``.
"""

import pytest

import tsalign
from tsalign import Aligner, AlignmentRange

REFERENCE = "ACGTACGT"
QUERY = "ACGACGT"
# The alignment of REFERENCE and QUERY with default settings.
EXPECTED_CIGAR = "3=1D4="
EXPECTED_COST = 3


def test_align_defaults():
    result = tsalign.align(REFERENCE, QUERY)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.stats()["cost"] == EXPECTED_COST


def test_aligner_align_defaults():
    result = Aligner().align(REFERENCE, QUERY)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.alignments()


def test_names():
    result = Aligner().align(REFERENCE, QUERY, reference_name="ref1", query_name="qry1")
    assert result is not None
    assert result.stats()["sequences"]["reference_name"] == "ref1"
    assert result.stats()["sequences"]["query_name"] == "qry1"


@pytest.mark.parametrize("align", [tsalign.align, Aligner().align])
def test_cost_limit(align):
    # A cost limit above the alignment cost still finds the alignment ...
    result = align(REFERENCE, QUERY, cost_limit=EXPECTED_COST + 1)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    # ... but one below it does not.
    assert align(REFERENCE, QUERY, cost_limit=EXPECTED_COST - 1) is None


@pytest.mark.parametrize("align", [tsalign.align, Aligner().align])
def test_memory_limit(align):
    # A generous memory limit does not change the result.
    result = align(REFERENCE, QUERY, memory_limit=1_000_000_000)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    # A tiny memory limit makes the search give up.
    assert align(REFERENCE, QUERY, memory_limit=1) is None


def test_range():
    reference = "NNN" + REFERENCE + "NNN"
    with_range = Aligner().align(
        reference,
        QUERY,
        range=AlignmentRange(reference_start=3, reference_end=3 + len(REFERENCE)),
    )
    with_limits = tsalign.align(
        reference,
        QUERY,
        reference_start=3,
        reference_limit=3 + len(REFERENCE),
    )
    assert with_range is not None and with_limits is not None
    assert with_range.cigar() == with_limits.cigar()
    assert with_range.stats()["cost"] == EXPECTED_COST


@pytest.mark.parametrize("align", [tsalign.align, Aligner().align])
def test_gap_characters(align):
    # Gap characters are removed from both sequences before aligning.
    result = align("ACG-TACGT", "AC-GACGT-", gap_characters="-")
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.stats()["cost"] == EXPECTED_COST


def test_gap_characters_with_range():
    # Coordinates refer to the sequences after removing the gap characters.
    result = Aligner().align(
        "NN-NACGT-ACGTNNN",
        "ACGACGT",
        gap_characters="-",
        range=AlignmentRange(reference_start=3, reference_end=3 + len(REFERENCE)),
    )
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.stats()["cost"] == EXPECTED_COST


def test_gap_characters_combined_with_limits():
    # gap_characters, cost_limit and memory_limit must each reach the right parameter.
    aligner = Aligner()

    def run(cost_limit, memory_limit):
        return aligner.align(
            "ACG-TACGT",
            "ACGACGT",
            gap_characters="-",
            cost_limit=cost_limit,
            memory_limit=memory_limit,
        )

    result = run(cost_limit=EXPECTED_COST + 1, memory_limit=1_000_000_000)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert run(cost_limit=EXPECTED_COST - 1, memory_limit=1_000_000_000) is None
    assert run(cost_limit=1_000, memory_limit=1) is None
