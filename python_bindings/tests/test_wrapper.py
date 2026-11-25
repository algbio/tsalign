"""End-to-end tests of the pure-Python wrapper in ``tsalign/__init__.py``.

These go through ``tsalign.align`` and ``tsalign.Aligner.align`` (not the compiled
extension directly), so they catch mismatches between the wrapper and the Rust
signature of ``_tsalign.Aligner.align``.
"""

import pathlib

import pytest

import tsalign
from tsalign import Aligner, AlignmentRange

# Not a reverse complement palindrome like ACGTACGT, for which a template switch would cost the same as the deletion.
REFERENCE = "AAGTCCAG"
QUERY = "AAGCCAG"
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


def test_aligner_default_costs():
    # The built-in default costs equal sample_tsa_config/config.tsa, like the CLI's default.
    sample_config = pathlib.Path(__file__).resolve().parents[2] / "sample_tsa_config" / "config.tsa"
    # test_files/twin_10_ts_only.fa: a template switch, so that the template switch costs matter.
    reference = "GGAGGAGGTT"
    query = "AACCTCCTCC"
    default = Aligner().align(reference, query)
    explicit = Aligner(costs_file=sample_config).align(reference, query)
    assert default is not None and explicit is not None
    assert default.cigar() == explicit.cigar()
    assert default.stats()["cost"] == explicit.stats()["cost"]


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
def test_cost_limit_is_inclusive(align):
    # Issue #250: an alignment whose cost equals the cost limit is returned.
    result = align(REFERENCE, QUERY, cost_limit=EXPECTED_COST)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.stats()["cost"] == EXPECTED_COST


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
    result = align("AAG-TCCAG", "AA-GCCAG-", gap_characters="-")
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert result.stats()["cost"] == EXPECTED_COST


def test_gap_characters_with_range():
    # Coordinates refer to the sequences after removing the gap characters.
    result = Aligner().align(
        "NN-NAAGT-CCAGNNN",
        "AAGCCAG",
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
            "AAG-TCCAG",
            "AAGCCAG",
            gap_characters="-",
            cost_limit=cost_limit,
            memory_limit=memory_limit,
        )

    result = run(cost_limit=EXPECTED_COST + 1, memory_limit=1_000_000_000)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR
    assert run(cost_limit=EXPECTED_COST - 1, memory_limit=1_000_000_000) is None
    assert run(cost_limit=1_000, memory_limit=1) is None


def test_lower_bound_chaining_strategy_is_rejected():
    # Issue #248: the lower bound chaining strategy is disabled.
    with pytest.raises(ValueError, match="lower_bound"):
        Aligner(chaining_strategy="lower_bound")
    with pytest.raises(ValueError, match="lower_bound"):
        tsalign.align(REFERENCE, QUERY, chaining_strategy="lower_bound")


def test_none_chaining_strategy_is_accepted():
    result = Aligner(chaining_strategy="none").align(REFERENCE, QUERY)
    assert result is not None
    assert result.cigar() == EXPECTED_CIGAR


# Issue #249: invalid inputs raise ValueError instead of a Rust panic (pyo3_runtime.PanicException).
@pytest.mark.parametrize(
    "kwargs",
    [
        {"reference_limit": len(REFERENCE) + 1},
        {"query_start": len(QUERY) + 1},
        {"reference_start": 5, "reference_limit": 4},
        {"range": AlignmentRange(0, len(REFERENCE), 0, len(QUERY) + 1)},
        # The range refers to the sequences without gap characters.
        {"gap_characters": "T", "reference_limit": len(REFERENCE)},
    ],
)
@pytest.mark.parametrize("align", [tsalign.align, Aligner().align])
def test_invalid_range_raises_value_error(align, kwargs):
    with pytest.raises(ValueError, match="range"):
        align(REFERENCE, QUERY, **kwargs)


@pytest.mark.parametrize("align", [tsalign.align, Aligner().align])
def test_non_alphabet_character_raises_value_error(align):
    with pytest.raises(ValueError, match="Reference contains non-alphabet character"):
        align("AAGTXCAG", QUERY)
    with pytest.raises(ValueError, match="Query contains non-alphabet character"):
        align(REFERENCE, "AAG-CCAG")
    # Gap characters are removed before checking the alphabet.
    assert align(REFERENCE, "AAG-CCAG", gap_characters="-") is not None
