"""Smoke tests of ``Alignment.viz_template_switches()`` (issue #279).

The visualisation is printed by Rust to the stdout file descriptor, hence it is captured with ``capfd``.
"""

import pathlib

import pytest

import tsalign
from tsalign import Aligner, ChainAligner, TemplateSwitchEntranceOp

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
SMALL_CONFIG = REPOSITORY / "test_files" / "config" / "small" / "config.tsa"
CHAIN_CONFIG = REPOSITORY / "test_files" / "config" / "chainalign" / "config.tsa"

# The query contains the reverse complement of CCG at positions 4..7 of the reference.
INVERSION_REFERENCE = "AAAACCGTTTAGCTAGCTTGCA"
INVERSION_QUERY = "AAAACGGTTTAGCTAGCTTGCA"
# twin_10_ts_only.fa: the query is the reverse complement of the reference.
TS_REFERENCE = "GGAGGAGGTT"
TS_QUERY = "AACCTCCTCC"


def template_switch_count(result):
    return sum(isinstance(op, TemplateSwitchEntranceOp) for _, op in result.alignments())


def visualise(result, capfd):
    capfd.readouterr()
    assert result.viz_template_switches() is None
    return capfd.readouterr().out


@pytest.mark.parametrize(
    "reference, query", [(INVERSION_REFERENCE, INVERSION_QUERY), (TS_REFERENCE, TS_QUERY)]
)
def test_classic_alignment_with_template_switch(reference, query, capfd):
    result = Aligner(costs_file=SMALL_CONFIG).align(reference, query)
    assert result is not None
    count = template_switch_count(result)
    assert count >= 1

    out = visualise(result, capfd)
    assert out.count("Switch process:") == count
    assert out.count("Direction: reverse") == count
    assert "Parent: reference" in out
    assert "Child: query" in out
    # The visualisation does not show a no-TS alignment.
    assert "No-ts alignment:" not in out


def test_alignment_without_template_switches(capfd):
    result = Aligner(no_ts=True).align(INVERSION_REFERENCE, INVERSION_QUERY)
    assert result is not None
    assert template_switch_count(result) == 0
    assert visualise(result, capfd) == ""


def test_alignment_of_equal_sequences(capfd):
    result = tsalign.align("ACGTACGT", "ACGTACGT")
    assert result is not None
    assert visualise(result, capfd) == ""


def test_names_are_shown(capfd):
    result = Aligner(costs_file=SMALL_CONFIG).align(
        TS_REFERENCE, TS_QUERY, reference_name="my_reference", query_name="my_query"
    )
    assert result is not None
    out = visualise(result, capfd)
    assert "my_reference" in out
    assert "my_query" in out


def test_chain_aligner(capfd):
    result = ChainAligner(costs_file=CHAIN_CONFIG, k=4).align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    count = template_switch_count(result)
    assert count >= 1
    out = visualise(result, capfd)
    assert out.count("Switch process:") == count
