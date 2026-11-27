"""Tests of tschainalign in the Python bindings (``tsalign.ChainAligner``, issue #269).

The bindings are built in debug mode by ``maturin develop``, and preprocessing is expensive,
so the inputs are small, and the sequences up to length 16 share a single preprocessing per configuration.

The comparison with the CLI builds the ``tsalign`` binary with cargo (debug profile).
Set ``TSALIGN_CLI`` to the path of a ``tsalign`` binary to use that one instead.
"""

import json
import os
import pathlib
import re
import shutil
import subprocess

import pytest

import tsalign
from tsalign import AlignmentRange, ChainAligner

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
TEST_FILES = REPOSITORY / "test_files"
# The standard config of tschainalign: only reverse template switches.
CHAIN_CONFIG = TEST_FILES / "config" / "chainalign" / "config.tsa"
# A config that tschainalign does not support (e.g. its length costs).
UNSUPPORTED_CONFIG = REPOSITORY / "sample_tsa_config" / "config.tsa"

# Inputs of up to 100 characters, so that the debug build preprocesses them quickly.
CLI_COMPARISON_FILES = [
    "twin_10_ts_only.fa",
    "twin_10_ts.fa",
    "twin_10_0.1.fa",
    "twin_in_place_inversion.fa",
    "twin_100_0.01.fa",
]

# twin_10_ts_only.fa: the query is the reverse complement of the reference.
TS_REFERENCE = "GGAGGAGGTT"
TS_QUERY = "AACCTCCTCC"
# The cost of twin_10_ts_only.fa with the chainalign config, as computed by the CLI.
TS_COST = 2


def read_pair_fasta(path):
    """Returns the reference and the query of a pair fasta file."""
    records = []
    for line in pathlib.Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            records.append([])
        elif line:
            records[-1].append(line)
    assert len(records) == 2, path
    return tuple("".join(record).upper() for record in records)


@pytest.fixture(scope="module")
def exact_aligner():
    """An aligner with exact anchors that keeps its lower bounds between the tests."""
    return ChainAligner(costs_file=CHAIN_CONFIG, k=4)


@pytest.fixture(scope="module")
def inexact_aligner():
    return ChainAligner(costs_file=CHAIN_CONFIG, k=5, max_anchor_mutations=1)


@pytest.fixture(scope="module")
def tsalign_cli(tmp_path_factory):
    """Path to the tsalign CLI binary."""
    if "TSALIGN_CLI" in os.environ:
        return pathlib.Path(os.environ["TSALIGN_CLI"])
    if shutil.which("cargo") is None:
        pytest.skip("cargo is not available to build the tsalign CLI")
    output = subprocess.run(
        ["cargo", "build", "-p", "tsalign", "--bin", "tsalign", "--message-format=json"],
        cwd=REPOSITORY,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    for line in output.decode().splitlines():
        message = json.loads(line)
        if message.get("reason") == "compiler-artifact" and message.get("executable"):
            if message["target"]["name"] == "tsalign":
                return pathlib.Path(message["executable"])
    raise AssertionError("cargo did not report the tsalign executable")


def cli_align(cli, pair_fasta, tmp_path, *args):
    """Runs the CLI with tschainalign and returns the cost of the alignment."""
    output = tmp_path / "statistics.toml"
    subprocess.run(
        [
            str(cli),
            "align",
            "--alignment-method",
            "a-star-chain-ts",
            "-c",
            str(CHAIN_CONFIG.parent),
            "--alphabet",
            "dna",
            "-l",
            "warn",
            "-p",
            str(pair_fasta),
            "-o",
            str(output),
            *args,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    match = re.search(r"^cost = ([0-9.]+)$", output.read_text(), re.MULTILINE)
    assert match is not None, output.read_text()
    return float(match.group(1))


def test_align_finds_template_switch(exact_aligner):
    result = exact_aligner.align(TS_REFERENCE, TS_QUERY, reference_name="ref1", query_name="qry1")
    assert result is not None
    stats = result.stats()
    assert stats["cost"] == TS_COST
    assert stats["chain_cost"] == TS_COST
    assert stats["template_switch_amount"] == 1
    assert stats["sequences"]["reference_name"] == "ref1"
    assert stats["sequences"]["query_name"] == "qry1"
    assert result.cigar()
    ops = result.alignments()
    assert any(isinstance(op, tsalign.TemplateSwitchEntranceOp) for _, op in ops)
    assert any(isinstance(op, tsalign.TemplateSwitchExitOp) for _, op in ops)


def test_chain_statistics(exact_aligner, inexact_aligner):
    result = exact_aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    chain_statistics = result.chain_statistics()
    assert chain_statistics == result.stats()["chain_statistics"]
    assert chain_statistics["chaining_rounds"] >= 1
    assert chain_statistics["anchor_pruning"] == "exact"

    # The classic aligner has no chain statistics.
    classic = tsalign.align(TS_REFERENCE, TS_QUERY)
    assert classic is not None
    assert classic.chain_statistics() is None


def test_viz_template_switches(exact_aligner, capfd):
    result = exact_aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    result.viz_template_switches()
    assert capfd.readouterr().out


@pytest.mark.parametrize("max_anchor_mutations", [0, 1])
@pytest.mark.parametrize("file_name", CLI_COMPARISON_FILES)
def test_same_cost_as_cli(tsalign_cli, tmp_path, file_name, max_anchor_mutations):
    pair_fasta = TEST_FILES / file_name
    # The CLI's align does not create the cache directory.
    cache_directory = tmp_path / "cache"
    cache_directory.mkdir()
    cli_cost = cli_align(
        tsalign_cli,
        pair_fasta,
        tmp_path,
        "--max-anchor-mutations",
        str(max_anchor_mutations),
        "--cache-directory",
        str(cache_directory),
    )

    # k is inferred from the sequence lengths as in the CLI.
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, max_anchor_mutations=max_anchor_mutations)
    reference, query = read_pair_fasta(pair_fasta)
    result = aligner.align(reference, query)
    assert result is not None
    assert result.stats()["cost"] == cli_cost


def test_python_uses_cli_cache(tsalign_cli, tmp_path):
    cache_directory = tmp_path / "cache"
    cache_directory.mkdir()
    pair_fasta = TEST_FILES / "twin_10_ts_only.fa"
    cli_cost = cli_align(
        tsalign_cli, pair_fasta, tmp_path, "-k", "4", "--cache-directory", str(cache_directory)
    )
    assert list(cache_directory.glob("*.tsc"))

    aligner = ChainAligner(
        costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory, force_no_preprocessing=True
    )
    result = aligner.align(*read_pair_fasta(pair_fasta))
    assert result is not None
    assert result.stats()["cost"] == cli_cost


def test_python_uses_cli_preprocess_cache(tsalign_cli, tmp_path):
    cache_directory = tmp_path / "cache"
    subprocess.run(
        [
            str(tsalign_cli),
            "preprocess",
            "-l",
            "warn",
            "-c",
            str(CHAIN_CONFIG.parent),
            "--cache-directory",
            str(cache_directory),
            "-k",
            "4",
            "16",
        ],
        check=True,
    )
    aligner = ChainAligner(
        costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory, force_no_preprocessing=True
    )
    result = aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    assert result.stats()["cost"] == TS_COST


def test_cli_uses_python_cache(tsalign_cli, tmp_path):
    cache_directory = tmp_path / "cache"
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory)
    result = aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    cache_files = list(cache_directory.glob("*.tsc"))
    assert len(cache_files) == 1
    # The same file name scheme as the CLI: <cost hash>-<k>-<max anchor mutations>-<lower bound>-<max_n>.tsc
    assert re.fullmatch(r"[0-9a-f]{40}-4-0-mmr-16\.tsc", cache_files[0].name)

    cli_cost = cli_align(
        tsalign_cli,
        TEST_FILES / "twin_10_ts_only.fa",
        tmp_path,
        "-k",
        "4",
        "--cache-directory",
        str(cache_directory),
        "--force-no-preprocessing",
    )
    assert cli_cost == result.stats()["cost"]


def test_preprocess_writes_cache_like_cli(tmp_path):
    cache_directory = tmp_path / "cache"
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory)
    aligner.preprocess(16)
    # Like `tsalign preprocess -k 4 16`, for max_n = 16, 8, 4, 2, 1.
    names = sorted(path.name.split("-", 1)[1] for path in cache_directory.glob("*.tsc"))
    assert names == sorted(f"4-0-mmr-{max_n}.tsc" for max_n in [1, 2, 4, 8, 16])


def test_lower_bounds_are_kept_in_memory(tmp_path):
    cache_directory = tmp_path / "cache"
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory)
    first = aligner.align(TS_REFERENCE, TS_QUERY)
    assert first is not None

    # Without the cache files, the aligner must reuse the lower bounds it keeps in memory.
    shutil.rmtree(cache_directory)
    second = aligner.align(TS_REFERENCE, TS_QUERY)
    assert second is not None
    assert second.stats()["cost"] == first.stats()["cost"]
    assert not cache_directory.exists()


def test_preprocess_in_memory():
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, k=4, force_no_preprocessing=True)
    # Nothing is preprocessed and there is no cache directory.
    with pytest.raises(FileNotFoundError):
        aligner.align(TS_REFERENCE, TS_QUERY)
    aligner.preprocess(16)
    result = aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    assert result.stats()["cost"] == TS_COST


def test_force_no_preprocessing_without_cache_file(tmp_path):
    aligner = ChainAligner(
        costs_file=CHAIN_CONFIG, k=4, cache_directory=tmp_path, force_no_preprocessing=True
    )
    with pytest.raises(FileNotFoundError, match="force_no_preprocessing"):
        aligner.align(TS_REFERENCE, TS_QUERY)


def test_cache_file_with_other_parameters_is_rejected(tmp_path):
    ChainAligner(costs_file=CHAIN_CONFIG, k=4, cache_directory=tmp_path).align(TS_REFERENCE, TS_QUERY)
    # Rename the cache file to the name of k = 5, so that it does not match its content.
    (cache_file,) = tmp_path.glob("*.tsc")
    cache_file.rename(cache_file.with_name(cache_file.name.replace("-4-0-", "-5-0-")))
    with pytest.raises(ValueError, match="different costs, k or max_anchor_mutations"):
        ChainAligner(costs_file=CHAIN_CONFIG, k=5, cache_directory=tmp_path).align(
            TS_REFERENCE, TS_QUERY
        )


def test_cost_limit(exact_aligner):
    # The cost limit is inclusive.
    assert exact_aligner.align(TS_REFERENCE, TS_QUERY, cost_limit=TS_COST + 1) is not None
    result = exact_aligner.align(TS_REFERENCE, TS_QUERY, cost_limit=TS_COST)
    assert result is not None
    assert result.stats()["cost"] == TS_COST
    assert exact_aligner.align(TS_REFERENCE, TS_QUERY, cost_limit=TS_COST - 1) is None
    # Huge limits are no limits.
    assert exact_aligner.align(TS_REFERENCE, TS_QUERY, cost_limit=2**40) is not None


def test_iteration_limit(exact_aligner):
    assert exact_aligner.align(TS_REFERENCE, TS_QUERY, iteration_limit=0) is None
    result = exact_aligner.align(TS_REFERENCE, TS_QUERY, iteration_limit=1_000_000)
    assert result is not None
    assert result.stats()["cost"] == TS_COST


def test_gap_characters_and_range(exact_aligner):
    result = exact_aligner.align("GGA-GGAGGTT", "AA-CCTCCTCC-", gap_characters="-")
    assert result is not None
    assert result.stats()["cost"] == TS_COST

    reference = "ACGT" + TS_REFERENCE
    with_range = exact_aligner.align(
        reference, TS_QUERY, range=AlignmentRange(reference_start=4, reference_end=len(reference))
    )
    with_limits = exact_aligner.align(reference, TS_QUERY, reference_start=4)
    assert with_range is not None and with_limits is not None
    assert with_range.stats()["reference_offset"] == 4
    assert with_range.cigar() == with_limits.cigar()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"anchor_pruning": "none"},
        {"anchor_pruning": "exact"},
        {"anchor_pruning": "inexact"},
        {"anchor_pruning": "inexact", "inexact_pruning_rules": "all,-secondary-exact-overlap"},
        {"anchor_pruning": "inexact", "inexact_pruning_rules": ["primary-drop-guards", "secondary"]},
        {"anchor_pruning": "inexact", "exact_overlap_diagonal_distance": 2},
        {"inexact_chaining_lower_bound": "max-match-run"},
        {
            "max_chaining_successors": 3,
            "max_exact_cost_function_cost": 4,
            "chaining_closed_list": "fx-hash-map",
            "chaining_open_list": "std-heap",
        },
    ],
)
def test_parameters(kwargs):
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, k=5, max_anchor_mutations=1, **kwargs)
    result = aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    chain_statistics = result.chain_statistics()
    assert chain_statistics["anchor_pruning"] == kwargs.get("anchor_pruning", "exact")
    if kwargs.get("anchor_pruning") == "inexact":
        rules = chain_statistics["inexact_pruning_rules"].split(",")
        if kwargs.get("inexact_pruning_rules") == "all,-secondary-exact-overlap":
            assert "secondary-exact-overlap" not in rules
            assert "primary-exact-overlap" in rules
        elif "inexact_pruning_rules" in kwargs:
            assert "primary-drop-guards" in rules
            assert "secondary-exact-overlap" in rules
            assert "primary-exact-overlap" not in rules


def test_inexact_anchors(inexact_aligner):
    result = inexact_aligner.align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    assert result.stats()["cost"] == TS_COST


def test_module_level_align():
    result = tsalign.align(
        TS_REFERENCE, TS_QUERY, method="a-star-chain-ts", costs_file=CHAIN_CONFIG, k=4, cost_limit=100
    )
    assert result is not None
    assert result.stats()["cost"] == TS_COST
    assert result.chain_statistics() is not None
    # The classic method is the default, and can be selected explicitly.
    classic = tsalign.align("ACGTACGT", "ACGACGT", method="a-star-template-switch")
    assert classic is not None
    assert classic.cigar() == tsalign.align("ACGTACGT", "ACGACGT").cigar()


def test_costs_string_and_config_directory():
    costs = CHAIN_CONFIG.read_text()
    from_string = ChainAligner(costs=costs, k=4).align(TS_REFERENCE, TS_QUERY)
    # Like the CLI's --configuration-directory, a directory containing config.tsa is accepted.
    from_directory = ChainAligner(costs_file=CHAIN_CONFIG.parent, k=4).align(TS_REFERENCE, TS_QUERY)
    assert from_string is not None and from_directory is not None
    assert from_string.stats()["cost"] == from_directory.stats()["cost"] == TS_COST


# Invalid parameters raise ValueError instead of a Rust panic (pyo3_runtime.PanicException), as in the CLI.
@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"k": 0}, "max_anchor_mutations = 0 must be less than k = 0"),
        ({"k": 65}, "too large for exact anchors"),
        ({"k": -1}, "k"),
        ({"k": 3, "max_anchor_mutations": 3}, "must be less than k"),
        ({"k": 300, "max_anchor_mutations": 1}, "too large for inexact anchors"),
        ({"max_anchor_mutations": -1}, "max_anchor_mutations"),
        ({"max_anchor_mutations": 256}, "max_anchor_mutations"),
        ({"anchor_pruning": "foo"}, "anchor_pruning"),
        ({"max_anchor_mutations": 1, "anchor_pruning": "foo"}, "anchor_pruning"),
        ({"chaining_closed_list": "foo"}, "chaining_closed_list"),
        ({"chaining_open_list": "foo"}, "chaining_open_list"),
        ({"max_chaining_successors": -1}, "max_chaining_successors"),
        ({"max_exact_cost_function_cost": -1}, "max_exact_cost_function_cost"),
        # Options for inexact anchors are rejected with exact anchors, as in the CLI.
        ({"anchor_pruning": "none"}, "no effect with exact anchors"),
        ({"inexact_chaining_lower_bound": "max-match-run"}, "no effect with exact anchors"),
        # The inexact pruning options require anchor_pruning="inexact".
        ({"max_anchor_mutations": 1, "inexact_pruning_rules": "all"}, "require anchor_pruning='inexact'"),
        (
            {"max_anchor_mutations": 1, "anchor_pruning": "exact", "exact_overlap_diagonal_distance": 1},
            "require anchor_pruning='inexact'",
        ),
        (
            {"max_anchor_mutations": 1, "anchor_pruning": "inexact", "inexact_pruning_rules": "foo"},
            "unknown inexact pruning rule 'foo'",
        ),
        (
            {"max_anchor_mutations": 1, "anchor_pruning": "inexact", "exact_overlap_diagonal_distance": -1},
            "exact_overlap_diagonal_distance",
        ),
        # Disabled in the CLI because it is broken.
        ({"max_anchor_mutations": 1, "inexact_chaining_lower_bound": "windowed-min-mutations"}, "disabled"),
        ({"max_anchor_mutations": 1, "inexact_chaining_lower_bound": "foo"}, "inexact_chaining_lower_bound"),
        ({"arena_index": "u16"}, "arena_index"),
        ({"arena_index": "U64"}, "arena_index"),
    ],
)
def test_invalid_parameters_raise_value_error(kwargs, match):
    with pytest.raises(ValueError, match=match):
        ChainAligner(costs_file=CHAIN_CONFIG, **kwargs)


def test_default_costs(tmp_path):
    # The built-in default costs are test_files/config/chainalign/config.tsa.
    assert tsalign.CHAIN_ALIGNER_DEFAULT_COSTS == CHAIN_CONFIG.read_text()
    cache_directory = tmp_path / "cache"
    result = ChainAligner(k=4, cache_directory=cache_directory).align(TS_REFERENCE, TS_QUERY)
    assert result is not None
    assert result.stats()["cost"] == TS_COST
    # The cache file name depends on the costs, so the default aligner wrote the cache of the explicit config.
    explicit = ChainAligner(
        costs_file=CHAIN_CONFIG, k=4, cache_directory=cache_directory, force_no_preprocessing=True
    ).align(TS_REFERENCE, TS_QUERY)
    assert explicit is not None
    assert explicit.stats()["cost"] == result.stats()["cost"]
    assert explicit.cigar() == result.cigar()


def test_module_level_align_default_costs():
    result = tsalign.align(TS_REFERENCE, TS_QUERY, method="a-star-chain-ts", k=4)
    assert result is not None
    assert result.stats()["cost"] == TS_COST


def test_cli_default_config_equals_python_default(tsalign_cli, tmp_path):
    # Without -c, the CLI uses the same built-in default config for tschainalign.
    output = tmp_path / "statistics.toml"
    subprocess.run(
        [
            str(tsalign_cli),
            "align",
            "--alignment-method",
            "a-star-chain-ts",
            "--alphabet",
            "dna",
            "-l",
            "warn",
            "-k",
            "4",
            "--cache-directory",
            str(tmp_path),
            "-p",
            str(TEST_FILES / "twin_10_ts_only.fa"),
            "-o",
            str(output),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        # Not the repository root, which contains the classic default sample_tsa_config.
        cwd=tmp_path,
    )
    match = re.search(r"^cost = ([0-9.]+)$", output.read_text(), re.MULTILINE)
    assert match is not None, output.read_text()
    python_cost = ChainAligner(k=4).align(TS_REFERENCE, TS_QUERY).stats()["cost"]
    assert float(match.group(1)) == python_cost == TS_COST


def test_invalid_costs_raise_value_error():
    with pytest.raises(ValueError, match="at most one"):
        ChainAligner(costs=CHAIN_CONFIG.read_text(), costs_file=CHAIN_CONFIG)
    # tschainalign does not support the sample config (e.g. its length costs).
    with pytest.raises(ValueError, match="not supported by tschainalign"):
        ChainAligner(costs_file=UNSUPPORTED_CONFIG)
    with pytest.raises(ValueError):
        ChainAligner(costs="not a config")


@pytest.mark.parametrize(
    "reference, query, kwargs, match",
    [
        ("ACNTACGT", "ACGT", {}, "supports only the characters A, C, G and T"),
        ("ACGT", "ACGX", {}, "supports only the characters A, C, G and T"),
        ("", "", {}, "at least one non-empty sequence"),
        ("ACGT", "ACGT", {"reference_limit": 5}, "range"),
        ("ACGT", "ACGT", {"query_start": 3, "query_limit": 2}, "range"),
        ("ACGT", "ACGT", {"gap_characters": "ä"}, "ASCII"),
        ("ACGT", "ACGT", {"cost_limit": -1}, "cost_limit"),
        ("ACGT", "ACGT", {"iteration_limit": -1}, "iteration_limit"),
    ],
)
def test_invalid_inputs_raise_value_error(exact_aligner, reference, query, kwargs, match):
    with pytest.raises(ValueError, match=match):
        exact_aligner.align(reference, query, **kwargs)


def test_inferred_k_is_validated():
    # k is inferred from the sequence lengths and max_anchor_mutations, and the inferred k is too large.
    aligner = ChainAligner(costs_file=CHAIN_CONFIG, max_anchor_mutations=200)
    with pytest.raises(ValueError, match="inferred"):
        aligner.align(TS_REFERENCE, TS_QUERY)


def test_invalid_method():
    with pytest.raises(ValueError, match="method"):
        tsalign.align(TS_REFERENCE, TS_QUERY, method="foo")


# The arena index type only limits how many nodes the chaining can store, so it does not change the alignment.
@pytest.mark.parametrize("max_anchor_mutations, k", [(0, 4), (1, 5)])
def test_arena_index_u64_gives_same_alignment(max_anchor_mutations, k):
    reference, query = read_pair_fasta(TEST_FILES / "twin_100_0.01.fa")
    results = [
        ChainAligner(costs_file=CHAIN_CONFIG, k=k, max_anchor_mutations=max_anchor_mutations, **kwargs).align(
            reference, query
        )
        for kwargs in [{}, {"arena_index": "u32"}, {"arena_index": "u64"}]
    ]
    results.append(
        tsalign.align(
            reference,
            query,
            method="a-star-chain-ts",
            costs_file=CHAIN_CONFIG,
            k=k,
            max_anchor_mutations=max_anchor_mutations,
            arena_index="u64",
        )
    )
    assert all(result is not None for result in results)
    assert len({result.cigar() for result in results}) == 1
    assert len({result.stats()["cost"] for result in results}) == 1
    assert len({result.chain_statistics()["chaining_rounds"] for result in results}) == 1
