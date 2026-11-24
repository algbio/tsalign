# Python Bindings for [`lib_tsalign`](https://crates.io/crates/lib_tsalign)

[![PyPI](https://img.shields.io/pypi/v/tsalign)](https://pypi.org/project/tsalign/)

Python bindings for the template switch aligner.  Aligns two DNA sequences
while detecting template switches — short-range translocations where a query
region is copied from (or aligns to) a different location, possibly on the
reverse complement strand.

## Installation

```bash
pip install tsalign
```

## Quick start

```python
import tsalign

result = tsalign.align("ACGTACGT", "ACGACGT")
print(result.cigar())   # compact alignment string
print(result.stats())   # cost, duration, node counts, …
```

## Aligner options

Create an `Aligner` once and reuse it for many sequences:

```python
aligner = tsalign.Aligner(
    min_length_strategy="preprocess_lookahead",  # default: "lookahead"
    total_length_strategy="maximise",            # default: "maximise"
    no_ts=False,                                 # set True for plain gap-affine
)

result = aligner.align("ACGTACGT", "ACGACGT")
```

## Custom cost configuration

Costs are specified in `.tsa` format.  Use `sample_tsa_config/config.tsa`
as a starting point and consult the [main repository README](../README.md#modifying-the-cost-function) for a description
of each parameter.

```python
aligner = tsalign.Aligner(costs_file="sample_tsa_config/config.tsa")
result = aligner.align("ACGTACGT", "ACGACGT")
```

You can also pass the cost string directly:

```python
with open("my_costs.tsa") as f:
    cost_str = f.read()
aligner = tsalign.Aligner(costs=cost_str)
```

## Restricting the alignment range

Use `AlignmentRange` to align only a window of the input sequences:

```python
from tsalign import Aligner, AlignmentRange

aligner = Aligner()
result = aligner.align(
    "NNNACGTACGTNNN",
    "ACGACGT",
    range=AlignmentRange(reference_start=3, reference_end=11),
)
print(result.cigar())
```

Individual start/limit keyword arguments are also accepted when `range` is
not provided:

```python
result = aligner.align(
    "NNNACGTACGTNNN",
    "ACGACGT",
    reference_start=3,
    reference_limit=11,
)
```

## Skipping gap characters

Characters listed in `gap_characters` are removed from both sequences
before aligning.  Range coordinates refer to the sequences without these characters.

```python
result = aligner.align("ACG-TACGT", "AC-GACGT", gap_characters="-")
```

## Working with alignment operations

`alignment.alignments()` returns a typed list of `(count, op)` pairs:

```python
from tsalign import align, TemplateSwitchEntranceOp, TemplateSwitchExitOp

result = align(reference, query)
for count, op in result.alignments():
    if isinstance(op, TemplateSwitchEntranceOp):
        print(f"Template switch: {op.direction}, descendant={op.descendant}, offset={op.first_offset}")
    elif isinstance(op, TemplateSwitchExitOp):
        print(f"Exit, anti-descendant gap: {op.anti_descendant_gap}")
    else:
        # SimpleAlignmentOp — a basic edit in the primary or secondary track
        print(f"{count}x {op.kind}")
```

## Visualisation

```python
result = tsalign.align(reference, query)
result.viz_template_switches()   # prints ASCII art to stdout
```

## Limiting search resources

```python
result = aligner.align(
    reference,
    query,
    cost_limit=100,       # return None if every alignment costs more than this (inclusive limit)
    memory_limit=500_000, # return None if memory exceeds roughly this number of bytes
)
if result is None:
    print("No alignment found within limits")
```

## Accepted sequence types

Any object whose `str()` representation is a valid DNA string (ACGTN) is
accepted — including `Bio.Seq`:

```python
from Bio.Seq import Seq
result = tsalign.align(Seq("ACGTACGT"), Seq("ACGACGT"))
```

## tschainalign (chaining aligner)

`ChainAligner` exposes tschainalign, the alignment method
`--alignment-method a-star-chain-ts` of the CLI.  It chains k-mer anchors
(exact, or with up to `max_anchor_mutations` mutations) and aligns the gaps
between them exactly.  It supports only the characters `A`, `C`, `G` and `T`,
and only cost configurations without flanks.  Its default costs are
[`test_files/config/chainalign/config.tsa`](../test_files/config/chainalign/config.tsa),
built into the module as `tsalign.CHAIN_ALIGNER_DEFAULT_COSTS` (also the CLI's
default for `a-star-chain-ts` without `-c`); the classic default costs of
`Aligner` are not supported by tschainalign.

```python
aligner = tsalign.ChainAligner(
    costs_file="my_costs.tsa",       # optional: or a directory containing config.tsa, or costs="..."
    max_anchor_mutations=1,          # default: 0 (exact anchors)
    cache_directory="tsalign-cache", # optional, shared with the CLI
)
result = aligner.align(reference, query, cost_limit=100, iteration_limit=1000)
if result is not None:
    print(result.cigar(), result.stats()["cost"])
    print(result.chain_statistics())  # anchor counts, pruning, chaining rounds, times
```

The result is the same `Alignment` type as for `Aligner`.  Its statistics
additionally contain `chain_cost` and `chain_statistics`.  `align` returns
`None` if the cost limit (on the chain cost, inclusive) or the iteration limit
(number of chaining rounds) is reached.  It has no `memory_limit`.

The other parameters are named as the CLI options, with the same defaults:
`k` (default: inferred from the sequence lengths), `inexact_chaining_lower_bound`,
`max_chaining_successors`, `max_exact_cost_function_cost`,
`chaining_closed_list`, `chaining_open_list`, `anchor_pruning`
(`"none"`, `"exact"`, `"inexact"`), `inexact_pruning_rules` (a string like the
CLI's `--inexact-pruning-rules`, or a list of rule names),
`exact_overlap_diagonal_distance` and `force_no_preprocessing`.
Invalid or unsupported parameters raise `ValueError`.

### Preprocessing and caching

tschainalign needs lower bounds that are expensive to compute
(preprocessing).  They depend on the costs, `k`, `max_anchor_mutations` and
the length of the longer sequence rounded up to a power of two.  A
`ChainAligner` keeps all lower bounds it has computed in memory, so reuse
one aligner for many alignments, and fix `k` if the sequence lengths vary.

With `cache_directory`, the lower bounds are also loaded from and stored in
the same cache files as `tsalign preprocess` and `tsalign align --cache-directory`
use, so Python and the CLI can share them.  `aligner.preprocess(max_length)`
precomputes them for all sequences up to `max_length`, like
`tsalign preprocess`.  With `force_no_preprocessing=True`, `align` raises
`FileNotFoundError` instead of preprocessing.

The GIL is released while preprocessing and aligning, so a `ChainAligner`
can be used from several threads.

`tsalign.align(reference, query, method="a-star-chain-ts")`
creates a temporary `ChainAligner`, which preprocesses on every call unless a
`cache_directory` is given.
