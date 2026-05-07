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
    chaining_strategy="lower_bound",             # default: "none"
    total_length_strategy="maximise",            # default: "maximise"
    no_ts=False,                                 # set True for plain gap-affine
)

result = aligner.align("ACGTACGT", "ACGACGT")
```

## Custom cost configuration

Costs are specified in `.tsa` format.  Use `sample_tsa_config/config.tsa`
as a starting point and consult the main repository README for a description
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

## Working with alignment operations

`alignment.alignments()` returns a typed list of `(count, op)` pairs:

```python
from tsalign import align, TemplateSwitchEntranceOp, TemplateSwitchExitOp

result = align(reference, query)
for count, op in result.alignments():
    if isinstance(op, TemplateSwitchEntranceOp):
        print(f"Template switch: {op.direction}, primary={op.primary}, offset={op.first_offset}")
    elif isinstance(op, TemplateSwitchExitOp):
        print(f"Exit, anti-primary gap: {op.anti_primary_gap}")
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
    cost_limit=100,       # return None if cost would exceed this
    memory_limit=500_000, # return None if memory exceeds this number of bytes
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
