# Test matrix plan: inexact anchor pruning (issue #236)

This is a plan for evaluating the anchor pruning modes and the individual inexact pruning rules of tschainalign
on real data in [template-switch-experiments](https://github.com/algbio/template-switch-experiments) (Snakemake, cluster).
It is written against the structure of its `Snakefile` (rules `tsalign`, `tsalign_statistics_files`,
`tsalign_preprocess`, the experiment functions such as `tsalign_bench_statistics`, and the `ari_email` /
`heli_email` datasets).

## What is evaluated

tschainalign has three anchor pruning modes (`--anchor-pruning`, only for `--max-anchor-mutations >= 1`):

| Mode | Meaning |
| --- | --- |
| `none` | Keep all anchors. The reference for costs. |
| `exact` (default) | Remove dominated primary anchors where this provably does not change the chain cost (m = 1 only). |
| `inexact` | Remove primary and secondary anchors with heuristic rules, for any m ≥ 1. May increase the cost and lose template switches. |

The rules of `inexact` are switched with `--inexact-pruning-rules <comma-separated list>` (default `all`).
Each rule exists for primary and secondary anchors; the names are stable:

| Rule | Proposal | What it removes |
| --- | --- | --- |
| `{primary,secondary}-prefix-domination` | P2/S2 | An anchor that is a shorter anchor with the same start plus one column, which is not cheaper than that column. Dominators may be dominated themselves (cascade). |
| `{primary,secondary}-suffix-domination` | P2/S2 mirrored | The same with the same end and the first column. |
| `{primary,secondary}-drop-guards` | P1/S1 | Modifier of the two domination rules: also prune zero-cost anchors, with `gap_extend > gap_open`, and with lower bounds other than max-match-run. Alone it removes nothing. |
| `{primary,secondary}-edge-substitution` | P3/S3 | Gap-free anchors (both sides of length k) with a substitution in the first or last column. |
| `{primary,secondary}-exact-overlap` | P4/S4 | Inexact anchors that overlap an exact anchor (in both sequences) whose diagonal is within `--exact-overlap-diagonal-distance` (default 1) of the anchor's diagonal range. Secondary: same template switch kind. |

Keywords: `all`, `none`, `primary`, `secondary`; a leading `-` removes rules, e.g. `all,-secondary-exact-overlap`.
The precise semantics are in the module documentation of `lib_ts_chainalign/src/anchors/inexact_pruning.rs`.

Notes for the CLI:

* `--anchor-pruning`, `--inexact-pruning-rules` and `--exact-overlap-diagonal-distance` are tschainalign options
  and are rejected for `--max-anchor-mutations 0` (#249). m = 0 cells must not pass them.
* `--inexact-pruning-rules` and `--exact-overlap-diagonal-distance` are rejected unless `--anchor-pruning inexact`.
* `--no-anchor-domination-pruning` still works as a deprecated alias of `--anchor-pruning none`.
* The CLI only offers the `max-match-run` lower bound (`windowed-min-mutations` is disabled), so the lower bound kind
  is not a factor, and the lower-bound part of `drop-guards` has no effect from the CLI. With the bench config
  (`gap_extend <= gap_open`, positive costs), `drop-guards` therefore changes nothing on this data; it can be
  left out of the matrix or kept as a single sanity cell.

## Output that the experiments read

tschainalign writes a `[chain_statistics]` table into the `-o` TOML file (`AlignmentStatistics::chain_statistics`):

```toml
[chain_statistics]
anchor_pruning = "inexact"
inexact_pruning_rules = "primary-prefix-domination,..."   # empty unless inexact
primary_anchors = 1234                 # before pruning
primary_anchors_after_pruning = 456
secondary_anchors = 7890               # sum over the four template switch kinds
secondary_anchors_after_pruning = 1234
anchor_seconds = 0.01                  # anchor generation
pruning_seconds = 0.001
cost_function_seconds = 0.02           # chaining cost function initialisation
chaining_seconds = 0.5                 # A* chaining, all rounds
evaluation_seconds = 0.1               # chain evaluation (gap alignments), all rounds
chaining_rounds = 12
chaining_opened_nodes = 12345
chaining_closed_nodes = 12000

[chain_statistics.pruned_anchors]      # removed per rule; each removed anchor counted once
primary-prefix-domination = 10         # (exact mode: primary-exact-domination)
...

[chain_statistics.matched_anchors]     # matched per rule, independently (overlapping)
...
```

Together with the existing `cost`, `chain_cost`, `template_switch_amount`, `duration_seconds`, and the alignment
itself (to compute identity and template switch positions).

**Required change in `tsalign_statistics_files`:** it nests the raw output under `[statistics]` by rewriting the
table headers `[result]` and `[sequences]` with `sed`. Add the same for the new tables, otherwise they end up at the
top level of `statistics.toml`:

```sh
sed 's/\[result\]/[statistics.result]/g; s/\[sequences\]/[statistics.sequences]/g; s/\[chain_statistics/[statistics.chain_statistics/g'
```

Also add `anchor_pruning`, `inexact_pruning_rules` and `exact_overlap_diagonal_distance` to the `.extra` header
lines written by that rule.

## Factors

| Factor | Levels | Notes |
| --- | --- | --- |
| Dataset | `ari_email` (883 pairs, window `ari`), `heli_email` (443 positive pairs, window `heli`) | Real data, as in `tsalign_bench_statistics`. Config `bench`. |
| k | 6, 7, 8, 9, 10 (full); 7, 9 (reduced) | `tsalign_preliminary_bench_statistics` used 6–10; `bench` uses 4 and 7 for exact anchors. |
| m (`--max-anchor-mutations`) | 1, 2 (+ 0 as reference without pruning flags) | Inexact pruning only matters for m ≥ 1. |
| Pruning mode | `none`, `exact`, `inexact` | `exact` = today's default. |
| Rule set (inexact) | see below | |
| `--exact-overlap-diagonal-distance` | 0, 1, 2 | Only for rule sets that contain an `exact-overlap` rule; a separate small sweep. |
| Lower bound kind | `max-match-run` only | The only one the CLI offers. |
| Other chaining options | fixed: `--max-chaining-successors 10 --max-exact-cost-function-cost 0 --chaining-closed-list fx-hash-map --chaining-open-list linear-heap` | The choices of `bench`. |

### Rule sets

| Id | `--inexact-pruning-rules` | Purpose |
| --- | --- | --- |
| `all` | `all` | Everything. |
| `primary` / `secondary` | `primary` / `secondary` | Effect of pruning only one anchor type. |
| `only-<rule>` (8) | e.g. `primary-exact-overlap` | Each rule alone (without `drop-guards`). |
| `all-but-<rule>` (8) | e.g. `all,-secondary-exact-overlap` | Leave-one-out: contribution of a rule given the others. |
| `domination-<type>` (2) | `<type>-prefix-domination,<type>-suffix-domination` | Both domination rules of one anchor type. |
| `domination-unguarded-<type>` (2, optional) | `<type>-prefix-domination,<type>-suffix-domination,<type>-drop-guards` | Effect of `drop-guards` (expected: none with the bench config). |

That is 21 inexact rule sets (23 with the optional ones), plus `none` and `exact`: 23 pruning cells per (pair, k, m).
For m = 2, `exact` is identical to `none` (the exact rule only prunes for m = 1) and can be skipped.

The rule set id (not the comma-separated list) should be the Snakemake wildcard, mapped to the list in a dict, so
that paths stay short and contain no commas.

## Exact CLI per cell

All cells use the existing `tsalign` rule command line (`align --force-no-preprocessing -p … -o … --cache-directory …
--alignment-method a-star-chain-ts --skip-characters 'N-' --alphabet dna -c config/tsalign/bench -k {k}
--max-anchor-mutations {m} --inexact-chaining-lower-bound max-match-run --max-chaining-successors 10
--max-exact-cost-function-cost 0 --chaining-closed-list fx-hash-map --chaining-open-list linear-heap --rq-ranges …`),
plus the pruning flags:

| Cell | Extra flags |
| --- | --- |
| m = 0 reference | (none; passing any pruning flag is an error) |
| `none` | `--anchor-pruning none` |
| `exact` | `--anchor-pruning exact` (or nothing, it is the default) |
| inexact rule set `R` | `--anchor-pruning inexact --inexact-pruning-rules <list of R>` |
| overlap distance sweep | `--anchor-pruning inexact --inexact-pruning-rules <list of R> --exact-overlap-diagonal-distance {0,1,2}` |

Note that the classic-aligner options of the current command line (`--ts-node-ord-strategy`,
`--ts-min-length-strategy`, `--allow-ts-14-out-of-range`) are rejected for tschainalign since #249 if given
explicitly; the `tsalign` rule has to pass them only for `a-star-template-switch` (that is needed independently of this
matrix). Likewise `--memory-limit` is rejected for tschainalign.

Suggested Snakemake integration: add two wildcards `anchor_pruning` (`none|exact|inexact|default`) and
`pruning_rules` (rule set id, `default` when not inexact) plus `overlap_distance` (`default|0|1|2`) to the last
component of `ALIGNMENT_TSALIGN_SUBDIR`, and a `params` entry that builds the flags (empty for `default`, which also
keeps all existing experiments valid and their paths stable if the default values are inserted by `safe_format`).

## Metrics

Per run (from `statistics.toml` and the `time -v` output already collected):

* **Runtime**: `runtime_raw` (wall clock, `time -v`), `duration_seconds`, and the phases
  `chain_statistics.{anchor,pruning,cost_function,chaining,evaluation}_seconds`, `chaining_rounds`.
* **Memory**: `memory_raw` (max RSS).
* **Anchors**: `primary_anchors`, `secondary_anchors` before and after pruning, `pruned_anchors.*` (per rule),
  `matched_anchors.*`.
* **Cost**: `cost` and `chain_cost`.
* **Template switches**: `template_switch_amount`, and the template switch positions and lengths from the alignment.
* **Alignment identity**: matches / alignment columns, from the alignment.

Derived per (pair, k, m), relative to the baselines:

* Δ chain cost and Δ cost vs. `none` (the main quality measure; the chain cost never decreases with pruning, the
  reported cost can in rare cases, because of issue #221) and vs. `exact`.
* **TS lost**: `none` has ≥ 1 template switch, the cell has fewer; **TS gained** (the opposite); **TS changed**: same
  number, but different positions (compare the 1-2/3-4 jump points, with a tolerance of a few bases).
* Identity difference vs. `none`.
* Speed-up of runtime and memory vs. `none` and vs. `exact`, and the fraction of anchors removed.
* Optionally vs. m = 0 with the same k, to see whether inexact anchors with pruning still find more template switches
  than exact anchors.

How to compare:

* **Baseline for quality**: `none` with the same (pair, k, m). Report the fraction of pairs with Δ chain cost = 0,
  the distribution of Δ cost, and the TS lost/changed counts per rule set.
* **Baseline for speed**: `exact` (today's default) with the same (pair, k, m); for m = 2 `none` (= `exact`).
* For single rules, compare `only-<rule>` with `none`; for interactions, compare `all-but-<rule>` with `all`.
* Plots: per rule set, speed-up (x) vs. fraction of pairs with a cost increase (y), one panel per (dataset, m);
  and boxplots of Δ cost.
* Pairs where `none` hits the time or memory limit should be reported separately (the cost difference is
  then unknown, but they are where pruning matters most).

## Size

Full matrix: 1326 pairs × 5 k × 2 m × 23 cells ≈ **305,000 runs**, plus the m = 0 references (6,630) and the overlap
sweep. That is too large.

Suggested staging:

1. **Screening** (rules): 50 pairs per dataset (seeds 0–49, as in `tsalign_preliminary_bench_statistics`)
   × k ∈ {7, 9} × m ∈ {1, 2} × 23 cells = **9,200 runs** (minus `exact` for m = 2).
   Plus the overlap distance sweep only for `only-primary-exact-overlap`, `only-secondary-exact-overlap` and `all`
   with distances 0 and 2: 100 × 4 × 6 = 2,400 runs.
2. **Confirmation** (full data): all 1326 pairs × k ∈ {6, 7, 8, 9, 10} × m ∈ {1, 2} × the 3–4 most promising rule
   sets from stage 1 plus `none` and `exact` ≈ 1326 × 10 × 6 ≈ **80,000 runs**. If that is still too many, use
   k ∈ {7, 9} (≈ 32,000).

The `none` runs with m = 2 dominate the cost (they have the most anchors); `tsalign_mem_mb` already doubles the
memory per anchor mutation. Pruned cells need less, so they could get the same memory as `none` to keep the
setup simple. The local evaluation in the PR (see below) gives the relative runtimes.

Runtimes per run are short for most pairs (seconds), so the Snakemake/Slurm job overhead dominates. Consider running
all pruning cells of one (pair, k, m) in a single job (a small wrapper script that calls tsalign once per cell and
writes one statistics file per cell), which cuts the number of Slurm jobs by ~20×.

## Local evaluation

The PR that introduces the rules contains a small local evaluation on the repository's test files
(see the PR description), which can guide the choice of rule sets for stage 2.
