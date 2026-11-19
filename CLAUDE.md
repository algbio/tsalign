# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Working style

- **Hand off longer work to agents.** If a task takes more than a handful of tool calls
  (investigations across many files, long test runs, implementing and testing a feature, triaging issues, ...),
  delegate it to a subagent (or run it in the background) so that the main session stays responsive
  for the user. Give the agent a self-contained brief (goal, relevant files, constraints such as "do not commit"),
  and review its result before reporting or committing it.
- **Agents that change code work in their own git worktree** (`isolation: "worktree"`), so that concurrent agents
  and the main session do not see each other's (temporary) edits or mix them into commits.
  The agent commits on its worktree branch; review the result and merge or cherry-pick it.
  Agents that only read and investigate can use the shared working tree.
  Commit everything an agent needs before starting it, since a worktree only contains committed code.
  New worktrees are created from `main`, not from the current branch, so tell the agent to first
  `git reset --hard <current branch>` in its (unchanged) worktree. Worktrees live in `.claude/worktrees/` (gitignored).
- Long-running commands (e.g. the end-to-end tests, release builds) should run in the background.
- Temporary files (scratch data, logs, prototypes) go into `tmp/` in the repository root, which is gitignored.
  Never use the system `/tmp` (this is a shared server).
- Commit or push only when asked.
  Commit messages are a short imperative summary line ending with a period, followed by an explanatory body.
- Failing tests that document real bugs stay failing; do not mark them `#[ignore]` or weaken them to make CI green.

## Project overview

Rust workspace for aligning two DNA sequences while allowing for template switch mutations (TSMs, four-point model).

| Crate | Purpose |
| --- | --- |
| `lib_tsalign` | Classic aligner: A* over the full template switch alignment graph (exact). Also the alignment/result types and the `.tsa` config format. |
| `lib_ts_chainalign` | tschainalign: anchor (k-mer) based chaining aligner with precomputed lower bounds; supports inexact anchors (`max_anchor_mutations`). |
| `tsalign` | CLI (`tsalign align`, `tsalign preprocess`, `tsalign show`). |
| `lib_tsshow` | Visualisation of alignments (`tsalign show`). |
| `seed_chain` | Seed chaining utilities. |
| `python_bindings` | Python bindings (published to PyPI as `tsalign`). |
| `tsalign-tests` | Integration and end-to-end tests of the CLI. |

Configs are directories with a `config.tsa` (see `sample_tsa_config/`, `test_files/config/*`).
tschainalign requires the configs to have no flanks and uses only reverse TSMs; `test_files/config/chainalign` is the standard config for it.
Test inputs are pair fasta files in `test_files/`.

## Building and running

```sh
cargo build                                    # debug
cargo run --release -- align -p test_files/twin_a.fa
cargo run --release -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign \
    --alphabet dna -k 9 --max-anchor-mutations 1 --cache-directory tmp/cache -p <pair.fa>
```

- tschainalign caches preprocessed lower bounds as `.tsc` files; always pass `--cache-directory`
  (otherwise they land in the working directory). `tsalign preprocess -k K --max-anchor-mutations M <max_length>`
  precomputes them; it writes nothing for `max_n < k`.
- The `e2e` cargo profile is release with debug assertions and overflow checks.

## Testing

```sh
cargo test -p lib_tsalign                      # run in debug: a should_panic test needs debug assertions
cargo test -p lib_ts_chainalign
cargo test --no-fail-fast -p tsalign-tests     # end-to-end tests, ~10 minutes; run in the background
cargo clippy --workspace --all-features --tests
cargo fmt --check
```

- CI (`.github/workflows/ci.yml`) uses `RUSTFLAGS=-Dwarnings`, so clippy warnings fail the build.
  The end-to-end tests run in the separate job `rust_e2e_test`.
- Use `--no-fail-fast` for multi-binary test runs, otherwise cargo stops after the first failing test binary.
- End-to-end tests of tschainalign (`tsalign-tests/tests/chain_*.rs`, harness in `chain_common/mod.rs`) build the
  `tsalign` binary with the `e2e` profile and check each result against independent oracles
  (recomputed cost with lib_tsalign, classic aligner with and without TSMs, invariance of TSM-free costs over k and
  max anchor mutations, monotonicity, symmetry, and for inexact anchors the chain cost of a run without pruning
  dominated anchors, `--no-anchor-domination-pruning`). Failing inputs are written to `target/tsalign-e2e/failures/`.
  Environment variables: `TSALIGN_E2E_SEED`, `TSALIGN_E2E_RANDOM_CASES`, `TSALIGN_E2E_TIMEOUT_SECS`.
- `lib_ts_chainalign` has two known failing tests in `windowed_min_mutation_chaining` (tie-breaking); the windowed
  min mutations lower bound is disabled in the CLI.

## GitHub

- Repository: `algbio/tsalign`. Feature branches are named `<issue number>-<slug>`.
- The installed `gh` is old (2.4.0): `gh issue view` fails (deprecated Projects classic) and `gh issue close` has no
  `--comment`. Use `gh api repos/algbio/tsalign/issues/<n>` (and `.../comments`) instead.
