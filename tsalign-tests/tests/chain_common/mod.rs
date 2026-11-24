//! Shared harness for the end-to-end tests of tschainalign (`tsalign align --alignment-method a-star-chain-ts`).
//!
//! Every alignment runs in a separate `tsalign` process built once per test process with the `e2e` cargo profile
//! (release optimisations with debug assertions and overflow checks).
//! Each chaining result is checked against independent oracles:
//!
//! * the process succeeds and reports an alignment;
//! * recomputing the cost of the reported alignment with lib_tsalign gives the reported cost,
//!   and the alignment consumes both sequences completely;
//! * the reported cost is at most the reported chain cost, and equal to it for exact anchors;
//! * the cost is at least the classic A* aligner's optimum with template switches;
//! * for exact anchors, the cost is at most the classic optimum without template switches,
//!   and an alignment without template switches has exactly the cost of the classic optimum without template switches,
//!   regardless of `k`;
//!   for inexact anchors, an alignment without template switches costs at least the classic optimum without them;
//! * for exact anchors, the cost does not increase when the anchors become less restrictive in a way that keeps
//!   every template switch supported (see [`anchors_dominate`]);
//! * swapping reference and query does not change the chain cost;
//! * for inexact anchors, the chain cost equals the chain cost of a run with all anchors,
//!   i.e. without pruning dominated anchors (`--anchor-pruning none`).
//!   Pruning must not change the cost of the optimal chain.
//! * for inexact anchors, a run with the inexact anchor pruning mode (`--anchor-pruning inexact`, all rules)
//!   succeeds with a valid alignment, and its chain cost is at least the chain cost without pruning.
//!   It may be larger, since the inexact pruning may remove anchors that the optimal chain needs.
//!   This is the only check that detects too expensive results for inexact anchors:
//!   when a pruned anchor was required for a template switch, the result is typically the optimal alignment without
//!   template switches, which passes all other checks that are applied to inexact anchors.
//!
//! tschainalign finds alignments that are optimal with regard to the anchors, not globally optimal alignments.
//! For inexact anchors (`max_anchor_mutations >= 1`), the chain cost is the sum of the costs of the gaps and anchors,
//! which are aligned independently, while the reported cost is the cost of the final alignment as a whole,
//! which may be lower (adjacent gaps are merged, and anchors are re-aligned without the limit on the mutations),
//! see [issue #221](https://github.com/algbio/tsalign/issues/221).
//! The chaining minimises the chain cost, so the reported cost is not optimal among the alignments that
//! the anchors support, and depends on how ties between chains of equal chain cost are broken.
//! Moreover, the chain cost of an alignment may be larger than its cost, so a chain that represents the optimal
//! alignment without template switches may be more expensive than another chain.
//! Hence, for inexact anchors, the following checks are not applied:
//! the upper bound by the optimum without template switches, the exact cost of alignments without template switches
//! and their independence of `k` and `max_anchor_mutations`, monotonicity, and the symmetry of the reported cost
//! (the chain cost is still checked for symmetry).
//!
//! tschainalign checks the length of template switches only based on the anchors,
//! so the jump aligners may make a template switch shorter or longer than the config allows.
//! This is accepted for now: for such alignments, the cost is recomputed without geometry costs,
//! and the lower bound from the classic aligner with template switches is not checked.

#![allow(dead_code)]

use std::{
    collections::HashMap,
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        Arc, Mutex, OnceLock,
        atomic::{AtomicUsize, Ordering},
    },
    time::{Duration, Instant},
};

use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U32Cost;
use lib_tsalign::{
    alignment_result::{AlignmentType, alignment::Alignment},
    config::TemplateSwitchConfig,
};

/// The configuration directory used for all tests, relative to the repository root.
pub const CONFIG: &str = "test_files/config/chainalign";

/// A parameter combination: `(k, max_anchor_mutations)`.
pub type Params = (u32, u8);

/// The default parameter grid.
/// Runs with two anchor mutations are slow, so they are only used for short sequences.
pub fn default_grid(max_length: usize) -> Vec<Params> {
    let mut grid = vec![(4, 0), (6, 0), (9, 0), (12, 0), (6, 1), (9, 1), (12, 1)];
    if max_length <= 60 {
        grid.extend([(10, 2), (13, 2)]);
    }
    grid
}

/// A reduced grid for sequences longer than 100 characters.
/// Runs with anchor mutations take more than a minute on some of these, so only one such combination is used.
pub fn long_grid() -> Vec<Params> {
    vec![(6, 0), (9, 0), (12, 0), (12, 1)]
}

/// A test case: one pair of sequences aligned with each parameter combination of a grid.
#[derive(Debug, Clone)]
pub struct Case {
    pub name: String,
    pub reference: Vec<u8>,
    pub query: Vec<u8>,
    pub grid: Vec<Params>,
    /// Check that swapping reference and query does not change the cost (for the first parameter combination of each max_anchor_mutations).
    pub check_symmetry: bool,
    /// Check that the classic A* aligner with template switches does not find an alignment more expensive than tschainalign.
    pub check_ts_lower_bound: bool,
}

impl Case {
    pub fn new(
        name: impl Into<String>,
        reference: impl AsRef<[u8]>,
        query: impl AsRef<[u8]>,
    ) -> Self {
        let reference = reference.as_ref().to_vec();
        let query = query.as_ref().to_vec();
        let max_length = reference.len().max(query.len());
        let grid = if max_length > 100 {
            long_grid()
        } else {
            default_grid(max_length)
        };
        Self {
            name: name.into(),
            reference,
            query,
            grid,
            check_symmetry: true,
            check_ts_lower_bound: true,
        }
    }

    pub fn with_grid(mut self, grid: impl IntoIterator<Item = Params>) -> Self {
        self.grid = grid.into_iter().collect();
        self
    }

    pub fn without_symmetry(mut self) -> Self {
        self.check_symmetry = false;
        self
    }

    pub fn without_ts_lower_bound(mut self) -> Self {
        self.check_ts_lower_bound = false;
        self
    }

    pub fn max_length(&self) -> usize {
        self.reference.len().max(self.query.len())
    }

    pub fn fasta(&self) -> String {
        format!(
            ">reference\n{}\n>query\n{}\n",
            String::from_utf8_lossy(&self.reference),
            String::from_utf8_lossy(&self.query)
        )
    }
}

/// Read a pair fasta file from the repository, skipping `skip_characters`.
pub fn read_pair_fasta(path: &str, skip_characters: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let content = fs::read_to_string(repo_root().join(path))
        .unwrap_or_else(|error| panic!("Error reading {path}: {error}"));
    let mut records: Vec<Vec<u8>> = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('>') {
            records.push(Vec::new());
        } else if let Some(record) = records.last_mut() {
            record.extend(line.bytes().filter(|c| !skip_characters.contains(c)));
        }
    }
    assert_eq!(records.len(), 2, "{path} is not a pair fasta file");
    let query = records.pop().unwrap();
    let reference = records.pop().unwrap();
    (reference, query)
}

/// Run all checks for a case and panic with a report of all violations.
pub fn run_case(case: &Case) {
    if let Err(report) = check_case(case) {
        panic!("{report}");
    }
}

/// Run many cases in parallel and panic with a report of all failing cases.
pub fn run_cases(cases: &[Case]) {
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .clamp(1, 4);
    let next = AtomicUsize::new(0);
    let reports = Mutex::new(Vec::new());

    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    let Some(case) = cases.get(index) else {
                        break;
                    };
                    if let Err(report) = check_case(case) {
                        reports.lock().unwrap().push((index, report));
                    }
                }
            });
        }
    });

    let mut reports = reports.into_inner().unwrap();
    if !reports.is_empty() {
        reports.sort_by_key(|(index, _)| *index);
        let failed = reports.len();
        let report = reports
            .into_iter()
            .map(|(_, report)| report)
            .collect::<Vec<_>>()
            .join("\n");
        panic!("{failed} of {} cases failed:\n\n{report}", cases.len());
    }
}

/// Run all checks for a case and return a report of all violations.
pub fn check_case(case: &Case) -> Result<(), String> {
    let mut violations = Vec::new();
    let reference = case.reference.as_slice();
    let query = case.query.as_slice();

    let no_ts = run_classic(reference, query, false);
    let no_ts_cost = match &no_ts {
        RunOutcome::Success(result) => Some(result.cost),
        other => {
            violations.push(format!("classic aligner without TS: {other}"));
            None
        }
    };
    let ts_cost = if case.check_ts_lower_bound {
        match run_classic(reference, query, true) {
            RunOutcome::Success(result) => Some(result.cost),
            other => {
                violations.push(format!("classic aligner with TS: {other}"));
                None
            }
        }
    } else {
        None
    };

    let mut costs: Vec<(Params, u32, Option<u32>, u32)> = Vec::new();
    for &(k, max_anchor_mutations) in &case.grid {
        let params = format!("k = {k}, max_anchor_mutations = {max_anchor_mutations}");
        let result = match run_chain(reference, query, k, max_anchor_mutations) {
            RunOutcome::Success(result) => result,
            other => {
                violations.push(format!("{params}: {other}"));
                continue;
            }
        };

        violations.extend(
            check_alignment(reference, query, &result)
                .into_iter()
                .map(|violation| format!("{params}: {violation}")),
        );

        let is_exact = max_anchor_mutations == 0;
        match result.chain_cost {
            None => violations.push(format!("{params}: no chain cost reported")),
            // The cost of the alignment as a whole is at most the sum of the costs of its parts (issue #221).
            Some(chain_cost) if result.cost > chain_cost => violations.push(format!(
                "{params}: cost {} is larger than the chain cost {chain_cost}",
                result.cost
            )),
            // Exact anchors consist of matches only, so the alignment costs exactly as much as its parts.
            Some(chain_cost) if is_exact && result.cost != chain_cost => {
                violations.push(format!(
                    "{params}: cost {} differs from the chain cost {chain_cost}, although the anchors are exact",
                    result.cost
                ))
            }
            Some(_) => {}
        }

        if let Some(no_ts_cost) = no_ts_cost {
            // For inexact anchors, the chain cost of an alignment may be larger than its cost (issue #221),
            // so the chaining may prefer a more expensive alignment over the optimal alignment without template switches.
            if is_exact && result.cost > no_ts_cost {
                violations.push(format!(
                    "{params}: cost {} is larger than the cost {no_ts_cost} of the optimal alignment without template switches",
                    result.cost
                ));
            }
            if result.template_switch_amount == 0
                && (result.cost < no_ts_cost || (is_exact && result.cost > no_ts_cost))
            {
                violations.push(format!(
                    "{params}: alignment without template switches has cost {}, but the optimal alignment without template switches has cost {no_ts_cost}",
                    result.cost
                ));
            }
        }
        // tschainalign may use template switches whose length is outside the range allowed by the config,
        // see `check_alignment`. These may be cheaper than what the classic aligner can find.
        if let Some(ts_cost) = ts_cost
            && result.cost < ts_cost
            && template_switch_lengths_in_range(&result.alignment)
        {
            violations.push(format!(
                "{params}: cost {} is smaller than the cost {ts_cost} found by the classic aligner with template switches",
                result.cost
            ));
        }

        // Pruning dominated anchors must not change the cost of the optimal chain.
        // The chaining with all anchors is the oracle for the chain cost of the default (pruned) run.
        // Exact anchors are never pruned.
        // The inexact anchor pruning mode must not fail, and must not decrease the chain cost.
        if !is_exact {
            let unpruned = run_chain_unpruned(reference, query, k, max_anchor_mutations);
            violations.extend(
                check_chain_cost_against_unpruned_outcome(&unpruned, &result)
                    .into_iter()
                    .chain(check_inexact_pruning(
                        reference,
                        query,
                        k,
                        max_anchor_mutations,
                        &unpruned,
                    ))
                    .map(|violation| format!("{params}: {violation}")),
            );
        }

        costs.push((
            (k, max_anchor_mutations),
            result.cost,
            result.chain_cost,
            result.template_switch_amount,
        ));
    }

    // The following checks reason about the optimality of the reported alignment among the alignments supported by the anchors.
    // For inexact anchors, the chaining does not optimise the reported cost, but the chain cost, which may be larger
    // than the cost of the reported alignment, and the chain cost of the same alignment may differ between different
    // k and max_anchor_mutations (issue #221). So these checks are only applied to exact anchors.
    let exact_costs: Vec<_> = costs
        .iter()
        .filter(|&&((_, m), ..)| m == 0)
        .map(|&(params, cost, _, ts)| (params, cost, ts))
        .collect();

    // Alignments without template switches do not depend on k.
    // (This also follows from comparing against the classic aligner, but holds even if the classic aligner fails.)
    for &((k1, m1), cost1, ts1) in &exact_costs {
        for &((k2, m2), cost2, ts2) in &exact_costs {
            if (k1, m1) < (k2, m2) && ts1 == 0 && ts2 == 0 && cost1 != cost2 {
                violations.push(format!(
                    "alignments without template switches have different costs: {cost1} for k = {k1}, max_anchor_mutations = {m1} \
                     and {cost2} for k = {k2}, max_anchor_mutations = {m2}"
                ));
            }
        }
    }

    // Monotonicity: if every template switch supported by the anchors of (k1, m1) is also supported by the anchors of (k2, m2),
    // then the cost for (k2, m2) is at most the cost for (k1, m1).
    for &((k1, m1), cost1, _) in &exact_costs {
        for &((k2, m2), cost2, _) in &exact_costs {
            if (k1, m1) != (k2, m2) && anchors_dominate((k2, m2), (k1, m1)) && cost2 > cost1 {
                violations.push(format!(
                    "cost {cost2} for k = {k2}, max_anchor_mutations = {m2} is larger than cost {cost1} for k = {k1}, max_anchor_mutations = {m1}, \
                     although every template switch supported by anchors in the latter is supported by anchors in the former"
                ));
            }
        }
    }

    // Symmetry: swapping reference and query does not change the chain cost.
    // For inexact anchors, chains of equal chain cost may represent alignments of different costs (issue #221),
    // so the reported cost may depend on the order of the sequences.
    // For exact anchors, the reported cost equals the chain cost, which is checked above.
    if case.check_symmetry {
        let mut seen_mutations = Vec::new();
        for &((k, max_anchor_mutations), _, chain_cost, _) in &costs {
            if seen_mutations.contains(&max_anchor_mutations) {
                continue;
            }
            seen_mutations.push(max_anchor_mutations);

            match run_chain(query, reference, k, max_anchor_mutations) {
                RunOutcome::Success(result) => {
                    if result.chain_cost != chain_cost {
                        violations.push(format!(
                            "k = {k}, max_anchor_mutations = {max_anchor_mutations}: chain cost {} with swapped reference and query differs from chain cost {}",
                            display_optional(result.chain_cost),
                            display_optional(chain_cost),
                        ));
                    }
                }
                other => violations.push(format!(
                    "k = {k}, max_anchor_mutations = {max_anchor_mutations} with swapped reference and query: {other}"
                )),
            }
        }
    }

    if violations.is_empty() {
        Ok(())
    } else {
        Err(failure_report(case, &violations))
    }
}

/// Check the chain cost of a tschainalign result against the chain cost of the same alignment
/// without pruning dominated anchors.
///
/// Removing anchors can only remove chains, so the chain cost with pruning is never smaller.
/// Pruning is sound if the chain cost is not larger either.
/// The reported cost is not compared, since it depends on how ties between chains of equal chain cost are broken
/// (see [issue #221](https://github.com/algbio/tsalign/issues/221)).
pub fn check_chain_cost_against_unpruned(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
    result: &ChainResult,
) -> Option<String> {
    check_chain_cost_against_unpruned_outcome(
        &run_chain_unpruned(reference, query, k, max_anchor_mutations),
        result,
    )
}

fn check_chain_cost_against_unpruned_outcome(
    unpruned: &RunOutcome,
    result: &ChainResult,
) -> Option<String> {
    match unpruned {
        RunOutcome::Success(unpruned) => (result.chain_cost != unpruned.chain_cost).then(|| {
            format!(
                "chain cost {} differs from the chain cost {} without pruning dominated anchors \
                 (costs {} and {}, alignments {} and {})",
                display_optional(result.chain_cost),
                display_optional(unpruned.chain_cost),
                result.cost,
                unpruned.cost,
                result.alignment,
                unpruned.alignment,
            )
        }),
        // Without the oracle, the chain cost is not verified, which must not pass silently.
        other => Some(format!("without pruning dominated anchors: {other}")),
    }
}

/// Check a run with the inexact anchor pruning mode (`--anchor-pruning inexact`).
///
/// The inexact pruning may increase the chain cost and lose template switches, so the costs are not compared for
/// equality. The run must succeed with a valid alignment, and its chain cost must be at least the chain cost without
/// pruning, since removing anchors only removes chains.
/// The reported cost is not compared (it may be lower than the chain cost, see
/// [issue #221](https://github.com/algbio/tsalign/issues/221)).
fn check_inexact_pruning(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
    unpruned: &RunOutcome,
) -> Vec<String> {
    let result = match run_chain_inexact_pruning(reference, query, k, max_anchor_mutations) {
        RunOutcome::Success(result) => result,
        other => return vec![format!("with inexact anchor pruning: {other}")],
    };
    let mut violations: Vec<_> = check_alignment(reference, query, &result)
        .into_iter()
        .map(|violation| format!("with inexact anchor pruning: {violation}"))
        .collect();
    if let RunOutcome::Success(unpruned) = unpruned
        && let (Some(chain_cost), Some(unpruned_chain_cost)) =
            (result.chain_cost, unpruned.chain_cost)
        && chain_cost < unpruned_chain_cost
    {
        violations.push(format!(
            "chain cost {chain_cost} with inexact anchor pruning is smaller than the chain cost {unpruned_chain_cost} without pruning"
        ));
    }
    if result.chain_cost.is_none() {
        violations.push("with inexact anchor pruning: no chain cost reported".to_string());
    }
    violations
}

fn display_optional(value: Option<u32>) -> String {
    value.map_or_else(|| "<missing>".to_string(), |value| value.to_string())
}

/// The minimum length of the inner alignment of a template switch in [`CONFIG`].
const MIN_TS_LENGTH: u32 = 10;

/// Returns true if every template switch that is supported by the anchors for `weaker`
/// is also supported by the anchors for `stronger`.
///
/// A template switch is entered at a secondary anchor and left after a chain of secondary anchors
/// whose span must be at least the minimum template switch length.
/// Anchors of a chain must not overlap, except for exact anchors shifted by one.
/// Inexact anchors allow up to max_anchor_mutations edit operations.
///
/// * Same k, more mutations: the anchor set only grows.
/// * Exact anchors, smaller k: each exact k-mer anchor is covered by a chain of shifted exact k'-mer anchors with the same span.
/// * Smaller k, more mutations, k' >= minimum template switch length: a single k'-mer anchor
///   inside the first anchor of the weaker chain suffices to support the template switch.
///
/// In general, a smaller k may increase the cost if anchors are inexact:
/// e.g. an inner alignment of length 10 with one mismatch in the middle is supported by k = 10, m = 1, but not by k = 6, m = 1.
///
/// Note that supporting a template switch does not imply that the chaining finds the alignment for inexact anchors,
/// since the chain cost of an alignment may be larger than its cost, see [issue #221](https://github.com/algbio/tsalign/issues/221).
/// Hence, the monotonicity check in [`check_case`] uses this only for exact anchors.
fn anchors_dominate(stronger: Params, weaker: Params) -> bool {
    let ((k_strong, m_strong), (k_weak, m_weak)) = (stronger, weaker);
    if k_strong > k_weak || m_strong < m_weak {
        return false;
    }
    k_strong == k_weak || (m_strong == 0 && m_weak == 0) || k_strong >= MIN_TS_LENGTH
}

fn failure_report(case: &Case, violations: &[String]) -> String {
    let failure_directory = work_root().join("failures");
    let _ = fs::create_dir_all(&failure_directory);
    let file_name: String = case
        .name
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    let fasta_path = failure_directory.join(format!("{file_name}.fa"));
    let _ = fs::write(&fasta_path, case.fasta());

    let mut report = String::new();
    writeln!(report, "=== Case {} ===", case.name).unwrap();
    write!(report, "{}", case.fasta()).unwrap();
    writeln!(report, "Violations:").unwrap();
    for violation in violations {
        writeln!(report, "  - {violation}").unwrap();
    }
    writeln!(
        report,
        "Reproduce from the repository root with:\n  cargo run --profile e2e -- align --alignment-method a-star-chain-ts -c {CONFIG} --alphabet dna -k <K> --max-anchor-mutations <M> -p {}",
        fasta_path.display()
    )
    .unwrap();
    report
}

/// Check the reported alignment by recomputing its cost with lib_tsalign.
///
/// tschainalign checks the length of template switches only based on the anchors,
/// and ignores the parts of the inner alignment that are aligned by the jump aligners.
/// Hence, it may report template switches whose length is outside the range allowed by the config.
/// This is accepted for now, and for such alignments, the cost is recomputed without geometry costs.
fn check_alignment(reference: &[u8], query: &[u8], result: &ChainResult) -> Vec<String> {
    let mut violations = Vec::new();
    let reference_genome = VectorGenome::<DnaAlphabet>::from_slice_u8(reference).unwrap();
    let query_genome = VectorGenome::<DnaAlphabet>::from_slice_u8(query).unwrap();
    let ignore_geometry_cost = !template_switch_lengths_in_range(&result.alignment);

    let computed = std::panic::catch_unwind(|| {
        result.alignment.compute_cost_and_lengths(
            reference_genome.as_genome_subsequence(),
            query_genome.as_genome_subsequence(),
            0,
            0,
            tsa_config(),
            ignore_geometry_cost,
        )
    });
    match computed {
        Ok((cost, reference_length, query_length)) => {
            if cost != U32Cost::from(result.cost) {
                violations.push(format!(
                    "reported cost {} differs from the recomputed cost {cost} of the reported alignment {} \
                     (an infinite recomputed cost may also mean that tschainalign allows a template switch that the classic model forbids)",
                    result.cost, result.alignment
                ));
            }
            if (reference_length, query_length) != (reference.len(), query.len()) {
                violations.push(format!(
                    "reported alignment {} consumes {reference_length} reference and {query_length} query characters, but the sequences have lengths {} and {}",
                    result.alignment,
                    reference.len(),
                    query.len()
                ));
            }
        }
        Err(_) => violations.push(format!(
            "recomputing the cost of the reported alignment {} panicked",
            result.alignment
        )),
    }
    violations
}

/// Returns true if the inner alignments of all template switches have a length in the descendant
/// that is allowed by the config.
///
/// tschainalign checks the length only based on the anchors, and ignores the parts of the inner alignment
/// that are aligned by the jump aligners, so it may violate this.
pub fn template_switch_lengths_in_range(alignment: &Alignment) -> bool {
    let length_range = tsa_config().length_costs.zero_range().unwrap();
    // The length of the inner alignment of the current template switch in the descendant.
    let mut length = 0;
    for alignment_type in alignment.iter_flat_cloned() {
        match alignment_type {
            AlignmentType::TemplateSwitchEntrance { .. } => length = 0,
            AlignmentType::SecondaryMatch
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryInsertion => length += 1,
            AlignmentType::TemplateSwitchExit { .. } if !length_range.contains(&length) => {
                return false;
            }
            _ => {}
        }
    }
    true
}

fn tsa_config() -> &'static TemplateSwitchConfig<DnaAlphabet, U32Cost> {
    static CONFIG_CELL: OnceLock<TemplateSwitchConfig<DnaAlphabet, U32Cost>> = OnceLock::new();
    CONFIG_CELL.get_or_init(|| {
        let path = repo_root().join(CONFIG).join("config.tsa");
        TemplateSwitchConfig::read_plain(fs::File::open(&path).unwrap()).unwrap()
    })
}

/// The result of a successful alignment run.
#[derive(Debug, Clone)]
pub struct ChainResult {
    pub cost: u32,
    /// The cost of the chain that the alignment was built from (tschainalign only).
    /// It may be larger than `cost` for inexact anchors, see [issue #221](https://github.com/algbio/tsalign/issues/221).
    pub chain_cost: Option<u32>,
    pub template_switch_amount: u32,
    pub alignment: Alignment,
}

#[derive(Debug, Clone)]
pub enum RunOutcome {
    Success(ChainResult),
    Failure { message: String },
    Timeout { seconds: u64 },
}

impl std::fmt::Display for RunOutcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunOutcome::Success(result) => write!(
                f,
                "success with cost {} and {} template switches",
                result.cost, result.template_switch_amount
            ),
            RunOutcome::Failure { message } => write!(f, "failed: {message}"),
            RunOutcome::Timeout { seconds } => write!(f, "timed out after {seconds}s"),
        }
    }
}

/// Align with tschainalign.
pub fn run_chain(reference: &[u8], query: &[u8], k: u32, max_anchor_mutations: u8) -> RunOutcome {
    run_chain_with_pruning(reference, query, k, max_anchor_mutations, "exact")
}

/// Align with tschainalign without pruning dominated anchors (`--anchor-pruning none`).
///
/// This is the oracle for the chain cost of [`run_chain`].
pub fn run_chain_unpruned(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
) -> RunOutcome {
    run_chain_with_pruning(reference, query, k, max_anchor_mutations, "none")
}

/// Align with tschainalign with the inexact anchor pruning mode and all its rules (`--anchor-pruning inexact`).
pub fn run_chain_inexact_pruning(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
) -> RunOutcome {
    run_chain_with_pruning(reference, query, k, max_anchor_mutations, "inexact")
}

/// Align with tschainalign with the given `--anchor-pruning` mode.
fn run_chain_with_pruning(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
    anchor_pruning: &str,
) -> RunOutcome {
    let mut args: Vec<String> = [
        "align",
        "--alignment-method",
        "a-star-chain-ts",
        "--alphabet",
        "dna",
    ]
    .map(String::from)
    .into();
    args.extend(["-k".to_string(), k.to_string()]);
    args.extend([
        "--max-anchor-mutations".to_string(),
        max_anchor_mutations.to_string(),
    ]);

    let cache = match preprocessed_cache(k, max_anchor_mutations) {
        Ok(cache) => cache,
        Err(message) => return RunOutcome::Failure { message },
    };
    args.extend([
        "--cache-directory".to_string(),
        cache.display().to_string(),
        "--force-no-preprocessing".to_string(),
    ]);
    args.extend(["--anchor-pruning".to_string(), anchor_pruning.to_string()]);

    // The chaining with all anchors is slower, and its result is only used as an oracle.
    let timeout = if anchor_pruning == "none" {
        timeout() * UNPRUNED_TIMEOUT_FACTOR
    } else {
        timeout()
    };
    run_align(reference, query, args, timeout)
}

/// The factor by which the timeout is increased for runs without pruning dominated anchors.
///
/// Pruning dominated anchors makes the chaining faster, so the oracle run without pruning may exceed the normal
/// timeout although the run with pruning does not.
/// If the oracle run exceeds the increased timeout as well, the chain cost cannot be verified, which is a violation.
const UNPRUNED_TIMEOUT_FACTOR: u32 = 3;

/// Align with the classic A* aligner, with or without template switches.
pub fn run_classic(reference: &[u8], query: &[u8], template_switches: bool) -> RunOutcome {
    /// Maps (reference, query, template_switches) to the outcome.
    type OracleCache = Mutex<HashMap<(Vec<u8>, Vec<u8>, bool), RunOutcome>>;
    static ORACLE_CACHE: OnceLock<OracleCache> = OnceLock::new();
    let key = (reference.to_vec(), query.to_vec(), template_switches);
    let cache = ORACLE_CACHE.get_or_init(Default::default);
    if let Some(outcome) = cache.lock().unwrap().get(&key) {
        return outcome.clone();
    }

    let mut args: Vec<String> = [
        "align",
        "--alignment-method",
        "a-star-template-switch",
        "--alphabet",
        "dna",
    ]
    .map(String::from)
    .into();
    if !template_switches {
        args.push("--no-ts".to_string());
    }
    let outcome = run_align(reference, query, args, timeout());
    cache.lock().unwrap().insert(key, outcome.clone());
    outcome
}

fn run_align(
    reference: &[u8],
    query: &[u8],
    mut args: Vec<String>,
    timeout: Duration,
) -> RunOutcome {
    let input = unique_path("pair.fa");
    let output = input.with_extension("toml");
    fs::write(
        &input,
        format!(
            ">reference\n{}\n>query\n{}\n",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(query)
        ),
    )
    .unwrap();
    args.extend([
        "-c".to_string(),
        CONFIG.to_string(),
        "-l".to_string(),
        "warn".to_string(),
        "-p".to_string(),
        input.display().to_string(),
        "-o".to_string(),
        output.display().to_string(),
    ]);

    let outcome = match run_tsalign(&args, timeout) {
        Err(outcome) => outcome,
        Ok(()) => match parse_output(&output) {
            Ok(result) => RunOutcome::Success(result),
            Err(message) => RunOutcome::Failure { message },
        },
    };

    let _ = fs::remove_file(&input);
    let _ = fs::remove_file(&output);
    outcome
}

fn parse_output(path: &Path) -> Result<ChainResult, String> {
    let content = fs::read_to_string(path).map_err(|error| format!("reading {path:?}: {error}"))?;
    let value: toml::Table =
        toml::from_str(&content).map_err(|error| format!("parsing output: {error}"))?;

    let result_type = value
        .get("type")
        .and_then(|v| v.as_str())
        .unwrap_or("<missing>");
    if result_type != "WithTarget" {
        return Err(format!("no alignment found (result type {result_type})"));
    }

    let cost = value
        .get("cost")
        .and_then(|v| v.as_float())
        .ok_or("missing cost")?;
    if !cost.is_finite() || cost < 0.0 || cost.fract() != 0.0 {
        return Err(format!("invalid cost {cost}"));
    }
    let chain_cost = match value.get("chain_cost") {
        None => None,
        Some(chain_cost) => {
            let chain_cost = chain_cost.as_float().ok_or("chain cost is not a number")?;
            if !chain_cost.is_finite() || chain_cost < 0.0 || chain_cost.fract() != 0.0 {
                return Err(format!("invalid chain cost {chain_cost}"));
            }
            Some(chain_cost as u32)
        }
    };
    let template_switch_amount = value
        .get("template_switch_amount")
        .and_then(|v| v.as_float())
        .ok_or("missing template_switch_amount")?;
    let alignment: Alignment = toml::Value::Table(value)
        .try_into()
        .map_err(|error| format!("parsing alignment: {error}"))?;

    Ok(ChainResult {
        cost: cost as u32,
        chain_cost,
        template_switch_amount: template_switch_amount as u32,
        alignment,
    })
}

/// Run tsalign from the repository root and return its exit status and stderr, without interpreting them.
///
/// This is meant for checking how tsalign rejects invalid inputs.
pub fn run_tsalign_raw(args: &[String]) -> (std::process::ExitStatus, String) {
    let output = run_tsalign_output(args);
    (
        output.status,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// Run tsalign from the repository root and return its exit status, stdout and stderr, without interpreting them.
pub fn run_tsalign_output(args: &[String]) -> std::process::Output {
    let _permit = ProcessPermit::acquire();
    Command::new(tsalign_binary())
        .args(args)
        .current_dir(repo_root())
        .env("RUST_BACKTRACE", "0")
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// A fresh path in the work directory of this test process, ending with `suffix`.
pub fn scratch_path(suffix: &str) -> PathBuf {
    unique_path(suffix)
}

/// Run tsalign with a timeout. Returns `Err` with the failure outcome if it does not succeed.
fn run_tsalign(args: &[String], timeout: Duration) -> Result<(), RunOutcome> {
    let _permit = ProcessPermit::acquire();
    let stdout_path = unique_path("stdout.txt");
    let stderr_path = unique_path("stderr.txt");
    let mut child = Command::new(tsalign_binary())
        .args(args)
        .current_dir(repo_root())
        .env("RUST_BACKTRACE", "0")
        .stdin(Stdio::null())
        .stdout(fs::File::create(&stdout_path).unwrap())
        .stderr(fs::File::create(&stderr_path).unwrap())
        .spawn()
        .unwrap();

    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break Some(status);
        }
        if start.elapsed() > timeout {
            let _ = child.kill();
            let _ = child.wait();
            break None;
        }
        std::thread::sleep(Duration::from_millis(10));
    };

    let stderr = fs::read_to_string(&stderr_path).unwrap_or_default();
    let _ = fs::remove_file(&stdout_path);
    let _ = fs::remove_file(&stderr_path);

    match status {
        None => Err(RunOutcome::Timeout {
            seconds: timeout.as_secs(),
        }),
        Some(status) if status.success() => Ok(()),
        Some(status) => Err(RunOutcome::Failure {
            message: format!("{status}; {}", summarise_stderr(&stderr)),
        }),
    }
}

/// Extract the panic message from the stderr of a tsalign run, or return its last lines.
fn summarise_stderr(stderr: &str) -> String {
    let lines: Vec<&str> = stderr.lines().collect();
    if let Some(index) = lines.iter().position(|line| line.contains("panicked at")) {
        let message: Vec<&str> = lines[index..]
            .iter()
            .take_while(|line| !line.starts_with("note:") && !line.starts_with("stack backtrace"))
            .copied()
            .collect();
        message.join(" | ")
    } else {
        lines[lines.len().saturating_sub(5)..].join(" | ")
    }
}

/// Limits the number of concurrently running tsalign processes to the number of CPUs,
/// such that parallel tests do not cause timeouts by oversubscription.
struct ProcessPermit;

static RUNNING_PROCESSES: (Mutex<usize>, std::sync::Condvar) =
    (Mutex::new(0), std::sync::Condvar::new());

impl ProcessPermit {
    fn acquire() -> Self {
        let limit = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let (count, condvar) = &RUNNING_PROCESSES;
        let mut count = condvar
            .wait_while(count.lock().unwrap(), |count| *count >= limit)
            .unwrap();
        *count += 1;
        Self
    }
}

impl Drop for ProcessPermit {
    fn drop(&mut self) {
        let (count, condvar) = &RUNNING_PROCESSES;
        *count.lock().unwrap() -= 1;
        condvar.notify_one();
    }
}

fn timeout() -> Duration {
    let seconds = std::env::var("TSALIGN_E2E_TIMEOUT_SECS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(60);
    Duration::from_secs(seconds)
}

/// Preprocess the chaining lower bounds for `(k, max_anchor_mutations)` once per test process.
/// Returns the cache directory, or the reason why preprocessing failed.
pub fn preprocessed_cache(k: u32, max_anchor_mutations: u8) -> Result<PathBuf, String> {
    type CacheCell = Arc<OnceLock<Result<PathBuf, String>>>;
    static CACHES: OnceLock<Mutex<HashMap<Params, CacheCell>>> = OnceLock::new();
    let cell = CACHES
        .get_or_init(Default::default)
        .lock()
        .unwrap()
        .entry((k, max_anchor_mutations))
        .or_default()
        .clone();

    cell.get_or_init(|| {
        let directory = work_directory().join(format!("cache-{k}-{max_anchor_mutations}"));
        fs::create_dir_all(&directory).unwrap();
        let args: Vec<String> = vec![
            "preprocess".into(),
            "-l".into(),
            "warn".into(),
            "--alphabet".into(),
            "dna".into(),
            "-c".into(),
            CONFIG.into(),
            "--cache-directory".into(),
            directory.display().to_string(),
            "-k".into(),
            k.to_string(),
            "--max-anchor-mutations".into(),
            max_anchor_mutations.to_string(),
            "256".into(),
        ];
        run_tsalign(&args, timeout())
            .map(|()| directory)
            .map_err(|outcome| format!("preprocessing {outcome}"))
    })
    .clone()
}

pub fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

/// The cargo target directory, derived from the path of the running test executable (`<target>/<profile>/deps/<test>`).
fn target_directory() -> PathBuf {
    let executable = std::env::current_exe().unwrap();
    executable
        .ancestors()
        .find(|path| path.file_name().is_some_and(|name| name == "deps"))
        .and_then(|deps| deps.parent())
        .and_then(|profile| profile.parent())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| repo_root().join("target"))
}

fn work_root() -> PathBuf {
    target_directory().join("tsalign-e2e")
}

/// A fresh directory for this test process.
/// Lower bounds are always preprocessed anew, because cached lower bounds from older code versions may be wrong.
fn work_directory() -> &'static Path {
    static DIRECTORY: OnceLock<PathBuf> = OnceLock::new();
    DIRECTORY.get_or_init(|| {
        let directory = work_root().join(format!("run-{}", std::process::id()));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        directory
    })
}

fn unique_path(suffix: &str) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    work_directory().join(format!("{id}-{suffix}"))
}

/// Build the tsalign binary with the `e2e` profile once per test process.
fn tsalign_binary() -> &'static Path {
    static BINARY: OnceLock<PathBuf> = OnceLock::new();
    BINARY.get_or_init(|| {
        let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
        let status = Command::new(cargo)
            .args([
                "build",
                "--profile",
                "e2e",
                "-p",
                "tsalign",
                "--bin",
                "tsalign",
            ])
            .current_dir(repo_root())
            .status()
            .unwrap();
        assert!(status.success(), "Building tsalign failed");

        let binary = target_directory()
            .join("e2e")
            .join(format!("tsalign{}", std::env::consts::EXE_SUFFIX));
        assert!(binary.exists(), "tsalign binary not found at {binary:?}");
        binary
    })
}
