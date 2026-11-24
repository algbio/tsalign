use generic_a_star::cost::U32Cost;
use lib_tsalign::alignment_geometry::AlignmentRange;

use crate::{
    alignment::sequences::AlignmentSequences,
    anchors::{
        Anchors,
        inexact_pruning::{
            AnchorType, InexactPruningParameters, InexactPruningRule, InexactPruningRules,
            PruningContext, RawAnchor, RuleCounts, RuleKind, prune,
        },
    },
    chain_align::performance_parameters::{AlignmentPerformanceParameters, AnchorPruning},
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    dna_rc_fn,
};

fn costs() -> GapAffineCosts<U32Cost> {
    GapAffineCosts::new(2u32, 3, 1).into_costs()
}

/// The costs of `test_files/config/chainalign`.
fn chainalign_config_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        costs(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits {
            inter_jump_12: -300..301,
            intra_jump_12: -300..301,
            jump_34: isize::MIN..isize::MAX,
            length_23: 10..200,
            ancestor_gap: isize::MIN..isize::MAX,
        },
    )
}

fn anchor(a: std::ops::Range<usize>, b: std::ops::Range<usize>, cost: u32) -> RawAnchor<U32Cost> {
    RawAnchor {
        a_start: a.start,
        a_end: a.end,
        b_start: b.start,
        b_end: b.end,
        cost: cost.into(),
    }
}

/// The input of a pruning run on raw anchors.
struct Setup<'a> {
    seq_a: &'a [u8],
    seq_b: &'a [u8],
    k: usize,
    max_mutations: usize,
    costs: GapAffineCosts<U32Cost>,
    max_match_run_lower_bound: bool,
    rules: &'a [RuleKind],
    exact_overlap_diagonal_distance: usize,
}

impl<'a> Setup<'a> {
    fn new(seq_a: &'a [u8], seq_b: &'a [u8], k: usize, rules: &'a [RuleKind]) -> Self {
        Self {
            seq_a,
            seq_b,
            k,
            max_mutations: 1,
            costs: costs(),
            max_match_run_lower_bound: true,
            rules,
            exact_overlap_diagonal_distance: 1,
        }
    }

    fn pruned(&self, anchors: &[RawAnchor<U32Cost>]) -> (Vec<RawAnchor<U32Cost>>, RuleCounts) {
        let mut kept = anchors.to_vec();
        let rules: InexactPruningRules = self
            .rules
            .iter()
            .map(|kind| InexactPruningRule::new(AnchorType::Primary, *kind))
            .collect();
        let context = PruningContext {
            seq_a: self.seq_a,
            seq_b: self.seq_b,
            k: self.k,
            max_mutations: self.max_mutations,
            costs: &self.costs,
            max_match_run_lower_bound: self.max_match_run_lower_bound,
            rules: rules.rule_set(AnchorType::Primary),
            exact_overlap_diagonal_distance: self.exact_overlap_diagonal_distance,
        };
        let counts = prune(&mut kept, |anchor| *anchor, &context);
        assert_eq!(counts.total_removed(), anchors.len() - kept.len());
        (kept, counts)
    }

    fn kept(&self, anchors: &[RawAnchor<U32Cost>]) -> Vec<RawAnchor<U32Cost>> {
        self.pruned(anchors).0
    }
}

#[test]
fn prefix_domination_removes_anchor_with_trailing_column() {
    let seq_a = b"ACGTA";
    let seq_b = b"ACGT";
    let exact = anchor(0..4, 0..4, 0);
    let long = anchor(0..5, 0..4, 3);
    let setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::PrefixDomination]);
    assert_eq!(setup.kept(&[exact, long]), vec![exact]);
    // The mirrored rule does not apply: (1..5, 0..4) is not an anchor.
    let setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::SuffixDomination]);
    assert_eq!(setup.kept(&[exact, long]), vec![exact, long]);
}

/// Dominators may be dominated themselves, and all dominated anchors are removed at once.
#[test]
fn domination_cascades() {
    // With two mutations and k = 4: ACGT vs TACGTA.
    let seq_a = b"ACGT";
    let seq_b = b"TACGTA";
    let exact = anchor(0..4, 1..5, 0);
    let leading_gap = anchor(0..4, 0..5, 3);
    let trailing_gap = anchor(0..4, 1..6, 3);
    let both_gaps = anchor(0..4, 0..6, 6);
    let anchors = [both_gaps, leading_gap, exact, trailing_gap];
    let mut setup = Setup::new(
        seq_a,
        seq_b,
        4,
        &[RuleKind::PrefixDomination, RuleKind::SuffixDomination],
    );
    setup.max_mutations = 2;
    // `both_gaps` is dominated by `leading_gap` (prefix) and `trailing_gap` (suffix), which are both dominated by `exact`.
    let (kept, counts) = setup.pruned(&anchors);
    assert_eq!(kept, vec![exact]);
    assert_eq!(counts.removed(RuleKind::PrefixDomination), 2);
    assert_eq!(counts.removed(RuleKind::SuffixDomination), 1);
    assert_eq!(counts.matched(RuleKind::SuffixDomination), 2);

    // Only prefix domination: `leading_gap` has no dominator with the same start.
    setup.rules = &[RuleKind::PrefixDomination];
    assert_eq!(setup.kept(&anchors), vec![leading_gap, exact]);
}

#[test]
fn suffix_domination_removes_anchor_with_leading_column() {
    // A gap at the start: -ACGT vs AACGT.
    let seq_a = b"ACGT";
    let seq_b = b"AACGT";
    let exact = anchor(0..4, 1..5, 0);
    let long = anchor(0..4, 0..5, 3);
    let setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::SuffixDomination]);
    assert_eq!(setup.kept(&[exact, long]), vec![exact]);
    let setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::PrefixDomination]);
    assert_eq!(setup.kept(&[exact, long]), vec![exact, long]);
    // The dominator is too expensive.
    let setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::SuffixDomination]);
    let expensive = anchor(0..4, 1..5, 1);
    assert_eq!(setup.kept(&[expensive, long]), vec![expensive, long]);
}

/// The domination rules keep the cases excluded by the proof of the exact pruning, unless `drop-guards` is enabled.
#[test]
fn drop_guards_extends_the_domination_rules() {
    let seq_a = b"ACGTA";
    let seq_b = b"ACGT";
    let exact = anchor(0..4, 0..4, 0);
    let long = anchor(0..5, 0..4, 3);
    let domination = [RuleKind::PrefixDomination];
    let unguarded = [RuleKind::PrefixDomination, RuleKind::DropGuards];

    // gap_extend > gap_open.
    let mut setup = Setup::new(seq_a, seq_b, 4, &domination);
    setup.costs = GapAffineCosts::new(2u32, 3, 4).into_costs();
    assert_eq!(setup.kept(&[exact, long]), vec![exact, long]);
    setup.rules = &unguarded;
    let (kept, counts) = setup.pruned(&[exact, long]);
    assert_eq!(kept, vec![exact]);
    assert_eq!(counts.removed(RuleKind::DropGuards), 1);
    assert_eq!(counts.removed(RuleKind::PrefixDomination), 0);

    // Another lower bound than max-match-run.
    let mut setup = Setup::new(seq_a, seq_b, 4, &domination);
    setup.max_match_run_lower_bound = false;
    assert_eq!(setup.kept(&[exact, long]), vec![exact, long]);
    setup.rules = &unguarded;
    assert_eq!(setup.kept(&[exact, long]), vec![exact]);

    // Zero-cost anchors (free gaps).
    let free_long = anchor(0..5, 0..4, 0);
    let mut setup = Setup::new(seq_a, seq_b, 4, &domination);
    setup.costs = GapAffineCosts::new(2u32, 0, 0).into_costs();
    assert_eq!(setup.kept(&[exact, free_long]), vec![exact, free_long]);
    setup.rules = &unguarded;
    assert_eq!(setup.kept(&[exact, free_long]), vec![exact]);

    // Alone, it removes nothing.
    let mut setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::DropGuards]);
    setup.max_match_run_lower_bound = false;
    assert_eq!(setup.kept(&[exact, long]), vec![exact, long]);
}

#[test]
fn edge_substitution_removes_substitution_at_either_end() {
    let rules = [RuleKind::EdgeSubstitution];
    let setup_for =
        |seq_a: &'static [u8], seq_b: &'static [u8]| Setup::new(seq_a, seq_b, 4, &rules);
    let substitution = anchor(0..4, 0..4, 2);

    // Last column.
    assert_eq!(setup_for(b"ACGT", b"ACGA").kept(&[substitution]), vec![]);
    // First column.
    assert_eq!(setup_for(b"ACGT", b"TCGT").kept(&[substitution]), vec![]);
    // Middle column.
    assert_eq!(
        setup_for(b"ACGT", b"AGGT").kept(&[substitution]),
        vec![substitution]
    );
    // The anchor is cheaper than its gap-free alignment, so it is not a substitution anchor.
    let cheaper = anchor(0..4, 0..4, 1);
    assert_eq!(setup_for(b"ACGT", b"ACGA").kept(&[cheaper]), vec![cheaper]);
    // Exact anchors and anchors with gaps are not substitution anchors.
    let exact = anchor(0..4, 0..4, 0);
    let gap = anchor(0..4, 0..3, 3);
    assert_eq!(
        setup_for(b"ACGT", b"ACGA").kept(&[exact, gap]),
        vec![exact, gap]
    );

    // With two mutations, one of the substitutions must be at an end.
    let two = anchor(0..5, 0..5, 4);
    let mut setup = setup_for(b"ACGTA", b"AGGTT");
    setup.k = 5;
    setup.max_mutations = 2;
    assert_eq!(setup.kept(&[two]), vec![]);
    let mut setup = setup_for(b"ACGTA", b"AGCTA");
    setup.k = 5;
    setup.max_mutations = 2;
    assert_eq!(setup.kept(&[two]), vec![two]);
}

#[test]
fn exact_overlap_removes_inexact_anchors_near_exact_anchors() {
    let seq_a = b"ACGTACGTACGT";
    let seq_b = b"ACGTACGTACGT";
    let exact = anchor(2..6, 2..6, 0);
    let same_diagonal = anchor(4..8, 4..8, 2);
    let neighbour_diagonal = anchor(3..7, 2..6, 2);
    let far_diagonal = anchor(4..8, 2..6, 2);
    // Diagonals 2 (start) to 3 (end).
    let gap_reaching_diagonal = anchor(4..9, 2..6, 3);
    let other_diagonal = anchor(6..10, 2..6, 2);
    let anchors = [
        exact,
        same_diagonal,
        neighbour_diagonal,
        far_diagonal,
        gap_reaching_diagonal,
        other_diagonal,
    ];

    let mut setup = Setup::new(seq_a, seq_b, 4, &[RuleKind::ExactOverlap]);
    assert_eq!(
        setup.kept(&anchors),
        vec![exact, far_diagonal, gap_reaching_diagonal, other_diagonal]
    );

    setup.exact_overlap_diagonal_distance = 0;
    assert_eq!(
        setup.kept(&anchors),
        vec![
            exact,
            neighbour_diagonal,
            far_diagonal,
            gap_reaching_diagonal,
            other_diagonal
        ]
    );

    setup.exact_overlap_diagonal_distance = 2;
    assert_eq!(setup.kept(&anchors), vec![exact, other_diagonal]);

    // Same diagonal, but no shared position.
    setup.exact_overlap_diagonal_distance = 1;
    let after = anchor(6..10, 6..10, 2);
    assert_eq!(setup.kept(&[exact, after]), vec![exact, after]);
    // Without exact anchors, nothing is removed.
    assert_eq!(
        setup.kept(&[same_diagonal, after]),
        vec![same_diagonal, after]
    );
}

/// An anchor that matches several rules is attributed to the first of them.
#[test]
fn removed_anchors_are_attributed_to_the_first_matching_rule() {
    let seq_a = b"ACGTA";
    let seq_b = b"ACGT";
    let exact = anchor(0..4, 0..4, 0);
    let long = anchor(0..5, 0..4, 3);
    let setup = Setup::new(
        seq_a,
        seq_b,
        4,
        &[RuleKind::ExactOverlap, RuleKind::PrefixDomination],
    );
    let (kept, counts) = setup.pruned(&[exact, long]);
    assert_eq!(kept, vec![exact]);
    assert_eq!(counts.matched(RuleKind::PrefixDomination), 1);
    assert_eq!(counts.matched(RuleKind::ExactOverlap), 1);
    assert_eq!(counts.removed(RuleKind::PrefixDomination), 1);
    assert_eq!(counts.removed(RuleKind::ExactOverlap), 0);
}

#[test]
fn rules_are_parsed_and_printed() {
    let all: InexactPruningRules = "all".parse().unwrap();
    assert_eq!(all, InexactPruningRules::all());
    assert_eq!(all, InexactPruningRules::default());
    assert_eq!(all.iter().count(), 10);
    assert!("none".parse::<InexactPruningRules>().unwrap().is_empty());
    assert_eq!(InexactPruningRules::none().to_string(), "none");

    let rules: InexactPruningRules = "primary-prefix-domination, secondary-exact-overlap"
        .parse()
        .unwrap();
    assert_eq!(
        rules.iter().collect::<Vec<_>>(),
        vec![
            InexactPruningRule::PrimaryPrefixDomination,
            InexactPruningRule::SecondaryExactOverlap
        ]
    );
    assert_eq!(
        rules.to_string(),
        "primary-prefix-domination,secondary-exact-overlap"
    );
    assert_eq!(rules.to_string().parse::<InexactPruningRules>(), Ok(rules));

    let primary: InexactPruningRules = "primary".parse().unwrap();
    assert!(primary.any_for(AnchorType::Primary));
    assert!(!primary.any_for(AnchorType::Secondary));
    assert_eq!(primary.iter().count(), 5);

    let without: InexactPruningRules = "all,-secondary,-primary-exact-overlap".parse().unwrap();
    assert_eq!(without.iter().count(), 4);
    assert!(!without.contains(InexactPruningRule::PrimaryExactOverlap));
    assert!(without.contains(InexactPruningRule::PrimaryEdgeSubstitution));

    assert!("primary-foo".parse::<InexactPruningRules>().is_err());

    // The names agree with the clap value names.
    use clap::ValueEnum;
    for rule in InexactPruningRule::ALL {
        assert_eq!(rule.to_possible_value().unwrap().get_name(), rule.name());
    }
}

fn random_sequences(seed: u64, amount: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut state = seed;
    let mut next = |limit: usize| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) % limit as u64) as usize
    };

    (0..amount)
        .map(|_| {
            let alphabet: &[u8] = if next(2) == 0 { b"AC" } else { b"ACGT" };
            let reference: Vec<u8> = (0..next(18) + 6)
                .map(|_| alphabet[next(alphabet.len())])
                .collect();
            let mut query: Vec<u8> = reference
                .iter()
                .flat_map(|&c| match next(10) {
                    0 => vec![alphabet[next(alphabet.len())]],
                    1 => vec![],
                    2 => vec![c, alphabet[next(alphabet.len())]],
                    _ => vec![c],
                })
                .collect();
            if query.is_empty() {
                query.push(b'A');
            }
            if next(2) == 0 && query.len() >= 10 {
                let length = 10 + next(query.len() - 9);
                let offset = next(query.len() - length + 1);
                let rc: Vec<u8> = query[offset..offset + length]
                    .iter()
                    .rev()
                    .copied()
                    .map(dna_rc_fn)
                    .collect();
                query[offset..offset + length].copy_from_slice(&rc);
            }
            (reference, query)
        })
        .collect()
}

/// The first line of the display of the anchors (primary anchors) and the other lines (secondary anchors).
fn anchor_lines(anchors: &Anchors<U32Cost>) -> (String, Vec<String>) {
    let display = format!("{anchors}");
    let mut lines = display.lines().map(str::to_string);
    let primary = lines.next().unwrap();
    (primary, lines.collect())
}

/// With one mutation, the inexact prefix domination with guards removes the same primary anchors as the exact pruning.
/// Rules for primary anchors leave secondary anchors untouched and vice versa.
#[test]
fn inexact_rules_on_generated_anchors() {
    let costs = chainalign_config_costs();
    let mut secondary_removed = 0;
    for (reference, query) in random_sequences(0x236, 20) {
        for k in [3u32, 5] {
            let sequences = AlignmentSequences::new_complete(reference.clone(), query.clone());
            let generate = || Anchors::new_inexact(&sequences, k, 1, &costs, &dna_rc_fn);
            let all = generate();
            let parameters = |rules: &str| InexactPruningParameters {
                rules: rules.parse().unwrap(),
                ..Default::default()
            };

            let mut exact = generate();
            exact.prune_dominated_primary_anchors(&sequences, k, 1, &costs);
            let mut inexact = generate();
            let (primary_counts, secondary_counts) = inexact.prune_inexact(
                &sequences,
                k,
                1,
                &costs,
                true,
                &parameters("primary-prefix-domination"),
                &dna_rc_fn,
            );
            assert_eq!(exact, inexact);
            assert_eq!(secondary_counts, RuleCounts::default());
            assert_eq!(
                primary_counts.total_removed(),
                all.primary_anchor_count() - inexact.primary_anchor_count()
            );

            let mut primary_only = generate();
            primary_only.prune_inexact(
                &sequences,
                k,
                1,
                &costs,
                true,
                &parameters("primary"),
                &dna_rc_fn,
            );
            assert_eq!(anchor_lines(&primary_only).1, anchor_lines(&all).1);

            let mut secondary_only = generate();
            let (primary_counts, secondary_counts) = secondary_only.prune_inexact(
                &sequences,
                k,
                1,
                &costs,
                true,
                &parameters("secondary"),
                &dna_rc_fn,
            );
            assert_eq!(anchor_lines(&secondary_only).0, anchor_lines(&all).0);
            assert_eq!(primary_counts, RuleCounts::default());
            assert_eq!(
                secondary_counts.total_removed(),
                all.secondary_anchor_count() - secondary_only.secondary_anchor_count()
            );
            secondary_removed += secondary_counts.total_removed();
        }
    }
    assert!(secondary_removed > 0);
}

fn lower_bounds(k: u32, max_anchor_mutations: u8, max_n: usize) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(
        max_n,
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        chainalign_config_costs(),
    )
}

/// Returns the chain cost and the number of template switches.
fn chain_cost(
    reference: &[u8],
    query: &[u8],
    lower_bounds: &ChainingLowerBounds<U32Cost>,
    anchor_pruning: AnchorPruning,
    rules: InexactPruningRules,
) -> (f64, f64) {
    let performance_parameters = AlignmentPerformanceParameters {
        max_exact_cost_function_cost: U32Cost::from(0u32),
        anchor_pruning,
        inexact_pruning: InexactPruningParameters {
            rules,
            ..Default::default()
        },
        ..Default::default()
    };
    let result = crate::align(
        reference.to_vec(),
        query.to_vec(),
        AlignmentRange::new_complete(reference.len(), query.len()),
        &performance_parameters,
        &dna_rc_fn,
        "reference",
        "query",
        lower_bounds,
    );
    let statistics = result.statistics();
    let chain_statistics = statistics.chain_statistics.as_ref().unwrap();
    assert_eq!(
        chain_statistics.primary_anchors + chain_statistics.secondary_anchors
            - chain_statistics.primary_anchors_after_pruning
            - chain_statistics.secondary_anchors_after_pruning,
        chain_statistics.pruned_anchors.values().sum::<u64>()
    );
    (
        statistics.chain_cost.unwrap().raw(),
        statistics.template_switch_amount.raw(),
    )
}

/// Removing anchors only removes chains, so the chain cost with inexact pruning is never smaller than without pruning.
#[test]
fn inexact_pruning_never_decreases_the_chain_cost() {
    // Chaining with small k is slow in debug builds, and two mutations even more so, so only a few cases are run.
    // The end-to-end tests check this on all their inputs.
    let cases = random_sequences(0xbeef, 8);
    for (k, max_anchor_mutations, amount) in [(5, 1, 8), (6, 2, 3)] {
        let lower_bounds = lower_bounds(k, max_anchor_mutations, 32);
        for (reference, query) in cases.iter().take(amount) {
            let unpruned = chain_cost(
                reference,
                query,
                &lower_bounds,
                AnchorPruning::None,
                InexactPruningRules::none(),
            )
            .0;
            for rules in ["all", "primary", "secondary", "none"] {
                let pruned = chain_cost(
                    reference,
                    query,
                    &lower_bounds,
                    AnchorPruning::Inexact,
                    rules.parse().unwrap(),
                )
                .0;
                assert!(
                    pruned >= unpruned,
                    "k = {k}, m = {max_anchor_mutations}, rules = {rules}: {pruned} < {unpruned}, reference = {}, query = {}",
                    String::from_utf8_lossy(reference),
                    String::from_utf8_lossy(query),
                );
                if rules == "none" {
                    assert_eq!(pruned, unpruned);
                }
            }
        }
    }
}

/// Without a guard for the length of template switches, the secondary rules remove the only secondary anchors that
/// are long enough for the template switch of this input, and it is lost (cost 5 -> 18).
/// This is accepted for the inexact mode (issue #236). The primary rules keep it.
#[test]
fn twin_in_place_inversion_loses_its_template_switch_with_secondary_rules() {
    let reference = b"AAAAAACCTTTTTCCCCCCCC";
    let query = b"AAAAAAGGAAAAAGGCCCCCC";
    let lower_bounds = lower_bounds(9, 1, 32);
    let run = |anchor_pruning, rules: &str| {
        chain_cost(
            reference,
            query,
            &lower_bounds,
            anchor_pruning,
            rules.parse().unwrap(),
        )
    };
    assert_eq!(run(AnchorPruning::None, "none"), (5.0, 1.0));
    assert_eq!(run(AnchorPruning::Exact, "none"), (5.0, 1.0));
    assert_eq!(run(AnchorPruning::Inexact, "primary"), (5.0, 1.0));
    assert_eq!(run(AnchorPruning::Inexact, "all"), (18.0, 0.0));
}
