use generic_a_star::cost::U32Cost;
use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;

use crate::{
    alignment::sequences::AlignmentSequences,
    anchors::{
        Anchors,
        domination::{is_pruning_sound, prune_dominated_primary_anchors},
        primary::PrimaryAnchor,
    },
    chain_align::performance_parameters::AlignmentPerformanceParameters,
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    dna_rc_fn,
};

fn primary_costs() -> GapAffineCosts<U32Cost> {
    GapAffineCosts::new(2u32, 3, 1).into_costs()
}

/// The costs of `test_files/config/chainalign`.
fn chainalign_config_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        primary_costs(),
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

fn anchor(
    seq1: std::ops::Range<usize>,
    seq2: std::ops::Range<usize>,
    cost: u32,
) -> PrimaryAnchor<U32Cost> {
    PrimaryAnchor::new_from_ranges(seq1, seq2, cost.into())
}

fn pruned(
    mut anchors: Vec<PrimaryAnchor<U32Cost>>,
    seq1: &[u8],
    seq2: &[u8],
    k: usize,
) -> Vec<PrimaryAnchor<U32Cost>> {
    prune_dominated_primary_anchors(&mut anchors, seq1, seq2, k, &primary_costs());
    anchors
}

#[test]
fn pruning_is_restricted_to_the_proven_cases() {
    assert!(is_pruning_sound(1, &primary_costs()));
    // No candidates with exact anchors, and no proof for more than one mutation.
    assert!(!is_pruning_sound(0, &primary_costs()));
    assert!(!is_pruning_sound(2, &primary_costs()));
    assert!(!is_pruning_sound(3, &primary_costs()));
    // Merging the removed column with a following gap must not increase the cost.
    assert!(is_pruning_sound(
        1,
        &GapAffineCosts::new(2u32, 3, 3).into_costs::<U32Cost>()
    ));
    assert!(!is_pruning_sound(
        1,
        &GapAffineCosts::new(2u32, 3, 4).into_costs::<U32Cost>()
    ));
}

/// An anchor with a trailing gap is dominated by the exact anchor it extends.
#[test]
fn trailing_gap_is_dominated_by_exact_anchor() {
    let seq1 = b"ACGTA";
    let seq2 = b"ACGT";
    let exact = anchor(0..4, 0..4, 0);
    let long = anchor(0..5, 0..4, 3);
    assert_eq!(pruned(vec![exact, long], seq1, seq2, 4), vec![exact]);

    // The same with the gap in the other sequence.
    let long = anchor(0..4, 0..5, 3);
    assert_eq!(pruned(vec![exact, long], seq2, seq1, 4), vec![exact]);
}

/// An anchor with a gap in its middle is dominated by the anchor without its last match.
#[test]
fn trailing_match_is_dominated_by_shorter_inexact_anchor() {
    // ACGGT vs AC-GT.
    let seq1 = b"ACGGT";
    let seq2 = b"ACGT";
    let short = anchor(0..4, 0..3, 3);
    let long = anchor(0..5, 0..4, 3);
    assert_eq!(pruned(vec![short, long], seq1, seq2, 4), vec![short]);
}

#[test]
fn anchor_without_dominator_is_kept() {
    let seq1 = b"ACGGT";
    let seq2 = b"ACGT";
    let long = anchor(0..5, 0..4, 3);

    // No other anchor.
    assert_eq!(pruned(vec![long], seq1, seq2, 4), vec![long]);
    // The shorter anchor starts elsewhere.
    let other_start = anchor(1..5, 0..3, 3);
    assert_eq!(
        pruned(vec![other_start, long], seq1, seq2, 4),
        vec![other_start, long]
    );
    // The shorter anchor is too expensive.
    let expensive = anchor(0..4, 0..3, 4);
    assert_eq!(
        pruned(vec![expensive, long], seq1, seq2, 4),
        vec![expensive, long]
    );
    // The last column would be a substitution, which makes the extension too expensive.
    let seq2 = b"ACGA";
    let short = anchor(0..4, 0..3, 3);
    assert_eq!(pruned(vec![short, long], seq1, seq2, 4), vec![short, long]);
}

/// Only anchors with a side of length k + 1 are pruned, even if a shorter anchor could be extended to them.
/// This ensures that dominators are never pruned themselves, and that anchors with zero cost are never pruned.
#[test]
fn only_anchors_with_a_side_of_length_k_plus_one_are_pruned() {
    let seq1 = b"ACGTA";
    let seq2 = b"ACGTA";
    // With k = 5, the anchor (0..4, 0..5) can be extended by a gap to (0..5, 0..5),
    // which would dominate an anchor at this range with cost 6. But the latter has no side of length 6.
    let short = anchor(0..4, 0..5, 3);
    let substitutions = anchor(0..5, 0..5, 6);
    assert_eq!(
        pruned(vec![short, substitutions], seq1, seq2, 5),
        vec![short, substitutions]
    );

    // A (hypothetical) anchor of length k + 1 on both sides has no dominator candidate without a side of length k + 1
    // other than the one that removes a column from both sides.
    let exact = anchor(0..4, 0..4, 0);
    let long_a = anchor(0..5, 0..4, 3);
    let both_long = anchor(0..5, 0..5, 1);
    assert_eq!(
        pruned(vec![exact, long_a, both_long], seq1, seq2, 4),
        vec![exact]
    );
    let both_long = anchor(0..5, 0..5, 3);
    let long_b = anchor(0..4, 0..5, 0);
    assert_eq!(
        pruned(vec![long_b, both_long], seq1, seq2, 4),
        vec![long_b, both_long]
    );
}

/// Anchors with zero cost are never pruned, even if they have a side of length k + 1 (which requires free gaps).
#[test]
fn zero_cost_anchors_are_kept() {
    let seq1 = b"ACGTA";
    let seq2 = b"ACGT";
    let exact = anchor(0..4, 0..4, 0);
    let long = anchor(0..5, 0..4, 0);
    let mut anchors = vec![exact, long];
    let free_gaps = GapAffineCosts::new(2u32, 0, 0).into_costs::<U32Cost>();
    assert_eq!(
        prune_dominated_primary_anchors(&mut anchors, seq1, seq2, 4, &free_gaps),
        0
    );
    assert_eq!(anchors, vec![exact, long]);
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
            // Sometimes replace a part of the query by its reverse complement, to create template switches.
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

/// With one allowed mutation, exactly the primary anchors with a side of length k + 1 are pruned,
/// because the anchors are complete.
#[test]
fn all_long_anchors_are_pruned_from_complete_anchors() {
    let costs = chainalign_config_costs();
    for (reference, query) in random_sequences(0xd0d0, 20) {
        for k in [3usize, 5] {
            let sequences = AlignmentSequences::new_complete(reference.clone(), query.clone());
            let all = Anchors::new_inexact(&sequences, k as u32, 1, &costs, &dna_rc_fn);
            let mut pruned_anchors =
                Anchors::new_inexact(&sequences, k as u32, 1, &costs, &dna_rc_fn);
            let pruned_amount =
                pruned_anchors.prune_dominated_primary_anchors(&sequences, k as u32, 1, &costs);

            let expected: Vec<_> = all
                .enumerate_primaries()
                .map(|(_, anchor)| anchor)
                .filter(|anchor| anchor.range().len_a() != k + 1 && anchor.range().len_b() != k + 1)
                .collect();
            let actual: Vec<_> = pruned_anchors
                .enumerate_primaries()
                .map(|(_, anchor)| anchor)
                .collect();
            assert_eq!(expected, actual);
            assert_eq!(
                pruned_amount,
                all.enumerate_primaries().count() - actual.len()
            );
            // Anchors with zero cost are never pruned.
            assert_eq!(
                all.enumerate_primaries()
                    .filter(|(_, anchor)| anchor.cost() == U32Cost::from(0u32))
                    .count(),
                actual
                    .iter()
                    .filter(|anchor| anchor.cost() == U32Cost::from(0u32))
                    .count()
            );
            // Secondary anchors are never pruned.
            assert_eq!(
                format!("{all}").lines().skip(1).collect::<Vec<_>>(),
                format!("{pruned_anchors}")
                    .lines()
                    .skip(1)
                    .collect::<Vec<_>>()
            );

            // Nothing is pruned with two allowed mutations.
            let mut two_mutations =
                Anchors::new_inexact(&sequences, k as u32, 2, &costs, &dna_rc_fn);
            assert_eq!(
                two_mutations.prune_dominated_primary_anchors(&sequences, k as u32, 2, &costs),
                0
            );
        }
    }
}

/// Returns the chain cost of the alignment with and without pruning dominated anchors.
fn chain_costs(
    reference: &[u8],
    query: &[u8],
    lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> (f64, f64) {
    let chain_cost = |prune_dominated_anchors| {
        let performance_parameters = AlignmentPerformanceParameters {
            max_exact_cost_function_cost: U32Cost::from(0u32),
            prune_dominated_anchors,
            ..Default::default()
        };
        crate::align(
            reference.to_vec(),
            query.to_vec(),
            AlignmentRange::new_complete(reference.len(), query.len()),
            &performance_parameters,
            &dna_rc_fn,
            "reference",
            "query",
            lower_bounds,
        )
        .statistics()
        .chain_cost
        .unwrap()
        .raw()
    };
    (chain_cost(true), chain_cost(false))
}

fn lower_bounds(k: u32, max_n: usize) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(
        max_n,
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations: 1,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        chainalign_config_costs(),
    )
}

/// Pruning dominated secondary anchors without restriction made the template switch of this input infeasible,
/// because the only secondary anchor that is long enough for the minimum length of a template switch
/// (10 = k + 1) is dominated. The cost increased from 5 to 18 (nine substitutions).
/// Secondary anchors must not be pruned.
#[test]
fn twin_in_place_inversion_keeps_its_template_switch() {
    let reference = b"AAAAAACCTTTTTCCCCCCCC";
    let query = b"AAAAAAGGAAAAAGGCCCCCC";
    let lower_bounds = lower_bounds(9, 32);
    assert_eq!(chain_costs(reference, query, &lower_bounds), (5.0, 5.0));
}

/// Pruning dominated anchors does not change the cost of the optimal chain.
///
/// This is a small smoke test, since chaining with one mutation and small k is slow without optimisations.
/// The end-to-end tests compare the chain costs with and without pruning for every case with inexact anchors.
#[test]
fn pruning_does_not_change_the_chain_cost() {
    let cases = random_sequences(0xfeed, 12);
    for k in [5, 6] {
        let lower_bounds = lower_bounds(k, 32);
        for (reference, query) in &cases {
            let (pruned_cost, unpruned_cost) = chain_costs(reference, query, &lower_bounds);
            assert_eq!(
                pruned_cost,
                unpruned_cost,
                "k = {k}, reference = {}, query = {}",
                String::from_utf8_lossy(reference),
                String::from_utf8_lossy(query),
            );
        }
    }
}
