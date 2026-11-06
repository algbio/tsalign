use std::collections::BTreeMap;

use generic_a_star::cost::U32Cost;
use num_traits::bounds::UpperBounded;

use crate::{
    alignment::{
        Alignment, AlignmentType,
        coordinates::{PrimaryAlignmentCoordinates, SpecificSecondaryAlignmentCoordinates},
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    alignment_history::{
        history_vec::UnsignedIntAlignmentHistoryVec,
        test_util::{Rng, alignment_segments, contains_anchor},
    },
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
    windowed_min_mutation_chaining::ts_12_jump::Ts12JumpAligner,
};

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => unimplemented!("Unsupported character {c}"),
    }
}

#[test]
fn start_end() {
    let seq1 = b"AAGG".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(0, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 5, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        50,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -1,
                    ts_kind: TsKind::TS12
                }
            ),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(4u8));
}

#[test]
fn partial_alignment() {
    let seq1 = b"AAGG".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(1, 1);
    let end = SpecificSecondaryAlignmentCoordinates::new(1, 4, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        50,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -1,
                    ts_kind: TsKind::TS12
                }
            ),
            (1, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(4u8));
}

#[test]
fn gap_directions() {
    let seq1 = b"CCCCCCACCAACAAAAAA".to_vec();
    let seq2 = b"AAAAAACAAGGGGGGAGG".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(10u8),
            U32Cost::from(1u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(10u8),
            U32Cost::from(1u8),
            U32Cost::from(1u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(9, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 18, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        50,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (4, AlignmentType::Match),
            (1, AlignmentType::GapA),
            (2, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -9,
                    ts_kind: TsKind::TS12
                }
            ),
            (2, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (4, AlignmentType::Match),
            (1, AlignmentType::GapA),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(6u8));
}

#[test]
fn max_match_run_0() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"AAAAAAAACCTCCTCC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(8, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        1,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (16, AlignmentType::GapA),
            (
                1,
                AlignmentType::TsStart {
                    jump: -8,
                    ts_kind: TsKind::TS12
                }
            ),
        ]
    );
    assert_eq!(cost, U32Cost::from(20u8));
}

#[test]
fn max_match_run_1() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"AAAAAAAACCTCCTCC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(8, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        2,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (5, AlignmentType::GapA),
            (1, AlignmentType::Match),
            (7, AlignmentType::GapA),
            (1, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -10,
                    ts_kind: TsKind::TS12
                }
            ),
            (1, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(18u8));
}

#[test]
fn max_match_run_2() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"AAAAAAAACCCCCCCC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(30u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(60u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(8, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -8,
                    ts_kind: TsKind::TS12
                }
            ),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(14u8));
}

#[test]
fn only_jump() {
    let seq1 = b"ACGTACGTAC".to_vec();
    let seq2 = b"ACGTACGTAC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(5, 6);
    let end = SpecificSecondaryAlignmentCoordinates::new(3, 6, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(
            1,
            AlignmentType::TsStart {
                jump: -2,
                ts_kind: TsKind::TS12
            }
        ),]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn only_jump_start() {
    let seq1 = b"ACGTACGTAC".to_vec();
    let seq2 = b"ACGTACGTAC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(0, 0);
    let end = SpecificSecondaryAlignmentCoordinates::new(0, 0, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(
            1,
            AlignmentType::TsStart {
                jump: 0,
                ts_kind: TsKind::TS12
            }
        ),]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn only_jump_end() {
    let seq1 = b"ACGTACGTAC".to_vec();
    let seq2 = b"ACGTACGTAC".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = PrimaryAlignmentCoordinates::new(10, 10);
    let end = SpecificSecondaryAlignmentCoordinates::new(10, 10, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(
            1,
            AlignmentType::TsStart {
                jump: 0,
                ts_kind: TsKind::TS12
            }
        ),]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

const PARETO_TEST_COSTS: [(u8, u8, u8); 5] =
    [(2, 3, 1), (1, 1, 1), (3, 2, 1), (1, 4, 2), (2, 5, 0)];

fn random_alignment_costs(rng: &mut Rng) -> AlignmentCosts<U32Cost> {
    let mut random_gap_affine_costs = || {
        let (substitution, gap_open, gap_extend) =
            PARETO_TEST_COSTS[rng.below(PARETO_TEST_COSTS.len())];
        GapAffineCosts::new(
            U32Cost::from(substitution),
            U32Cost::from(gap_open),
            U32Cost::from(gap_extend),
        )
    };
    let primary_costs = random_gap_affine_costs();
    let secondary_costs = random_gap_affine_costs();

    AlignmentCosts {
        primary_costs,
        secondary_costs,
        ts_base_cost: U32Cost::from(rng.below(4) as u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    }
}

/// Returns the lengths of the ancestor and the descendant of the template switch.
fn ts_lengths(ts_kind: TsKind, seq1: &[u8], seq2: &[u8]) -> (usize, usize) {
    match ts_kind {
        TsKind::TS11 => (seq1.len(), seq1.len()),
        TsKind::TS12 => (seq1.len(), seq2.len()),
        TsKind::TS21 => (seq2.len(), seq1.len()),
        TsKind::TS22 => (seq2.len(), seq2.len()),
    }
}

/// Returns the minimum cost of each coordinate whose minimum cost is at most `cost_limit`.
fn min_costs_up_to<Coordinates: Ord + Copy>(
    targets: &[(Coordinates, U32Cost)],
    cost_limit: U32Cost,
) -> BTreeMap<Coordinates, U32Cost> {
    let mut result = BTreeMap::new();
    for (coordinates, cost) in targets {
        let entry = result.entry(*coordinates).or_insert(*cost);
        *entry = (*entry).min(*cost);
    }
    result.retain(|_, cost| *cost <= cost_limit);
    result
}

fn assert_anchor_free(
    alignment: &Alignment,
    anchor_k: u8,
    max_anchor_mutations: u8,
    description: &str,
) {
    for segment in alignment_segments(alignment) {
        assert!(
            !contains_anchor(&segment, anchor_k.into(), max_anchor_mutations.into()),
            "{description}"
        );
    }
}

#[test]
fn pareto_pruning_does_not_change_alignments() {
    let mut rng = Rng::new(3);

    for _ in 0..60 {
        let len = 4 + rng.below(12);
        let seq1 = rng.dna(len);
        let mutation_percent = 10 + rng.below(30);
        let seq2 = rng.mutate(&seq1, mutation_percent);
        let anchor_k = [3, 4, 5][rng.below(3)];
        let max_anchor_mutations = rng.below(usize::from(anchor_k).min(3)) as u8;
        let cost_table = random_alignment_costs(&mut rng);
        let ts_kind = TsKind::iter().nth(rng.below(4)).unwrap();
        let (ancestor_len, descendant_len) = ts_lengths(ts_kind, &seq1, &seq2);
        let start = PrimaryAlignmentCoordinates::new(
            rng.below(seq1.len() / 3 + 1),
            rng.below(seq2.len() / 3 + 1),
        );
        let end = SpecificSecondaryAlignmentCoordinates::new(
            rng.below(ancestor_len / 2 + 1),
            descendant_len - rng.below(descendant_len / 3 + 1),
            ts_kind,
        );
        let cost_limit = U32Cost::from(rng.below(12) as u8);
        let description = format!(
            "seq1 = {}, seq2 = {}, start = {start}, end = {end}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}",
            String::from_utf8_lossy(&seq1),
            String::from_utf8_lossy(&seq2),
        );

        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let mut aligners = [true, false].map(|pareto_pruning| {
            let mut aligner = Ts12JumpAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
                &sequences,
                &cost_table,
                &rc_fn,
                anchor_k,
                max_anchor_mutations,
            );
            aligner.set_pareto_pruning(pareto_pruning);
            aligner
        });

        let [
            (cost, alignment, targets),
            (unpruned_cost, _, unpruned_targets),
        ] = aligners.each_mut().map(|aligner| {
            let mut targets = Vec::new();
            let (cost, alignment) = aligner.align(start, end, &mut targets);
            (cost, alignment, targets)
        });

        assert_eq!(cost, unpruned_cost, "{description}");
        if cost != U32Cost::max_value() {
            assert_anchor_free(&alignment, anchor_k, max_anchor_mutations, &description);
            assert_eq!(
                min_costs_up_to(&targets, cost),
                min_costs_up_to(&unpruned_targets, cost),
                "{description}"
            );
        }

        // Reuse the aligners to check that no state leaks between alignments.
        let [targets, unpruned_targets] = aligners.each_mut().map(|aligner| {
            let mut targets = Vec::new();
            aligner.align_until_cost_limit(start, end, cost_limit, &mut targets);
            targets
        });
        assert_eq!(
            min_costs_up_to(&targets, cost_limit),
            min_costs_up_to(&unpruned_targets, cost_limit),
            "{description}, cost_limit = {cost_limit}"
        );
    }
}
