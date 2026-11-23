use generic_a_star::cost::U32Cost;
use num_traits::{Bounded, Zero};

use crate::{
    alignment::{
        coordinates::{
            AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            range::{AnySecondaryAlignmentRange, PrimaryAlignmentRange},
        },
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    anchors::Anchors,
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
};

fn dna_rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => panic!("Unsupported character: {c}"),
    }
}

fn create_lower_bounds(max_match_run: u32) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(
        20,
        ChainingLowerBoundConfig {
            anchor_k: max_match_run + 1,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        AlignmentCosts::new(
            GapAffineCosts::new(2u8.into(), 3u8.into(), 1u8.into()),
            GapAffineCosts::new(4u8.into(), 6u8.into(), 2u8.into()),
            U32Cost::from(2u8).into(),
            TsLimits::new_unlimited(),
        ),
    )
}

fn create_chaining_cost_function(
    sequences: &AlignmentSequences,
    max_match_run: u32,
) -> (Anchors<U32Cost>, ChainingCostFunction<U32Cost>) {
    let anchors = Anchors::new_exact(sequences, max_match_run.checked_add(1).unwrap(), &dna_rc_fn);
    let cost_function = ChainingCostFunction::new_from_lower_bounds(
        &create_lower_bounds(max_match_run),
        &anchors,
        sequences,
        U32Cost::zero(),
        &dna_rc_fn,
    );
    (anchors, cost_function)
}

#[test]
fn start_end_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(5, 5),
        PrimaryAlignmentCoordinates::new(5, 5),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(cost_function.start_to_end().is_zero());
}

#[test]
fn start_anchor_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        cost_function
            .primary_from_start(
                anchors
                    .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                        sequences.primary_start(),
                        5
                    ))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn anchor_end_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        cost_function
            .primary_to_end(
                anchors
                    .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                        PrimaryAlignmentCoordinates::new(4, 4),
                        5
                    ))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn start_end_indirect_lt_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 8);
    assert!(cost_function.start_to_end().is_zero());
}

#[test]
fn start_end_indirect_geq_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 7);
    assert!(!cost_function.start_to_end().is_zero());
}

#[test]
fn start_anchor_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        !cost_function
            .primary_from_start(
                anchors
                    .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                        PrimaryAlignmentCoordinates::new(2, 2),
                        5
                    ))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn anchor_end_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        !cost_function
            .primary_to_end(
                anchors
                    .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                        PrimaryAlignmentCoordinates::new(3, 3),
                        5
                    ))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn anchor_anchor_direct_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.primary(
            anchors
                .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                    PrimaryAlignmentCoordinates::new(1, 1),
                    3
                ))
                .unwrap(),
            anchors
                .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                    PrimaryAlignmentCoordinates::new(4, 4),
                    3
                ))
                .unwrap(),
        ),
        U32Cost::max_value(),
    );
}

#[test]
fn anchor_anchor_indirect_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.primary(
            anchors
                .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                    PrimaryAlignmentCoordinates::new(1, 1),
                    3
                ))
                .unwrap(),
            anchors
                .primary_index_from_range(PrimaryAlignmentRange::new_equal_length(
                    PrimaryAlignmentCoordinates::new(5, 5),
                    3
                ))
                .unwrap(),
        ),
        2u8.into(),
    );
}

#[test]
fn anchor_anchor_direct_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.secondary(
            anchors
                .secondary_index_from_range(
                    AnySecondaryAlignmentRange::new_equal_length(
                        AnySecondaryAlignmentCoordinates::new(9, 1),
                        3
                    )
                    .into_specific(TsKind::TS12)
                )
                .unwrap(),
            anchors
                .secondary_index_from_range(
                    AnySecondaryAlignmentRange::new_equal_length(
                        AnySecondaryAlignmentCoordinates::new(6, 4),
                        3
                    )
                    .into_specific(TsKind::TS12)
                )
                .unwrap(),
            TsKind::TS12,
        ),
        U32Cost::max_value(),
    );
}

#[test]
fn anchor_anchor_indirect_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.secondary(
            anchors
                .secondary_index_from_range(
                    AnySecondaryAlignmentRange::new_equal_length(
                        AnySecondaryAlignmentCoordinates::new(9, 1),
                        3
                    )
                    .into_specific(TsKind::TS12)
                )
                .unwrap(),
            anchors
                .secondary_index_from_range(
                    AnySecondaryAlignmentRange::new_equal_length(
                        AnySecondaryAlignmentCoordinates::new(5, 5),
                        3
                    )
                    .into_specific(TsKind::TS12)
                )
                .unwrap(),
            TsKind::TS12,
        ),
        4u8.into(),
    );
}

/// Costs of the `test_files/config/chainalign` config with unlimited template switch geometry.
fn chainalign_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u8.into(), 3u8.into(), 1u8.into()),
        GapAffineCosts::new(4u8.into(), 3u8.into(), 2u8.into()),
        U32Cost::from(2u8).into(),
        TsLimits::new_unlimited(),
    )
}

/// Create the chaining cost function like `crate::align` does for exact anchors with `max_exact_cost_function_cost = 0`.
fn create_chainalign_cost_function(
    sequences: &AlignmentSequences,
    k: u32,
) -> (Anchors<U32Cost>, ChainingCostFunction<U32Cost>) {
    let anchors = Anchors::new_exact(sequences, k, &dna_rc_fn);
    let lower_bounds = ChainingLowerBounds::new(
        sequences.seq1().len().max(sequences.seq2().len()),
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        chainalign_costs(),
    );
    let cost_function = ChainingCostFunction::new_from_lower_bounds(
        &lower_bounds,
        &anchors,
        sequences,
        U32Cost::zero(),
        &dna_rc_fn,
    );
    (anchors, cost_function)
}

/// Reference `CC` + `GATTACAGGCTTAC`, query `CC` + reverse complement: the optimal alignment is
/// `2=[TS]` with cost 2, where the template switch starts directly after the `CC` at the start.
///
/// For k = 4 and k = 6, the `CC` is shorter than k, so there is no primary anchor before the template switch,
/// and the template switch must be entered by a 12-jump from the start of the alignment.
/// Its first secondary anchor is `SA((16, 16 - k], [2, 2 + k))` of kind 21 (descendant seq1, ancestor seq2).
/// The 12-jump from the start to this anchor is `2=` followed by the jump, which costs only the template switch
/// base cost 2, so the cost function entry must be exactly 2.
///
/// The 12-jump aligner used to require a non-match in the jump because the descendant gap (2) is non-empty,
/// so the exact value was not found within the cost limit, and the entry was set to the lower bound,
/// which is raised to `max_exact_cost_function_cost + ts_base_cost + 1 = 3`.
#[test]
fn jump_12_from_start_after_matches_shorter_than_k() {
    let sequences = AlignmentSequences::new_complete(
        b"CCGATTACAGGCTTAC".to_vec(),
        b"CCGTAAGCCTGTAATC".to_vec(),
    );
    for k in [4, 6] {
        let (anchors, cost_function) = create_chainalign_cost_function(&sequences, k);
        let index = anchors
            .secondary_index_from_range(
                AnySecondaryAlignmentRange::new_equal_length(
                    AnySecondaryAlignmentCoordinates::new(16, 2),
                    k as usize,
                )
                .into_specific(TsKind::TS21),
            )
            .unwrap();
        assert_eq!(
            (
                cost_function.jump_12_from_start(index, TsKind::TS21),
                cost_function.is_jump_12_from_start_exact(index, TsKind::TS21)
            ),
            (U32Cost::from(2u8), true),
            "k = {k}"
        );
    }
}

/// Reference `GATTACAGGCTTAC` + `CC`, query reverse complement + `CC`: the optimal alignment is
/// `[TS]2=` with cost 2, where the template switch ends directly before the `CC` at the end.
///
/// For k = 4 and k = 6, the `CC` is shorter than k, so there is no primary anchor after the template switch,
/// and the template switch must be exited by a 34-jump to the end of the alignment.
/// Its last secondary anchor is `SA((k, 0], [14 - k, 14))` of kind 21 (descendant seq1, ancestor seq2).
/// The 34-jump from this anchor to the end is the jump followed by `2=`, which costs 0
/// (the template switch base cost is counted in the 12-jump), so the cost function entry must be exactly 0.
///
/// The 34-jump aligner used to require a non-match in the jump because the descendant gap (2) is non-empty,
/// so the exact value was not found within the cost limit, and the entry was set to the lower bound,
/// which is raised to `max_exact_cost_function_cost + 1 = 1`.
#[test]
fn jump_34_to_end_before_matches_shorter_than_k() {
    let sequences = AlignmentSequences::new_complete(
        b"GATTACAGGCTTACCC".to_vec(),
        b"GTAAGCCTGTAATCCC".to_vec(),
    );
    for k in [4, 6] {
        let (anchors, cost_function) = create_chainalign_cost_function(&sequences, k);
        let index = anchors
            .secondary_index_from_range(
                AnySecondaryAlignmentRange::new_equal_length(
                    AnySecondaryAlignmentCoordinates::new(k as usize, 14 - k as usize),
                    k as usize,
                )
                .into_specific(TsKind::TS21),
            )
            .unwrap();
        assert_eq!(
            (
                cost_function.jump_34_to_end(index, TsKind::TS21),
                cost_function.is_jump_34_to_end_exact(index, TsKind::TS21)
            ),
            (U32Cost::from(0u8), true),
            "k = {k}"
        );
    }
}
