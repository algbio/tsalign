use generic_a_star::cost::U32Cost;
use num_traits::{Bounded, Zero};

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates, sequences::AlignmentSequences, ts_kind::TsKind,
    },
    anchors::Anchors,
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
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
        max_match_run,
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
) -> (Anchors, ChainingCostFunction<U32Cost>) {
    let anchors = Anchors::new(sequences, max_match_run.checked_add(1).unwrap(), &dna_rc_fn);
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
fn test_start_end_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(5, 5),
        AlignmentCoordinates::new_primary(5, 5),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(cost_function.start_to_end().is_zero());
}

#[test]
fn test_start_anchor_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        cost_function
            .primary_from_start(
                anchors
                    .primary_index_from_start_coordinates(sequences.primary_start())
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn test_anchor_end_direct() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        cost_function
            .primary_to_end(
                anchors
                    .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(4, 4))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn test_start_end_indirect_lt_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 8);
    assert!(cost_function.start_to_end().is_zero());
}

#[test]
fn test_start_end_indirect_geq_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (_, cost_function) = create_chaining_cost_function(&sequences, 7);
    assert!(!cost_function.start_to_end().is_zero());
}

#[test]
fn test_start_anchor_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        !cost_function
            .primary_from_start(
                anchors
                    .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(2, 2))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn test_anchor_end_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 4);
    assert!(
        !cost_function
            .primary_to_end(
                anchors
                    .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(3, 3))
                    .unwrap()
            )
            .is_zero()
    );
}

#[test]
fn test_anchor_anchor_direct_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.primary(
            anchors
                .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(1, 1))
                .unwrap(),
            anchors
                .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(4, 4))
                .unwrap(),
        ),
        U32Cost::max_value(),
    );
}

#[test]
fn test_anchor_anchor_indirect_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.primary(
            anchors
                .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(1, 1))
                .unwrap(),
            anchors
                .primary_index_from_start_coordinates(AlignmentCoordinates::new_primary(5, 5))
                .unwrap(),
        ),
        2u8.into(),
    );
}

#[test]
fn test_anchor_anchor_direct_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.secondary(
            anchors
                .secondary_index_from_start_coordinates(AlignmentCoordinates::new_secondary(
                    9,
                    1,
                    TsKind::TS12
                ))
                .unwrap(),
            anchors
                .secondary_index_from_start_coordinates(AlignmentCoordinates::new_secondary(
                    6,
                    4,
                    TsKind::TS12
                ))
                .unwrap(),
            TsKind::TS12,
        ),
        U32Cost::max_value(),
    );
}

#[test]
fn test_anchor_anchor_indirect_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        AlignmentCoordinates::new_primary(1, 1),
        AlignmentCoordinates::new_primary(9, 9),
    );
    let (anchors, cost_function) = create_chaining_cost_function(&sequences, 2);
    assert_eq!(
        cost_function.secondary(
            anchors
                .secondary_index_from_start_coordinates(AlignmentCoordinates::new_secondary(
                    9,
                    1,
                    TsKind::TS12
                ))
                .unwrap(),
            anchors
                .secondary_index_from_start_coordinates(AlignmentCoordinates::new_secondary(
                    5,
                    5,
                    TsKind::TS12
                ))
                .unwrap(),
            TsKind::TS12,
        ),
        4u8.into(),
    );
}
