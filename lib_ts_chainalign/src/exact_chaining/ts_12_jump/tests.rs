use generic_a_star::cost::U32Cost;

use crate::{
    alignment::{
        AlignmentType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
    exact_chaining::ts_12_jump::Ts12JumpAligner,
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
fn test_start_end() {
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

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_secondary(0, 5, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
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
fn test_partial_alignment() {
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

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_secondary(1, 4, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
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
fn test_gap_directions() {
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

    let start = AlignmentCoordinates::new_primary(9, 0);
    let end = AlignmentCoordinates::new_secondary(0, 18, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
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
fn test_max_match_run_0() {
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

    let start = AlignmentCoordinates::new_primary(8, 0);
    let end = AlignmentCoordinates::new_secondary(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 0);
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
fn test_max_match_run_1() {
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

    let start = AlignmentCoordinates::new_primary(8, 0);
    let end = AlignmentCoordinates::new_secondary(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 1);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (14, AlignmentType::GapA),
            (
                1,
                AlignmentType::TsStart {
                    jump: -8,
                    ts_kind: TsKind::TS12
                }
            ),
            (1, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(18u8));
}

#[test]
fn test_max_match_run_2() {
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

    let start = AlignmentCoordinates::new_primary(8, 0);
    let end = AlignmentCoordinates::new_secondary(0, 16, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
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
fn test_only_jump() {
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

    let start = AlignmentCoordinates::new_primary(5, 6);
    let end = AlignmentCoordinates::new_secondary(3, 6, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
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
fn test_only_jump_start() {
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

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_secondary(0, 0, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
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
fn test_only_jump_end() {
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

    let start = AlignmentCoordinates::new_primary(10, 10);
    let end = AlignmentCoordinates::new_secondary(10, 10, TsKind::TS12);
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
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
