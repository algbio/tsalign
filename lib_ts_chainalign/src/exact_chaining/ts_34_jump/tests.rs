use generic_a_star::cost::U32Cost;

use crate::{
    alignment::{
        AlignmentType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
    exact_chaining::ts_34_jump::Ts34JumpAligner,
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
    let seq2 = b"TTACG".to_vec();
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

    let start: AlignmentCoordinates = AlignmentCoordinates::new_secondary(2, 0, TsKind::TS12);
    let end = AlignmentCoordinates::new_primary(4, 5);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (1, AlignmentType::TsEnd { jump: 1 }),
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn test_partial_alignment() {
    let seq1 = b"AAGG".to_vec();
    let seq2 = b"TTACG".to_vec();
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

    let start: AlignmentCoordinates = AlignmentCoordinates::new_secondary(1, 1, TsKind::TS12);
    let end = AlignmentCoordinates::new_primary(3, 4);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::TsEnd { jump: 1 }),
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
        ]
    );
    assert_eq!(cost, U32Cost::from(2u8));
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

    let start = AlignmentCoordinates::new_secondary(18, 0, TsKind::TS21);
    let end = AlignmentCoordinates::new_primary(18, 9);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (4, AlignmentType::Match),
            (1, AlignmentType::GapA),
            (2, AlignmentType::Match),
            (1, AlignmentType::TsEnd { jump: -9 }),
            (2, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (4, AlignmentType::Match),
            (1, AlignmentType::GapA),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(4u8));
}

#[test]
fn test_max_match_run_0() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"CCCCCCCCAAAAAAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(8, 0, TsKind::TS12);
    let end = AlignmentCoordinates::new_primary(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 0);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::TsEnd { jump: 8 }),
            (16, AlignmentType::GapA),
        ]
    );
    assert_eq!(cost, U32Cost::from(18u8));
}

#[test]
fn test_max_match_run_1() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"CCCCCCCCAAAAAAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(8, 0, TsKind::TS12);
    let end = AlignmentCoordinates::new_primary(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 1);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::TsEnd { jump: -2 }),
            (5, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (4, AlignmentType::GapA),
            (1, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(13u8));
}

#[test]
fn test_max_match_run_2() {
    let seq1 = b"GGAGGAGGAACAACAA".to_vec();
    let seq2 = b"CCCCCCCCAAAAAAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(8, 0, TsKind::TS12);
    let end = AlignmentCoordinates::new_primary(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (1, AlignmentType::TsEnd { jump: 8 }),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(4u8));
}

#[test]
fn test_only_jump() {
    let seq1 = b"ACATCTGCAA".to_vec();
    let seq2 = b"ACGCAGATAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(2, 8, TsKind::TS21);
    let end = AlignmentCoordinates::new_primary(8, 8);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 6 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}

#[test]
fn test_only_jump_start() {
    let seq1 = b"ACATCTGCAA".to_vec();
    let seq2 = b"ACGCAGATAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(0, 0, TsKind::TS21);
    let end = AlignmentCoordinates::new_primary(0, 0);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 0 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}

#[test]
fn test_only_jump_end() {
    let seq1 = b"ACATCTGCAA".to_vec();
    let seq2 = b"ACGCAGATAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(1u8),
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

    let start = AlignmentCoordinates::new_secondary(10, 10, TsKind::TS21);
    let end = AlignmentCoordinates::new_primary(10, 10);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 0 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}
