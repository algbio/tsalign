use generic_a_star::cost::U32Cost;
use num_traits::bounds::UpperBounded;

use crate::{
    alignment::{
        AlignmentType,
        coordinates::{PrimaryAlignmentCoordinates, SpecificSecondaryAlignmentCoordinates},
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
    max_match_run_chaining::ts_34_jump::Ts34JumpAligner,
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

    let start = SpecificSecondaryAlignmentCoordinates::new(2, 0, TsKind::TS12);
    let end = PrimaryAlignmentCoordinates::new(4, 5);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
        .without_enforced_non_match_next_to_anchors();
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
fn partial_alignment() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(1, 1, TsKind::TS12);
    let end = PrimaryAlignmentCoordinates::new(3, 4);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
        .without_enforced_non_match_next_to_anchors();
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

    let start = SpecificSecondaryAlignmentCoordinates::new(18, 0, TsKind::TS21);
    let end = PrimaryAlignmentCoordinates::new(18, 9);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
        .without_enforced_non_match_next_to_anchors();
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
fn max_match_run_0() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(8, 0, TsKind::TS12);
    let end = PrimaryAlignmentCoordinates::new(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 0, 0)
        .without_enforced_non_match_next_to_anchors();
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    // Multiple alignments are optimal.
    let expected_alignments = [
        vec![
            (1, AlignmentType::TsEnd { jump: 6 }),
            (2, AlignmentType::Substitution),
            (14, AlignmentType::GapA),
        ],
        vec![
            (1, AlignmentType::TsEnd { jump: 7 }),
            (1, AlignmentType::Substitution),
            (15, AlignmentType::GapA),
        ],
        vec![
            (1, AlignmentType::TsEnd { jump: 8 }),
            (16, AlignmentType::GapA),
        ],
    ];
    assert!(
        expected_alignments.contains(&alignment.alignment),
        "{:?}",
        alignment.alignment
    );
    assert_eq!(cost, U32Cost::from(18u8));
}

#[test]
fn max_match_run_1() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(8, 0, TsKind::TS12);
    let end = PrimaryAlignmentCoordinates::new(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 1, 0)
        .without_enforced_non_match_next_to_anchors();
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
fn max_match_run_2() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(8, 0, TsKind::TS12);
    let end = PrimaryAlignmentCoordinates::new(16, 16);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
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
fn only_jump() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(2, 8, TsKind::TS21);
    let end = PrimaryAlignmentCoordinates::new(8, 8);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 6 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}

#[test]
fn only_jump_start() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(0, 0, TsKind::TS21);
    let end = PrimaryAlignmentCoordinates::new(0, 0);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 0 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}

#[test]
fn only_jump_end() {
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

    let start = SpecificSecondaryAlignmentCoordinates::new(10, 10, TsKind::TS21);
    let end = PrimaryAlignmentCoordinates::new(10, 10);
    let mut aligner = Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::TsEnd { jump: 0 })]
    );
    assert_eq!(cost, U32Cost::from(0u8));
}

/// With exact anchors, the parts of the alignment that touch an anchor must be non-matches next to the anchor.
/// The jump itself counts as non-match, and the start and end of the alignment are no anchors.
#[test]
fn enforce_non_match_next_to_anchors() {
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
    let is_non_match_or_jump = |alignment_type: &AlignmentType| {
        matches!(
            alignment_type,
            AlignmentType::Substitution
                | AlignmentType::GapA
                | AlignmentType::GapB
                | AlignmentType::TsEnd { .. }
        )
    };

    for (seq1, seq2) in [
        (b"ACCGTA".to_vec(), b"TACGGT".to_vec()),
        (b"AAAAAA".to_vec(), b"TTTTT".to_vec()),
        (b"ACGTAC".to_vec(), b"GTACGT".to_vec()),
    ] {
        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let primary_end = sequences.primary_end();

        for primary_a in 0..=seq1.len() {
            for primary_b in 0..=seq2.len() {
                for ancestor in 0..=seq1.len() {
                    for descendant in 0..=seq2.len() {
                        let primary = PrimaryAlignmentCoordinates::new(primary_a, primary_b);
                        let secondary = SpecificSecondaryAlignmentCoordinates::new(
                            ancestor,
                            descendant,
                            TsKind::TS12,
                        );
                        if primary_b < descendant {
                            continue;
                        }
                        let (start, end) = (secondary, primary);
                        let start_is_anchor = true;
                        let end_is_anchor = end != primary_end;

                        let mut enforced =
                            Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0);
                        let mut unenforced =
                            Ts34JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
                                .without_enforced_non_match_next_to_anchors();
                        let (cost, alignment) = enforced.align(start, end, &mut Vec::new());
                        let (unenforced_cost, _) = unenforced.align(start, end, &mut Vec::new());
                        let context = format!(
                            "{} vs {} from {start} to {end}: {alignment:?}",
                            String::from_utf8_lossy(&seq1),
                            String::from_utf8_lossy(&seq2)
                        );

                        assert!(cost >= unenforced_cost, "{context}");
                        if cost == U32Cost::max_value() {
                            continue;
                        }

                        let first = alignment.alignment.first().unwrap().1;
                        let last = alignment.alignment.last().unwrap().1;
                        if start_is_anchor {
                            assert!(is_non_match_or_jump(&first), "{context}");
                        }
                        if end_is_anchor {
                            assert!(is_non_match_or_jump(&last), "{context}");
                        }
                    }
                }
            }
        }
    }
}
