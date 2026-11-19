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
    max_match_run_chaining::ts_12_jump::Ts12JumpAligner,
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
        .without_enforced_non_match_next_to_anchors();
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
        .without_enforced_non_match_next_to_anchors();
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX, 0)
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 0, 0)
        .without_enforced_non_match_next_to_anchors();
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 1, 0)
        .without_enforced_non_match_next_to_anchors();
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

    // Multiple alignments are optimal.
    let expected_alignments = [
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
                    ts_kind: TsKind::TS12,
                },
            ),
            (1, AlignmentType::Match),
        ],
        vec![
            (1, AlignmentType::Match),
            (14, AlignmentType::GapA),
            (
                1,
                AlignmentType::TsStart {
                    jump: -8,
                    ts_kind: TsKind::TS12,
                },
            ),
            (1, AlignmentType::Match),
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
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
    let mut aligner = Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
        .without_enforced_non_match_next_to_anchors();
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
                | AlignmentType::TsStart { .. }
        )
    };

    for (seq1, seq2) in [
        (b"ACCGTA".to_vec(), b"TACGGT".to_vec()),
        (b"AAAAAA".to_vec(), b"TTTTT".to_vec()),
        (b"ACGTAC".to_vec(), b"GTACGT".to_vec()),
    ] {
        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let primary_start = sequences.primary_start();

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
                        if descendant < primary_b {
                            continue;
                        }
                        let (start, end) = (primary, secondary);
                        let start_is_anchor = start != primary_start;
                        let end_is_anchor = true;

                        let mut enforced =
                            Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0);
                        let mut unenforced =
                            Ts12JumpAligner::new(&sequences, &cost_table, &rc_fn, 2, 0)
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

/// Costs of the `test_files/config/chainalign` config with unlimited template switch geometry.
fn chainalign_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        U32Cost::from(2u8).into(),
        TsLimits::new_unlimited(),
    )
}

/// The `(max_match_run, enforcement_offset)` parameters of the aligners as used by the chaining cost function
/// for k = 4 with exact anchors and with one anchor mutation.
const K4_ALIGNER_PARAMETERS: [(u32, u8); 2] = [(3, 0), (2, 4)];

/// The reference is `CC` followed by an inner sequence, and the query is `CC` followed by the reverse complement
/// of the inner sequence.
/// The optimal alignment is `2=` followed by a template switch that covers the inner sequence completely (cost 2).
fn matches_then_template_switch() -> AlignmentSequences {
    AlignmentSequences::new_complete(b"CCGATTACAGGCTTAC".to_vec(), b"CCGTAAGCCTGTAATC".to_vec())
}

/// A 12-jump from the start of the alignment to the first secondary anchor of an optimal template switch,
/// with anchors of length k = 4 (exact, and with one mutation).
///
/// The secondary anchor `SA((16, 12], [2, 6))` of kind 21 starts directly at the template switch entrance:
/// the descendant (seq1) continues at position 2 in the secondary alignment,
/// and the ancestor (seq2) is read in reverse from position 16.
/// The optimal 12-jump is `2=` (the `CC` at the start, shorter than k so no primary anchor covers it)
/// followed directly by the jump, so it costs only the template switch base cost 2.
///
/// The start of the alignment is no anchor, so this primary part may consist of matches only.
/// `Ts12JumpAligner::align` used to enforce a non-match whenever the descendant gap is non-empty
/// (`enforce_non_match`, independently of `enforce_non_match_next_to_anchors`), also from the start,
/// which forced a deletion into the primary part (`1=1D` + jump, cost 5).
#[test]
fn matches_from_alignment_start_then_jump() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let start = sequences.primary_start();
    let end = SpecificSecondaryAlignmentCoordinates::new(16, 2, TsKind::TS21);

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligner = Ts12JumpAligner::new(
            &sequences,
            &costs,
            &rc_fn,
            max_match_run,
            enforcement_offset,
        );
        let (cost, alignment) = aligner.align(start, end, &mut Vec::new());

        assert_eq!(
            (cost, alignment.alignment),
            (
                U32Cost::from(2u8),
                vec![
                    (2, AlignmentType::Match),
                    (
                        1,
                        AlignmentType::TsStart {
                            jump: 14,
                            ts_kind: TsKind::TS21
                        }
                    ),
                ]
            ),
            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}"
        );
    }
}

/// Same as [`matches_from_alignment_start_then_jump`], but via `align_until_cost_limit`,
/// which the chaining cost function uses to compute the exact costs of 12-jumps from the start.
///
/// With cost limit 2 (the template switch base cost, i.e. `max_exact_cost_function_cost = 0`),
/// the anchor start `(16, 2)` must be reported with cost 2.
/// `fill_additional_targets` used to report only nodes with a non-match if the descendant gap is non-empty,
/// also from the start, so the optimal `2=` + jump was not reported.
#[test]
fn matches_from_alignment_start_then_jump_until_cost_limit() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let start = sequences.primary_start();
    let target = SpecificSecondaryAlignmentCoordinates::new(16, 2, TsKind::TS21).into_any();

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligner = Ts12JumpAligner::new(
            &sequences,
            &costs,
            &rc_fn,
            max_match_run,
            enforcement_offset,
        );
        let mut targets = Vec::new();
        aligner.align_until_cost_limit(
            start,
            sequences.secondary_end(TsKind::TS21),
            U32Cost::from(2u8),
            &mut targets,
        );

        let target_cost = targets
            .iter()
            .filter(|(coordinates, _)| *coordinates == target)
            .map(|(_, cost)| *cost)
            .min();
        assert_eq!(
            target_cost,
            Some(U32Cost::from(2u8)),
            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}: {targets:?}"
        );
    }
}

fn is_non_match(alignment_type: &AlignmentType) -> bool {
    matches!(
        alignment_type,
        AlignmentType::Substitution | AlignmentType::GapA | AlignmentType::GapB
    )
}

/// A 12-jump from an inner primary coordinate (the end of a primary anchor) to the start of a secondary anchor,
/// with a non-empty descendant gap and an all-match primary and secondary part.
///
/// From `(1, 1)`, an optimal 12-jump to `SA((14, ...], [4, ...))` of kind 21 is `1=` (the second `C`),
/// the jump to the template switch entrance `(16, 2)`, and `2=` in the secondary part, so it costs only
/// the template switch base cost 2.
/// (Another optimal one is `2=`, the jump to `(15, 3)` and `1=`.)
///
/// With inexact anchors (`enforcement_offset > 0`), there is no shifted anchor that could be chained for free
/// instead, so this all-match jump must be allowed (bugs (ii) and (iii) of the 203 branch).
/// With exact anchors (`enforcement_offset == 0`), a non-match is still enforced, also without
/// enforced non-matches next to anchors.
#[test]
fn matches_from_inner_start_then_jump_then_matches() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let start = PrimaryAlignmentCoordinates::new(1, 1);
    let end = SpecificSecondaryAlignmentCoordinates::new(14, 4, TsKind::TS21);

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligners = vec![(
            "default",
            Ts12JumpAligner::new(
                &sequences,
                &costs,
                &rc_fn,
                max_match_run,
                enforcement_offset,
            ),
        )];
        if enforcement_offset == 0 {
            aligners.push((
                "without enforced non-match next to anchors",
                Ts12JumpAligner::new(
                    &sequences,
                    &costs,
                    &rc_fn,
                    max_match_run,
                    enforcement_offset,
                )
                .without_enforced_non_match_next_to_anchors(),
            ));
        }

        for (name, mut aligner) in aligners {
            let (cost, alignment) = aligner.align(start, end, &mut Vec::new());
            let context = format!(
                "{name}, max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}: {alignment:?}"
            );

            let has_non_match = alignment
                .alignment
                .iter()
                .any(|(_, alignment_type)| is_non_match(alignment_type));
            if enforcement_offset > 0 {
                assert_eq!(cost, U32Cost::from(2u8), "{context}");
                assert!(!has_non_match, "{context}");
                assert!(
                    alignment
                        .alignment
                        .iter()
                        .any(|(_, alignment_type)| matches!(
                            alignment_type,
                            AlignmentType::TsStart { .. }
                        )),
                    "{context}"
                );
            } else {
                assert!(cost > U32Cost::from(2u8), "{context}");
                assert!(has_non_match, "{context}");
            }
        }
    }
}

/// Same as [`matches_from_inner_start_then_jump_then_matches`], but via `align_until_cost_limit`,
/// which the chaining cost function uses to compute the exact costs of 12-jumps between anchors.
///
/// With cost limit 2 (the template switch base cost), the anchor start `(14, 4)` must be reported with cost 2
/// for inexact anchors, and must not be reported for exact anchors.
#[test]
fn matches_from_inner_start_then_jump_then_matches_until_cost_limit() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let start = PrimaryAlignmentCoordinates::new(1, 1);
    let target = SpecificSecondaryAlignmentCoordinates::new(14, 4, TsKind::TS21).into_any();

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligner = Ts12JumpAligner::new(
            &sequences,
            &costs,
            &rc_fn,
            max_match_run,
            enforcement_offset,
        );
        let mut targets = Vec::new();
        aligner.align_until_cost_limit(
            start,
            sequences.secondary_end(TsKind::TS21),
            U32Cost::from(2u8),
            &mut targets,
        );

        let target_cost = targets
            .iter()
            .filter(|(coordinates, _)| *coordinates == target)
            .map(|(_, cost)| *cost)
            .min();
        assert_eq!(
            target_cost,
            (enforcement_offset > 0).then(|| U32Cost::from(2u8)),
            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}: {targets:?}"
        );
    }
}

/// A 12-jump with an empty descendant gap, whose target can be reached by the bare jump,
/// and also by alignments that contain only deletions (gaps in the descendant).
///
/// From `(2, 2)`, the target `(14, 2)` of kind 21 is reached by the bare jump with cost 2 (the base cost),
/// and e.g. by jumping to `(15, 2)` followed by a secondary deletion, or by a primary deletion followed by the jump,
/// both with cost 5.
/// `align` accepts any alignment as target if the descendant gap is empty, so `fill_additional_targets`
/// reports the deletion-only alignments as well (they used to be dropped).
/// Since the bare jump dominates them, the minimum reported cost is the same as the cost of `align`.
#[test]
fn empty_descendant_gap_with_only_deletions() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let start = PrimaryAlignmentCoordinates::new(2, 2);
    let end = SpecificSecondaryAlignmentCoordinates::new(14, 2, TsKind::TS21);
    let target = end.into_any();

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let context =
            format!("max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}");
        let mut aligner = Ts12JumpAligner::new(
            &sequences,
            &costs,
            &rc_fn,
            max_match_run,
            enforcement_offset,
        );
        let (cost, alignment) = aligner.align(start, end, &mut Vec::new());
        assert_eq!(
            (cost, alignment.alignment),
            (
                U32Cost::from(2u8),
                vec![(
                    1,
                    AlignmentType::TsStart {
                        jump: 12,
                        ts_kind: TsKind::TS21
                    }
                )]
            ),
            "{context}"
        );

        let mut targets = Vec::new();
        aligner.align_until_cost_limit(
            start,
            sequences.secondary_end(TsKind::TS21),
            U32Cost::from(5u8),
            &mut targets,
        );
        let target_costs: Vec<_> = targets
            .iter()
            .filter(|(coordinates, _)| *coordinates == target)
            .map(|(_, cost)| *cost)
            .collect();
        assert_eq!(
            target_costs.iter().min(),
            Some(&cost),
            "{context}: {targets:?}"
        );
        assert!(
            target_costs.contains(&U32Cost::from(5u8)),
            "{context}: {target_costs:?}"
        );
    }
}

/// The exact costs reported by `align_until_cost_limit` (via `fill_additional_targets`) must agree with
/// the costs of `align` to the same coordinates, otherwise the consistency checks of the chaining cost function fail.
/// This covers empty and non-empty descendant gaps, starts at the start of the alignment and inner starts,
/// and exact and inexact anchors.
#[test]
fn additional_targets_agree_with_align() {
    let sequences = matches_then_template_switch();
    let costs = chainalign_costs();
    let ts_kind = TsKind::TS21;
    let cost_limit = U32Cost::from(8u8);
    let (len1, len2) = (16, 16);

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        for without_enforced_non_match_next_to_anchors in [false, true] {
            let new_aligner = || {
                let aligner = Ts12JumpAligner::new(
                    &sequences,
                    &costs,
                    &rc_fn,
                    max_match_run,
                    enforcement_offset,
                );
                if without_enforced_non_match_next_to_anchors {
                    aligner.without_enforced_non_match_next_to_anchors()
                } else {
                    aligner
                }
            };

            for start in [
                PrimaryAlignmentCoordinates::new(0, 0),
                PrimaryAlignmentCoordinates::new(1, 1),
                PrimaryAlignmentCoordinates::new(2, 2),
                PrimaryAlignmentCoordinates::new(3, 2),
            ] {
                let mut targets = Vec::new();
                new_aligner().align_until_cost_limit(
                    start,
                    sequences.secondary_end(ts_kind),
                    cost_limit,
                    &mut targets,
                );

                for ancestor in 0..=len2 {
                    for descendant in start.a()..=len1 {
                        let end = SpecificSecondaryAlignmentCoordinates::new(
                            ancestor, descendant, ts_kind,
                        );
                        let reported_cost = targets
                            .iter()
                            .filter(|(coordinates, _)| *coordinates == end.into_any())
                            .map(|(_, cost)| *cost)
                            .min();
                        let (cost, alignment) = new_aligner().align(start, end, &mut Vec::new());
                        let expected_cost = (cost <= cost_limit).then_some(cost);
                        assert_eq!(
                            reported_cost, expected_cost,
                            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}, \
                             without_enforced_non_match_next_to_anchors = {without_enforced_non_match_next_to_anchors}, \
                             from {start} to {end}: {alignment:?}"
                        );
                    }
                }
            }
        }
    }
}
