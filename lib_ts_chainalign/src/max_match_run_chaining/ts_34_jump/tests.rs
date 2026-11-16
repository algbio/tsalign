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

/// The reference is an inner sequence followed by `CC`, and the query is the reverse complement
/// of the inner sequence followed by `CC`.
/// The optimal alignment is a template switch that covers the inner sequence completely followed by `2=` (cost 2).
fn template_switch_then_matches() -> AlignmentSequences {
    AlignmentSequences::new_complete(b"GATTACAGGCTTACCC".to_vec(), b"GTAAGCCTGTAATCCC".to_vec())
}

/// A 34-jump from the last secondary anchor of an optimal template switch to the end of the alignment,
/// with anchors of length k = 4 (exact, and with one mutation).
///
/// The secondary anchor `SA((4, 0], [10, 14))` of kind 21 ends directly at the template switch exit:
/// the descendant (seq1) is at position 14 and the ancestor (seq2) was read in reverse until position 0.
/// The optimal 34-jump is the jump back to `(14, 14)` followed by `2=` (the `CC` at the end,
/// shorter than k so no primary anchor covers it), so it costs 0
/// (the template switch base cost is counted in the 12-jump).
///
/// The end of the alignment is no anchor, so this primary part may consist of matches only.
/// But `Ts34JumpAligner::align` enforces a non-match whenever the descendant gap is non-empty
/// (`enforce_non_match`, independently of `enforce_non_match_next_to_anchors`),
/// which forces a deletion into the primary part (jump + `1D1=`, cost 3).
#[test]
fn jump_then_matches_to_alignment_end() {
    let sequences = template_switch_then_matches();
    let costs = chainalign_costs();
    let start = SpecificSecondaryAlignmentCoordinates::new(0, 14, TsKind::TS21);
    let end = sequences.primary_end();

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligner = Ts34JumpAligner::new(
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
                U32Cost::from(0u8),
                vec![
                    (1, AlignmentType::TsEnd { jump: 14 }),
                    (2, AlignmentType::Match),
                ]
            ),
            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}"
        );
    }
}

/// Same as [`jump_then_matches_to_alignment_end`], but via `align_until_cost_limit`,
/// which the chaining cost function uses to compute the exact costs of 34-jumps to the end.
///
/// With cost limit 0 (`max_exact_cost_function_cost = 0`), the end of the alignment must be reported with cost 0.
/// Currently, `fill_additional_targets` only reports nodes with a non-match if the descendant gap is non-empty,
/// so the optimal jump + `2=` is not reported.
#[test]
fn jump_then_matches_to_alignment_end_until_cost_limit() {
    let sequences = template_switch_then_matches();
    let costs = chainalign_costs();
    let start = SpecificSecondaryAlignmentCoordinates::new(0, 14, TsKind::TS21);
    let end = sequences.primary_end();

    for (max_match_run, enforcement_offset) in K4_ALIGNER_PARAMETERS {
        let mut aligner = Ts34JumpAligner::new(
            &sequences,
            &costs,
            &rc_fn,
            max_match_run,
            enforcement_offset,
        );
        let mut targets = Vec::new();
        aligner.align_until_cost_limit(start, U32Cost::from(0u8), &mut targets);

        let end_cost = targets
            .iter()
            .filter(|(coordinates, _)| *coordinates == end)
            .map(|(_, cost)| *cost)
            .min();
        assert_eq!(
            end_cost,
            Some(U32Cost::from(0u8)),
            "max_match_run = {max_match_run}, enforcement_offset = {enforcement_offset}: {targets:?}"
        );
    }
}
