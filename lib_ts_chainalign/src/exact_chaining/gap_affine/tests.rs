use generic_a_star::cost::U32Cost;
use num_traits::{Zero, bounds::UpperBounded};

use crate::alignment::AlignmentType;
use crate::alignment::coordinates::{
    PrimaryAlignmentCoordinates, SpecificSecondaryAlignmentCoordinates,
};
use crate::alignment::ts_kind::TsKind;
use crate::exact_chaining::gap_affine::{AlignmentCoordinates, GapAffineAligner};
use crate::panic_on_extend::PanicOnExtend;
use crate::{alignment::sequences::AlignmentSequences, costs::GapAffineCosts};

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
    let seq1 = b"ACGT".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(4, 5);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![(4, AlignmentType::Match), (1, AlignmentType::GapA)]
    );
    assert_eq!(cost, U32Cost::from(3u8));
}

#[test]
fn test_partial_alignment() {
    let seq1 = b"ACCGT".to_vec();
    let seq2 = b"ACGGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(4, 4);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match)
        ]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn test_gap_directions() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(11, 11);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::Match),
            (2, AlignmentType::GapB),
            (5, AlignmentType::Match),
            (2, AlignmentType::GapA),
            (1, AlignmentType::Match)
        ]
    );
    assert_eq!(cost, U32Cost::from(8u8));
}

#[test]
fn test_extremity_gaps() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(3, 3);
    let end = AlignmentCoordinates::new_primary(10, 10);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (2, AlignmentType::GapB),
            (5, AlignmentType::Match),
            (2, AlignmentType::GapA),
        ]
    );
    assert_eq!(cost, U32Cost::from(8u8));
}

#[test]
fn test_extremity_substitutions() {
    let seq1 = b"AGGGA".to_vec();
    let seq2 = b"TGGGT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(5, 5);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Substitution),
            (3, AlignmentType::Match),
            (1, AlignmentType::Substitution),
        ]
    );
    assert_eq!(cost, U32Cost::from(4u8));
}

#[test]
fn test_substitutions_as_gaps() {
    let seq1 = b"AAAAAAAAAAAAAAAAAAAA".to_vec();
    let seq2 = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(3u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(20, 20);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(
        alignment.alignment == vec![(20, AlignmentType::GapA), (20, AlignmentType::GapB),]
            || alignment.alignment == vec![(20, AlignmentType::GapB), (20, AlignmentType::GapA),]
    );
    assert_eq!(cost, U32Cost::from(44u8));
}

#[test]
fn test_max_match_run_0() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 0);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(
        alignment.alignment == vec![(8, AlignmentType::GapA), (8, AlignmentType::GapB),]
            || alignment.alignment == vec![(8, AlignmentType::GapB), (8, AlignmentType::GapA),]
    );
    assert_eq!(cost, U32Cost::from(20u8));
}

#[test]
fn test_max_match_run_1() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 1);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (1, AlignmentType::Match),
            (2, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::GapA),
        ]
    );
    assert_eq!(cost, U32Cost::from(12u8));
}

#[test]
fn test_max_match_run_2() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 2);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (2, AlignmentType::Substitution),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(cost, U32Cost::from(6u8));
}

#[test]
fn test_secondary_12() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS12);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS12);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::Substitution), (7, AlignmentType::Match)]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn test_secondary_21() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS21);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS21);
    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(7, AlignmentType::Match), (1, AlignmentType::Substitution)]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn test_start_end_direct() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(5, 5),
        PrimaryAlignmentCoordinates::new(5, 5),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_start(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_start() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert_eq!(min_cost, 0u8.into());
}

#[test]
fn test_start_anchor_direct() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(5, 5),
        PrimaryAlignmentCoordinates::new(6, 6),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_start(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_start(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_start() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert_eq!(min_cost, 0u8.into());
}

#[test]
fn test_anchor_end_direct() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(5, 5),
        PrimaryAlignmentCoordinates::new(6, 6),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let (cost, alignment) = aligner.align(
        sequences.primary_end(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 4);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_end(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_end() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert_eq!(min_cost, 0u8.into());
}

#[test]
fn test_start_end_indirect_lt_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 8);
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![(8, AlignmentType::Match)]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 8);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_start(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_end() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert_eq!(min_cost, 0u8.into());
}

#[test]
fn test_start_end_indirect_geq_k() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 7);
    let (cost, _alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, 7);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_start(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_end() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert!(!min_cost.is_zero());
}

#[test]
fn test_start_anchor_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let end = PrimaryAlignmentCoordinates::new(2, 2);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, _alignment) = aligner.align(
        sequences.primary_start(),
        end,
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        sequences.primary_start(),
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(
            |(target, cost)| {
                if target == end { Some(cost) } else { None }
            },
        )
        .min()
        .unwrap();
    assert!(!min_cost.is_zero());
}

#[test]
fn test_anchor_end_indirect() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(8, 8);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, _alignment) = aligner.align(
        start,
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        start,
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == sequences.primary_end() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert!(!min_cost.is_zero());
}

#[test]
fn test_anchor_anchor_direct_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(5, 5);
    let end = PrimaryAlignmentCoordinates::new(5, 5);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, U32Cost::max_value());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        start,
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(
            |(target, cost)| {
                if target == end { Some(cost) } else { None }
            },
        )
        .min()
        .unwrap_or_else(U32Cost::max_value);
    assert_eq!(min_cost, U32Cost::max_value());
}

#[test]
fn test_anchor_anchor_indirect_primary() {
    let seq1 = b"ATTTTTTTTA".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(5, 5);
    let end = PrimaryAlignmentCoordinates::new(6, 6);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, _alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_primary_targets = Vec::new();
    aligner.align_until_cost_limit(
        start,
        U32Cost::max_value(),
        &mut additional_primary_targets,
        &mut PanicOnExtend,
    );

    let min_cost = additional_primary_targets
        .into_iter()
        .filter_map(
            |(target, cost)| {
                if target == end { Some(cost) } else { None }
            },
        )
        .min()
        .unwrap();
    assert!(!min_cost.is_zero());
}

#[test]
fn test_anchor_anchor_direct_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);
    let end = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, U32Cost::max_value());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_secondary_targets = Vec::new();
    aligner.align_until_cost_limit(
        start,
        U32Cost::max_value(),
        &mut PanicOnExtend,
        &mut additional_secondary_targets,
    );

    let min_cost = additional_secondary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == end.into() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap_or_else(U32Cost::max_value);
    assert_eq!(min_cost, U32Cost::max_value());
}

#[test]
fn test_anchor_anchor_indirect_secondary() {
    let seq1 = b"GAAAAAAAAG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new(
        seq1,
        seq2,
        PrimaryAlignmentCoordinates::new(1, 1),
        PrimaryAlignmentCoordinates::new(9, 9),
    );
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = SpecificSecondaryAlignmentCoordinates::new(6, 4, TsKind::TS12);
    let end = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let (cost, _alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::new(&sequences, &cost_table, &rc_fn, u32::MAX);
    let mut additional_secondary_targets = Vec::new();
    aligner.align_until_cost_limit(
        start,
        U32Cost::max_value(),
        &mut PanicOnExtend,
        &mut additional_secondary_targets,
    );

    let min_cost = additional_secondary_targets
        .into_iter()
        .filter_map(|(target, cost)| {
            if target == end.into() {
                Some(cost)
            } else {
                None
            }
        })
        .min()
        .unwrap();
    assert!(!min_cost.is_zero());
}
