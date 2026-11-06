use std::collections::BTreeMap;

use generic_a_star::cost::{AStarCost, U32Cost};
use num_traits::{Zero, bounds::UpperBounded};

use crate::alignment::coordinates::{
    AlignmentCoordinates, PrimaryAlignmentCoordinates, SpecificSecondaryAlignmentCoordinates,
};
use crate::alignment::ts_kind::TsKind;
use crate::alignment::{Alignment, AlignmentType};
use crate::alignment_history::history_vec::UnsignedIntAlignmentHistoryVec;
use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation,
    test_util::{Rng, alignment_cost, contains_anchor, ends_with_anchor},
};
use crate::panic_on_extend::PanicOnExtend;
use crate::windowed_min_mutation_chaining::gap_affine::GapAffineAligner;
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
fn start_end() {
    let seq1 = b"ACGT".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(4, 5);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(
        alignment.alignment,
        vec![(4, AlignmentType::Match), (1, AlignmentType::GapA)]
    );
    assert_eq!(cost, U32Cost::from(3u8));
}

#[test]
fn partial_alignment() {
    let seq1 = b"ACCGT".to_vec();
    let seq2 = b"ACGGTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(4, 4);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
fn gap_directions() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(11, 11);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
fn extremity_gaps() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(3, 3);
    let end = AlignmentCoordinates::new_primary(10, 10);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
fn extremity_substitutions() {
    let seq1 = b"AGGGA".to_vec();
    let seq2 = b"TGGGT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(5, 5);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
fn substitutions_as_gaps() {
    let seq1 = b"AAAAA".to_vec();
    let seq2 = b"TTTTT".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(3u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(5, 5);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(
        alignment.alignment == vec![(5, AlignmentType::GapA), (5, AlignmentType::GapB),]
            || alignment.alignment == vec![(5, AlignmentType::GapB), (5, AlignmentType::GapA),]
    );
    assert_eq!(cost, U32Cost::from(14u8));
}

#[test]
fn max_match_run_0() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        1,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(
        alignment.alignment == vec![(8, AlignmentType::GapA), (8, AlignmentType::GapB),]
            || alignment.alignment == vec![(8, AlignmentType::GapB), (8, AlignmentType::GapA),]
    );
    assert_eq!(cost, U32Cost::from(20u8));
}

#[test]
fn max_match_run_1() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        2,
        0,
    );
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
fn max_match_run_2() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
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
fn secondary_12() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS12);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS12);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(1, AlignmentType::Substitution), (7, AlignmentType::Match)]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn secondary_21() {
    let seq1 = b"GAAAAAAATG".to_vec();
    let seq2 = b"GTTTTTTTTG".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS21);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS21);
    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(
        alignment.alignment,
        vec![(7, AlignmentType::Match), (1, AlignmentType::Substitution)]
    );
    assert_eq!(cost, U32Cost::from(2u8));
}

#[test]
fn start_end_direct() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
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
fn start_anchor_direct() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_start(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
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
fn anchor_end_direct() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
    let (cost, alignment) = aligner.align(
        sequences.primary_end(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        5,
        0,
    );
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
fn start_end_indirect_lt_k() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        9,
        0,
    );
    let (cost, alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert_eq!(alignment.alignment, vec![(8, AlignmentType::Match)]);
    assert_eq!(cost, 0u8.into());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        9,
        0,
    );
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
fn start_end_indirect_geq_k() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        8,
        0,
    );
    let (cost, _alignment) = aligner.align(
        sequences.primary_start(),
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(!cost.is_zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        8,
        0,
    );
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
fn start_anchor_indirect() {
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

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, _alignment) = aligner.align(
        sequences.primary_start(),
        end,
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(cost.is_zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert!(min_cost.is_zero());
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
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(8, 8);

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, _alignment) = aligner.align(
        start,
        sequences.primary_end(),
        &mut Vec::new(),
        &mut PanicOnExtend,
    );

    assert!(cost.is_zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert!(min_cost.is_zero());
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
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(5, 5);
    let end = PrimaryAlignmentCoordinates::new(5, 5);

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, U32Cost::zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert_eq!(min_cost, U32Cost::zero());
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
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = PrimaryAlignmentCoordinates::new(5, 5);
    let end = PrimaryAlignmentCoordinates::new(6, 6);

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, _alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert!(cost.is_zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert!(min_cost.is_zero());
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
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);
    let end = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert_eq!(alignment.alignment, vec![]);
    assert_eq!(cost, U32Cost::zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert_eq!(min_cost, U32Cost::zero());
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
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = SpecificSecondaryAlignmentCoordinates::new(6, 4, TsKind::TS12);
    let end = SpecificSecondaryAlignmentCoordinates::new(5, 5, TsKind::TS12);

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
    let (cost, _alignment) = aligner.align(start, end, &mut PanicOnExtend, &mut Vec::new());

    assert!(cost.is_zero());

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u128>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        100,
        0,
    );
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
    assert!(min_cost.is_zero());
}

#[test]
fn start_end_3_0() {
    let seq1 = b"AAAAAAAA".to_vec();
    let seq2 = b"AAAAAAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = sequences.primary_start();
    let end = sequences.primary_end();

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u8>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
    let (cost, _alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(cost, U32Cost::from_primitive(8));

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u8>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        0,
    );
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
    assert_eq!(min_cost, U32Cost::from_primitive(8));
}

#[test]
fn start_end_3_1() {
    let seq1 = b"AAAAAAAA".to_vec();
    let seq2 = b"AAAAAAAA".to_vec();
    let sequences = AlignmentSequences::new_complete(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));
    let start = sequences.primary_start();
    let end = sequences.primary_end();

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u8>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        1,
    );
    let (cost, _alignment) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);

    assert_eq!(cost, U32Cost::from_primitive(14));

    let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u8>>::new(
        &sequences,
        &cost_table,
        &rc_fn,
        3,
        1,
    );
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
    assert_eq!(min_cost, U32Cost::from_primitive(14));
}

const PARETO_TEST_COSTS: [(u8, u8, u8); 5] =
    [(2, 3, 1), (1, 1, 1), (3, 2, 1), (1, 4, 2), (2, 5, 0)];

fn random_cost_table(rng: &mut Rng) -> GapAffineCosts<U32Cost> {
    let (substitution, gap_open, gap_extend) =
        PARETO_TEST_COSTS[rng.below(PARETO_TEST_COSTS.len())];
    GapAffineCosts::new(
        U32Cost::from(substitution),
        U32Cost::from(gap_open),
        U32Cost::from(gap_extend),
    )
}

fn primitive_cost_table(cost_table: &GapAffineCosts<U32Cost>) -> GapAffineCosts<u32> {
    GapAffineCosts::new(
        cost_table.substitution.as_primitive(),
        cost_table.gap_open.as_primitive(),
        cost_table.gap_extend.as_primitive(),
    )
}

/// Converts a gap-affine alignment into alignment history operations.
fn alignment_operations(alignment: &Alignment) -> Vec<AlignmentHistoryOperation> {
    alignment
        .iter_unpacked()
        .map(|alignment_type| match alignment_type {
            AlignmentType::Match => AlignmentHistoryOperation::Match,
            AlignmentType::Substitution => AlignmentHistoryOperation::Substitution,
            AlignmentType::GapA => AlignmentHistoryOperation::GapInA,
            AlignmentType::GapB => AlignmentHistoryOperation::GapInB,
            alignment_type => panic!("Unexpected alignment type {alignment_type:?}"),
        })
        .collect()
}

/// Checks that the alignment is a valid primary alignment from `start` to `end` with the given cost.
#[expect(clippy::too_many_arguments)]
fn assert_valid_primary_alignment(
    alignment: &Alignment,
    (seq1, seq2): (&[u8], &[u8]),
    start: (usize, usize),
    end: (usize, usize),
    anchor_k: u8,
    max_anchor_mutations: u8,
    cost_table: &GapAffineCosts<U32Cost>,
    cost: U32Cost,
) {
    let operations = alignment_operations(alignment);
    assert!(!contains_anchor(
        &operations,
        anchor_k.into(),
        max_anchor_mutations.into()
    ));
    assert_eq!(
        alignment_cost(&operations, &primitive_cost_table(cost_table)),
        cost.as_primitive()
    );

    let (mut a, mut b) = start;
    for operation in operations {
        match operation {
            AlignmentHistoryOperation::Match => assert_eq!(seq1[a], seq2[b]),
            AlignmentHistoryOperation::Substitution => assert_ne!(seq1[a], seq2[b]),
            _ => {}
        }
        a += operation.len_a();
        b += operation.len_b();
    }
    assert_eq!((a, b), end);
}

/// Computes the cost of an optimal alignment without anchors by enumerating all alignments.
fn brute_force_primary_cost(
    (seq1, seq2): (&[u8], &[u8]),
    start: (usize, usize),
    end: (usize, usize),
    anchor_k: usize,
    max_anchor_mutations: usize,
    cost_table: &GapAffineCosts<u32>,
) -> Option<u32> {
    #[expect(clippy::too_many_arguments)]
    fn recurse(
        alignment: &mut Vec<AlignmentHistoryOperation>,
        (a, b): (usize, usize),
        (seq1, seq2): (&[u8], &[u8]),
        end: (usize, usize),
        anchor_k: usize,
        max_anchor_mutations: usize,
        cost_table: &GapAffineCosts<u32>,
        result: &mut Option<u32>,
    ) {
        if (a, b) == end {
            let cost = alignment_cost(alignment, cost_table);
            *result = Some(result.map_or(cost, |result| result.min(cost)));
            return;
        }

        let mut operations = Vec::new();
        if a < end.0 && b < end.1 {
            operations.push(if seq1[a] == seq2[b] {
                AlignmentHistoryOperation::Match
            } else {
                AlignmentHistoryOperation::Substitution
            });
        }
        if a < end.0 {
            operations.push(AlignmentHistoryOperation::GapInB);
        }
        if b < end.1 {
            operations.push(AlignmentHistoryOperation::GapInA);
        }

        for operation in operations {
            alignment.push(operation);
            if !ends_with_anchor(alignment, anchor_k, max_anchor_mutations) {
                recurse(
                    alignment,
                    (a + operation.len_a(), b + operation.len_b()),
                    (seq1, seq2),
                    end,
                    anchor_k,
                    max_anchor_mutations,
                    cost_table,
                    result,
                );
            }
            alignment.pop();
        }
    }

    let mut result = None;
    recurse(
        &mut Vec::new(),
        start,
        (seq1, seq2),
        end,
        anchor_k,
        max_anchor_mutations,
        cost_table,
        &mut result,
    );
    result
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

#[test]
fn pareto_alignment_agrees_with_brute_force() {
    let mut rng = Rng::new(0);

    for _ in 0..400 {
        let len = rng.below(7);
        let seq1 = rng.dna(len);
        let mut seq2 = rng.mutate(&seq1, 30);
        seq2.truncate(7);
        let anchor_k = 2 + rng.below(3) as u8;
        let max_anchor_mutations = rng.below(usize::from(anchor_k).min(3)) as u8;
        let cost_table = random_cost_table(&mut rng);
        let start = (rng.below(seq1.len() / 2 + 1), rng.below(seq2.len() / 2 + 1));
        let end = (
            start.0 + rng.below(seq1.len() - start.0 + 1),
            start.1 + rng.below(seq2.len() - start.1 + 1),
        );
        let description = format!(
            "seq1 = {}, seq2 = {}, start = {start:?}, end = {end:?}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}",
            String::from_utf8_lossy(&seq1),
            String::from_utf8_lossy(&seq2),
        );

        let expected = brute_force_primary_cost(
            (&seq1, &seq2),
            start,
            end,
            anchor_k.into(),
            max_anchor_mutations.into(),
            &primitive_cost_table(&cost_table),
        );

        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
            &sequences,
            &cost_table,
            &rc_fn,
            anchor_k,
            max_anchor_mutations,
        );
        let (cost, alignment) = aligner.align(
            AlignmentCoordinates::new_primary(start.0, start.1),
            AlignmentCoordinates::new_primary(end.0, end.1),
            &mut Vec::new(),
            &mut PanicOnExtend,
        );

        if let Some(expected) = expected {
            assert_eq!(cost.as_primitive(), expected, "{description}");
            assert_valid_primary_alignment(
                &alignment,
                (&seq1, &seq2),
                start,
                end,
                anchor_k,
                max_anchor_mutations,
                &cost_table,
                cost,
            );
        } else {
            assert_eq!(cost, U32Cost::max_value(), "{description}");
        }
    }
}

#[test]
fn pareto_pruning_does_not_change_primary_alignments() {
    let mut rng = Rng::new(1);

    for _ in 0..150 {
        let len = 5 + rng.below(30);
        let seq1 = rng.dna(len);
        let mutation_percent = 5 + rng.below(25);
        let seq2 = rng.mutate(&seq1, mutation_percent);
        let anchor_k = [3, 4, 5, 6, 8][rng.below(5)];
        let max_anchor_mutations = rng.below(usize::from(anchor_k).min(4)) as u8;
        let cost_table = random_cost_table(&mut rng);
        let start = (rng.below(seq1.len() / 3 + 1), rng.below(seq2.len() / 3 + 1));
        let end = (
            seq1.len() - rng.below(seq1.len() / 3 + 1),
            seq2.len() - rng.below(seq2.len() / 3 + 1),
        );
        let cost_limit = U32Cost::from(rng.below(12) as u32);
        let description = format!(
            "seq1 = {}, seq2 = {}, start = {start:?}, end = {end:?}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}",
            String::from_utf8_lossy(&seq1),
            String::from_utf8_lossy(&seq2),
        );

        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let mut aligners = [true, false].map(|pareto_pruning| {
            let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
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
            let (cost, alignment) = aligner.align(
                AlignmentCoordinates::new_primary(start.0, start.1),
                AlignmentCoordinates::new_primary(end.0, end.1),
                &mut targets,
                &mut PanicOnExtend,
            );
            (cost, alignment, targets)
        });

        assert_eq!(cost, unpruned_cost, "{description}");
        if cost != U32Cost::max_value() {
            assert_valid_primary_alignment(
                &alignment,
                (&seq1, &seq2),
                start,
                end,
                anchor_k,
                max_anchor_mutations,
                &cost_table,
                cost,
            );
            assert_eq!(
                min_costs_up_to(&targets, cost),
                min_costs_up_to(&unpruned_targets, cost),
                "{description}"
            );
        }

        // Reuse the aligners to check that no state leaks between alignments.
        let [targets, unpruned_targets] = aligners.each_mut().map(|aligner| {
            let mut targets = Vec::new();
            aligner.align_until_cost_limit(
                AlignmentCoordinates::new_primary(start.0, start.1),
                cost_limit,
                &mut targets,
                &mut PanicOnExtend,
            );
            targets
        });
        assert_eq!(
            min_costs_up_to(&targets, cost_limit),
            min_costs_up_to(&unpruned_targets, cost_limit),
            "{description}, cost_limit = {cost_limit}"
        );
    }
}

#[test]
fn pareto_pruning_does_not_change_secondary_alignments() {
    let mut rng = Rng::new(2);

    for _ in 0..100 {
        let len = 5 + rng.below(25);
        let seq1 = rng.dna(len);
        let mutation_percent = 5 + rng.below(25);
        let seq2 = rng.mutate(&seq1, mutation_percent);
        let anchor_k = [3, 4, 5, 6][rng.below(4)];
        let max_anchor_mutations = rng.below(usize::from(anchor_k).min(3)) as u8;
        let cost_table = random_cost_table(&mut rng);
        let ts_kind = TsKind::iter().nth(rng.below(4)).unwrap();
        let (ancestor_len, descendant_len) = match ts_kind {
            TsKind::TS11 => (seq1.len(), seq1.len()),
            TsKind::TS12 => (seq1.len(), seq2.len()),
            TsKind::TS21 => (seq2.len(), seq1.len()),
            TsKind::TS22 => (seq2.len(), seq2.len()),
        };
        let start = (
            ancestor_len - rng.below(ancestor_len / 3 + 1),
            rng.below(descendant_len / 3 + 1),
        );
        let end = (
            rng.below(ancestor_len / 3 + 1),
            descendant_len - rng.below(descendant_len / 3 + 1),
        );
        let description = format!(
            "seq1 = {}, seq2 = {}, ts_kind = {ts_kind}, start = {start:?}, end = {end:?}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}",
            String::from_utf8_lossy(&seq1),
            String::from_utf8_lossy(&seq2),
        );

        let sequences = AlignmentSequences::new_complete(seq1.clone(), seq2.clone());
        let [
            (cost, alignment, targets),
            (unpruned_cost, _, unpruned_targets),
        ] = [true, false].map(|pareto_pruning| {
            let mut aligner = GapAffineAligner::<_, UnsignedIntAlignmentHistoryVec<u64>>::new(
                &sequences,
                &cost_table,
                &rc_fn,
                anchor_k,
                max_anchor_mutations,
            );
            aligner.set_pareto_pruning(pareto_pruning);
            let mut targets = Vec::new();
            let (cost, alignment) = aligner.align(
                AlignmentCoordinates::new_secondary(start.0, start.1, ts_kind),
                AlignmentCoordinates::new_secondary(end.0, end.1, ts_kind),
                &mut PanicOnExtend,
                &mut targets,
            );
            (cost, alignment, targets)
        });

        assert_eq!(cost, unpruned_cost, "{description}");
        if cost != U32Cost::max_value() {
            let operations = alignment_operations(&alignment);
            assert!(
                !contains_anchor(&operations, anchor_k.into(), max_anchor_mutations.into()),
                "{description}"
            );
            assert_eq!(
                min_costs_up_to(&targets, cost),
                min_costs_up_to(&unpruned_targets, cost),
                "{description}"
            );
        }
    }
}
