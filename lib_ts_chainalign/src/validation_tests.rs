//! The tschainalign API rejects invalid inputs with an error instead of panicking (issue #249).

use generic_a_star::{AStarError, arena::ArenaFullError, cost::U32Cost};
use lib_tsalign::{
    alignment_geometry::{AlignmentRange, InvalidAlignmentRange},
    arena_index::ArenaIndexType,
    error::SequenceRole,
};

use crate::{
    ChainAlignError,
    chain_align::{
        limits::AlignmentLimits, performance_parameters::AlignmentPerformanceParameters,
    },
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind, InvalidChainingLowerBoundConfig},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    dna_rc_fn,
};

fn costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    )
}

fn config(anchor_k: u32, max_anchor_mutations: u8) -> ChainingLowerBoundConfig {
    ChainingLowerBoundConfig {
        anchor_k,
        max_anchor_mutations,
        inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
    }
}

fn lower_bounds(max_n: usize) -> ChainingLowerBounds<U32Cost> {
    crate::preprocess(max_n, config(4, 0), costs()).unwrap()
}

fn align(
    reference: &[u8],
    query: &[u8],
    range: AlignmentRange,
    lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> Result<lib_tsalign::alignment_result::AlignmentResult<U32Cost>, ChainAlignError> {
    crate::align(
        reference.to_vec(),
        query.to_vec(),
        range,
        &AlignmentPerformanceParameters::default(),
        &dna_rc_fn,
        "reference",
        "query",
        lower_bounds,
    )
}

#[test]
fn preprocess_rejects_invalid_config() {
    for (anchor_k, max_anchor_mutations) in [(0, 0), (3, 3), (65, 0), (256, 1)] {
        assert!(
            crate::preprocess(16, config(anchor_k, max_anchor_mutations), costs()).is_err(),
            "k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
        );
    }
    assert_eq!(
        crate::preprocess(16, config(3, 3), costs()).err(),
        Some(InvalidChainingLowerBoundConfig::TooManyAnchorMutations {
            anchor_k: 3,
            max_anchor_mutations: 3
        })
    );
}

#[test]
fn accepts_valid_input() {
    let lower_bounds = lower_bounds(16);
    align(
        b"ACGTACGTAC",
        b"ACGTTCGTAC",
        AlignmentRange::new_complete(10, 10),
        &lower_bounds,
    )
    .unwrap();
    align(
        b"ACGTACGTACGTACGT",
        b"ACGTTCGTAC",
        AlignmentRange::new(2, 2, 8, 8),
        &lower_bounds,
    )
    .unwrap();
}

#[test]
fn rejects_range_outside_of_sequences() {
    let lower_bounds = lower_bounds(16);
    for range in [
        AlignmentRange::new(0, 0, 11, 10),
        AlignmentRange::new(0, 0, 10, 11),
        AlignmentRange::new(12, 0, 14, 10),
    ] {
        let error = align(b"ACGTACGTAC", b"ACGTTCGTAC", range.clone(), &lower_bounds).unwrap_err();
        assert!(
            matches!(
                error,
                ChainAlignError::InvalidRange(InvalidAlignmentRange::OutOfBounds { .. })
            ),
            "{range}: {error:?}"
        );
    }
}

#[test]
fn rejects_range_with_offset_after_limit() {
    let lower_bounds = lower_bounds(16);
    let error = align(
        b"ACGTACGTAC",
        b"ACGTTCGTAC",
        AlignmentRange::new(5, 0, 4, 10),
        &lower_bounds,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            ChainAlignError::InvalidRange(InvalidAlignmentRange::OffsetAfterLimit { .. })
        ),
        "{error:?}"
    );
}

#[test]
fn rejects_sequences_longer_than_preprocessed() {
    let lower_bounds = lower_bounds(8);
    let error = align(
        b"ACGTACGTAC",
        b"ACGTTCGT",
        AlignmentRange::new_complete(10, 8),
        &lower_bounds,
    )
    .unwrap_err();
    assert_eq!(
        error,
        ChainAlignError::SequenceTooLong {
            sequence: SequenceRole::Reference,
            length: 10,
            max_n: 8
        }
    );

    let error = align(
        b"ACGTACGT",
        b"ACGTTCGTAC",
        AlignmentRange::new_complete(8, 10),
        &lower_bounds,
    )
    .unwrap_err();
    assert_eq!(
        error,
        ChainAlignError::SequenceTooLong {
            sequence: SequenceRole::Query,
            length: 10,
            max_n: 8
        }
    );
}

#[test]
fn rejects_non_dna_characters() {
    let lower_bounds = lower_bounds(16);
    let error = align(
        b"ACGTNCGTAC",
        b"ACGTTCGTAC",
        AlignmentRange::new_complete(10, 10),
        &lower_bounds,
    )
    .unwrap_err();
    assert_eq!(
        error,
        ChainAlignError::NonDnaCharacter {
            sequence: SequenceRole::Reference,
            position: 4,
            character: 'N'
        }
    );

    let error = align(
        b"ACGTACGTAC",
        b"ACGTTCGTAa",
        AlignmentRange::new_complete(10, 10),
        &lower_bounds,
    )
    .unwrap_err();
    assert_eq!(
        error,
        ChainAlignError::NonDnaCharacter {
            sequence: SequenceRole::Query,
            position: 9,
            character: 'a'
        }
    );
}

/// A pseudo-random reference of the given length and a query with every 50th character substituted.
fn arena_test_sequences(length: usize) -> (Vec<u8>, Vec<u8>) {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let reference: Vec<u8> = (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect();
    let query = reference
        .iter()
        .enumerate()
        .map(|(index, &character)| {
            if index % 50 == 25 {
                match character {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                }
            } else {
                character
            }
        })
        .collect();
    (reference, query)
}

/// If the chaining search needs more nodes at the same time than the arena index type allows,
/// then the alignment fails with an error instead of panicking.
/// The same alignment succeeds with the default index type `u32`, and gives the same result with `u64`.
#[test]
fn full_arena_is_an_error() {
    // The chaining of these sequences needs more than 65535 nodes at the same time (but not 64 characters).
    let (reference, query) = arena_test_sequences(128);
    let max_n = reference.len().max(query.len()).next_power_of_two();
    let lower_bounds = lower_bounds(max_n);
    let range = AlignmentRange::new_complete(reference.len(), query.len());
    // Generating all successors at once fills the open list quickly.
    let performance_parameters = || AlignmentPerformanceParameters {
        max_successors: usize::MAX,
        ..AlignmentPerformanceParameters::default()
    };
    let align = |arena_index| {
        crate::align(
            reference.clone(),
            query.clone(),
            range.clone(),
            &AlignmentPerformanceParameters {
                arena_index,
                ..performance_parameters()
            },
            &dna_rc_fn,
            "reference",
            "query",
            &lower_bounds,
        )
    };

    let error = crate::align_with_index::<u16>(
        reference.clone(),
        query.clone(),
        range.clone(),
        &performance_parameters(),
        &AlignmentLimits::none(),
        &dna_rc_fn,
        "reference",
        "query",
        &lower_bounds,
    )
    .unwrap_err();
    assert_eq!(
        error,
        ChainAlignError::Search(AStarError::ArenaFull(ArenaFullError {
            capacity: 65_535,
            index_type: "u16"
        }))
    );
    assert_eq!(
        error.to_string(),
        "the A* node arena is full: it can hold at most 65535 nodes with u16 node indices"
    );

    let u32 = align(ArenaIndexType::U32).unwrap();
    let u64 = align(ArenaIndexType::U64).unwrap();
    assert_eq!(u32.cigar(), u64.cigar());
    assert_eq!(u32.statistics().cost, u64.statistics().cost);
}
