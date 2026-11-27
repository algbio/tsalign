//! The aligner API rejects invalid inputs with an error instead of panicking or silently ignoring them (issue #249).

use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::sequence::OwnedGenomeSequence,
};

use crate::{
    alignment_geometry::{AlignmentCoordinates, InvalidAlignmentRange},
    error::{AlignerInputError, SequenceRole},
};

use super::*;

fn align(
    aligner: &Aligner<DnaAlphabet>,
    reference: &[u8],
    query: &[u8],
    range: Option<AlignmentRange>,
    gap_characters: &[u8],
) -> Result<AlignmentResult<U64Cost>, AlignerInputError> {
    aligner.align(
        "reference",
        reference,
        "query",
        query,
        range,
        gap_characters,
        None,
        None,
        true,
        TSMUncertaintyRangeExtensionMode::None,
    )
}

#[test]
fn accepts_valid_input() {
    let aligner = Aligner::<DnaAlphabet>::default();
    let result = align(&aligner, b"ACGTACGT", b"ACGACGT", None, &[]).unwrap();
    assert!(matches!(result, AlignmentResult::WithTarget { .. }));

    // A range that covers the complete sequences, and an empty range at the end, are valid.
    align(
        &aligner,
        b"ACGT-ACGT",
        b"ACGACGT",
        Some(AlignmentRange::new_complete(8, 7)),
        b"-",
    )
    .unwrap();
    align(
        &aligner,
        b"ACGTACGT",
        b"ACGACGT",
        Some(AlignmentRange::new(8, 7, 8, 7)),
        &[],
    )
    .unwrap();
}

#[test]
fn rejects_non_alphabet_characters() {
    let aligner = Aligner::<DnaAlphabet>::default();

    let error = align(&aligner, b"ACGNACGT", b"ACGACGT", None, &[]).unwrap_err();
    assert!(
        matches!(
            error,
            AlignerInputError::NonAlphabetCharacter {
                sequence: SequenceRole::Reference,
                ..
            }
        ),
        "{error:?}"
    );
    assert!(
        error
            .to_string()
            .starts_with("Reference contains non-alphabet character")
    );

    let error = align(&aligner, b"ACGTACGT", b"ACG-ACGT", None, &[]).unwrap_err();
    assert!(
        matches!(
            error,
            AlignerInputError::NonAlphabetCharacter {
                sequence: SequenceRole::Query,
                ..
            }
        ),
        "{error:?}"
    );
    assert!(
        error
            .to_string()
            .starts_with("Query contains non-alphabet character")
    );

    // Gap characters are removed before checking the alphabet.
    align(&aligner, b"ACGTACGT", b"ACG-ACGT", None, b"-").unwrap();
}

#[test]
fn rejects_range_outside_of_sequences() {
    let aligner = Aligner::<DnaAlphabet>::default();

    for range in [
        AlignmentRange::new(0, 0, 9, 7),
        AlignmentRange::new(0, 0, 8, 8),
        AlignmentRange::new(9, 0, 9, 7),
        AlignmentRange::new(0, 100, 8, 200),
    ] {
        let error = align(&aligner, b"ACGTACGT", b"ACGACGT", Some(range.clone()), &[]).unwrap_err();
        assert!(
            matches!(
                error,
                AlignerInputError::InvalidRange(InvalidAlignmentRange::OutOfBounds { .. })
            ),
            "{range}: {error:?}"
        );
    }

    // The range refers to the sequences without gap characters.
    let error = align(
        &aligner,
        b"ACGT-ACGT",
        b"ACGACGT",
        Some(AlignmentRange::new_complete(9, 7)),
        b"-",
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            AlignerInputError::InvalidRange(InvalidAlignmentRange::OutOfBounds { .. })
        ),
        "{error:?}"
    );
}

#[test]
fn rejects_range_with_offset_after_limit() {
    let aligner = Aligner::<DnaAlphabet>::default();

    for range in [
        AlignmentRange::new(5, 0, 4, 7),
        AlignmentRange::new(0, 5, 8, 4),
    ] {
        let error = align(&aligner, b"ACGTACGT", b"ACGACGT", Some(range.clone()), &[]).unwrap_err();
        assert!(
            matches!(
                error,
                AlignerInputError::InvalidRange(InvalidAlignmentRange::OffsetAfterLimit { .. })
            ),
            "{range}: {error:?}"
        );
    }
}

/// The aligner does not support TSMs that start or end outside of the alignment range yet.
/// It used to panic in this case.
#[test]
fn rejects_unsupported_ts_14_out_of_range() {
    let mut aligner = Aligner::<DnaAlphabet>::default();
    aligner.set_ts_14_out_of_range_strategy(Ts14OutOfRangeStrategy::Allow);

    assert_eq!(
        aligner.validate_settings(),
        Err(AlignerInputError::Ts14OutOfRangeNotSupported)
    );
    for range in [
        AlignmentRange::new(2, 2, 6, 5),
        // An invalid range is reported as well, but the unsupported setting first.
        AlignmentRange::new(10, 0, 12, 7),
    ] {
        let error = align(&aligner, b"ACGTACGT", b"ACGACGT", Some(range), &[]).unwrap_err();
        assert_eq!(error, AlignerInputError::Ts14OutOfRangeNotSupported);
    }
}

/// Computing the additional TSM starts and ends used to panic for ranges outside of the sequences.
#[test]
fn additional_tsm_starts_and_ends_reject_invalid_ranges() {
    use crate::aligner::ts_14_out_of_range::AdditionalExplicitTSMStartsAndEndsError;

    for range in [
        AlignmentRange::new(10, 0, 12, 7),
        AlignmentRange::new(0, 8, 8, 9),
        AlignmentRange::new(0, 0, 9, 7),
    ] {
        let error =
            AdditionalExplicitTSMStartsAndEnds::new("AC-GTACGT", "ACGACGT", &range, &['-'], false)
                .err()
                .unwrap();
        assert!(
            matches!(
                error,
                AdditionalExplicitTSMStartsAndEndsError::InvalidRange(
                    InvalidAlignmentRange::OutOfBounds { .. }
                )
            ),
            "{range}: {error:?}"
        );
    }

    // An empty range at the end is valid.
    AdditionalExplicitTSMStartsAndEnds::new(
        "AC-GTACGT",
        "ACGACGT",
        &AlignmentRange::new(8, 7, 8, 7),
        &['-'],
        false,
    )
    .unwrap();

    let error = AdditionalExplicitTSMStartsAndEnds::new(
        "AC|GTA|CGT",
        "ACG|ACGT",
        &AlignmentRange::new(2, 3, 5, 7),
        &[],
        true,
    )
    .err()
    .unwrap();
    assert_eq!(
        error,
        AdditionalExplicitTSMStartsAndEndsError::MissingEmbeddedRangeDelimiter(SequenceRole::Query)
    );
}

#[test]
fn align_sequences_rejects_invalid_range() {
    let aligner = Aligner::<DnaAlphabet>::default();
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACGTACGT").unwrap();
    let query = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACGACGT").unwrap();

    let error = aligner
        .align_sequences(
            "reference",
            reference.as_genome_subsequence(),
            "query",
            query.as_genome_subsequence(),
            AlignmentRange::new(0, 0, 8, 8),
            AdditionalExplicitTSMStartsAndEnds::default(),
            None,
            None,
            true,
            TSMUncertaintyRangeExtensionMode::None,
        )
        .unwrap_err();
    assert!(
        matches!(
            error,
            AlignerInputError::InvalidRange(InvalidAlignmentRange::OutOfBounds { .. })
        ),
        "{error:?}"
    );
}

/// Additional TSM starts and ends are only used if TSMs may start and end outside of the range.
#[test]
fn align_sequences_rejects_ignored_additional_tsm_starts_and_ends() {
    let aligner = Aligner::<DnaAlphabet>::default();
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACGTACGT").unwrap();
    let query = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACGACGT").unwrap();
    let additional_tsm_starts_and_ends = || AdditionalExplicitTSMStartsAndEnds {
        explicit_tsm_starts: vec![AlignmentCoordinates::new(1, 1)],
        explicit_tsm_ends: Vec::new(),
    };
    let align_sequences = |aligner: &Aligner<DnaAlphabet>| {
        aligner.align_sequences(
            "reference",
            reference.as_genome_subsequence(),
            "query",
            query.as_genome_subsequence(),
            AlignmentRange::new(2, 2, 6, 5),
            additional_tsm_starts_and_ends(),
            None,
            None,
            true,
            TSMUncertaintyRangeExtensionMode::None,
        )
    };

    let error = align_sequences(&aligner).unwrap_err();
    assert_eq!(
        error,
        AlignerInputError::IgnoredAdditionalExplicitTSMStartsAndEnds
    );
}

/// A pseudo-random reference of the given length and a query with every 20th character substituted.
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
            if index % 20 == 10 {
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

/// If the search needs more nodes at the same time than the arena index type allows,
/// then the aligner returns an error instead of panicking.
/// The same alignment succeeds with the default index type `u32`, and gives the same result with `u64`.
#[test]
fn full_arena_is_an_error() {
    // The search for these sequences needs more than 65535 nodes at the same time (but not for 50 characters).
    let (reference, query) = arena_test_sequences(70);
    let align_with = |arena_index| {
        let mut aligner = Aligner::<DnaAlphabet>::default();
        aligner.set_arena_index(arena_index);
        align(&aligner, &reference, &query, None, &[])
    };

    let error = align_with(ArenaIndexType::U16).unwrap_err();
    assert_eq!(
        error,
        AlignerInputError::Search(generic_a_star::AStarError::ArenaFull(
            generic_a_star::arena::ArenaFullError {
                capacity: 65_535,
                index_type: "u16"
            }
        ))
    );
    assert_eq!(
        error.to_string(),
        "the A* node arena is full: it can hold at most 65535 nodes with u16 node indices"
    );

    let default = align(
        &Aligner::<DnaAlphabet>::default(),
        &reference,
        &query,
        None,
        &[],
    )
    .unwrap();
    let u32 = align_with(ArenaIndexType::U32).unwrap();
    let u64 = align_with(ArenaIndexType::U64).unwrap();
    assert_eq!(default.cigar(), u32.cigar());
    assert_eq!(default.cigar(), u64.cigar());
    assert_eq!(default.statistics().cost, u64.statistics().cost);
    assert_eq!(
        default.statistics().opened_nodes,
        u64.statistics().opened_nodes
    );
}
