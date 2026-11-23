use num_traits::real::Real;

use crate::a_star_aligner::alignment_geometry::AlignmentCoordinates;

use super::*;

/// This test case uses real data and used to panic.
#[test]
fn test_panic() {
    let mut aligner = Aligner::new();
    aligner.set_min_length_strategy(MinLengthStrategySelector::PreprocessFilter);
    let res = aligner.align(
            "reference",
            b"TACCGAGACTGCAGAAAGTGAAAGCTATACTAA",
            "query",
            b"TAACTTTTAATGCCAAATATTTTATCCAAATAGGAAATTGTTTTCCGGTAAAATTTAACAAAAGAACCAGTTTACCCCCTTCAATGATTTATTTTTCTTCTTAGATTGAACTCTCGGGTTAGATCTCATTTTAACTGAAATTTGGTAAAAAATCCATATTACGGTTCAAGCCTAACCGAGACTGCAGAAAGTGAAAGCTAAAAGCTAATTTTTTTTTTTTTTTTGTATTTCACACCTATCGCAATACATCCTGGACAACACTGTATATTGAAACATTTTTTGCCTACAGCAATGGGCCTATAATTTTTTCTCGGCATTAGCTCTACAATCCAATTCTATCCTGCTTCTTCTTGTAAACAGGGATAACTTTAACTAACATTCAGTTTGCTTGGGAAAGAACCGATTGATAATGTA",
            AlignmentRange::new_offset_limit(
                AlignmentCoordinates::new(27, 200),
                AlignmentCoordinates::new(33, 214)
            ).into(),
            &[],
            None,
            None,
            true,
            TSMUncertaintyRangeExtensionMode::None,
        );
    println!("{res:#?}");
    assert!(res.statistics().cost.is_sign_positive());
}

#[test]
fn test_193() {
    //        [-- 10 --][-- 11 ---]1
    let r = b"AAAAAAAAAACCCCCCCCCCCC";
    let q = b"AAAAAAAAAATTTTTTTTTTTC";

    let mut tsa = Aligner::<DnaAlphabetOrN>::default();

    let mut costs = TemplateSwitchConfig::default();
    costs.base_cost.qrr = 4u64.into();
    tsa.set_costs(costs.clone());
    let _tsa_res = tsa.align(
        "ref",
        r,
        "qry",
        q,
        None,
        &[],
        None,
        None,
        false,
        TSMUncertaintyRangeExtensionMode::None,
    );
}

/// Align with the `lower_bound` chaining strategy and check that the result is the same as with the `none` chaining strategy.
///
/// The chaining lower bound is only an optimisation, so it must not change the alignment.
fn assert_lower_bound_chaining_equals_no_chaining(
    reference: &[u8],
    query: &[u8],
    range: Option<AlignmentRange>,
) {
    let align = |chaining_strategy| {
        let mut aligner = Aligner::new();
        aligner.set_chaining_strategy(chaining_strategy);
        aligner.align(
            "reference",
            reference,
            "query",
            query,
            range.clone(),
            &[],
            None,
            None,
            false,
            TSMUncertaintyRangeExtensionMode::None,
        )
    };

    let expected = align(ChainingStrategySelector::None);
    let actual = align(ChainingStrategySelector::LowerBound);
    println!("none: {expected:#?}");
    println!("lower_bound: {actual:#?}");

    assert_eq!(expected.statistics().cost, actual.statistics().cost);
    assert_eq!(expected.cigar(), actual.cigar());
}

// Issue #245: the `lower_bound` chaining strategy used to panic when the reference was shorter than the chaining block size (20).
//
// The following tests only use references shorter than the block size.
// References of at least the block size precompute the chaining lower bound matrices,
// which is too expensive with the default config to run in a debug build.

#[test]
fn test_lower_bound_chaining_reference_shorter_than_block() {
    assert_lower_bound_chaining_equals_no_chaining(b"ACGTACGT", b"ACGAACGT", None);
}

#[test]
fn test_lower_bound_chaining_reference_one_shorter_than_block() {
    assert_lower_bound_chaining_equals_no_chaining(
        b"ACGTACGTTGCAACGTACG",
        b"ACGTACGATGCAACGTACGTTAGCA",
        None,
    );
}

#[test]
fn test_lower_bound_chaining_short_reference_with_range() {
    assert_lower_bound_chaining_equals_no_chaining(
        b"ACGTACGTTGCAACG",
        b"ACGTACCTTGCAACG",
        AlignmentRange::new_offset_limit(
            AlignmentCoordinates::new(4, 4),
            AlignmentCoordinates::new(10, 10),
        )
        .into(),
    );
}

#[test]
fn test_lower_bound_chaining_empty_query() {
    assert_lower_bound_chaining_equals_no_chaining(b"ACGTACGT", b"", None);
}

#[test]
fn test_lower_bound_chaining_empty_reference() {
    assert_lower_bound_chaining_equals_no_chaining(b"", b"ACGTACGT", None);
}

#[test]
fn test_lower_bound_chaining_empty_reference_and_query() {
    assert_lower_bound_chaining_equals_no_chaining(b"", b"", None);
}
