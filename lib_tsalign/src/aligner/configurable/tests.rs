use generic_a_star::cost::AStarCost;
use num_traits::real::Real;

use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;

use crate::aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection};
use crate::alignment_geometry::AlignmentCoordinates;
use crate::alignment_result::{AlignmentType, alignment::Alignment};

use compact_genome::{
    implementation::vec_sequence::VectorGenome, interface::sequence::OwnedGenomeSequence,
};

use super::*;

/// This test case uses real data and used to panic.
#[test]
fn panic() {
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
        ).unwrap();
    println!("{res:#?}");
    assert!(res.statistics().cost.is_sign_positive());
}

#[test]
fn issue_193() {
    //        [-- 10 --][-- 11 ---]1
    let r = b"AAAAAAAAAACCCCCCCCCCCC";
    let q = b"AAAAAAAAAATTTTTTTTTTTC";

    let mut tsa = Aligner::<DnaAlphabetOrN>::default();

    let mut costs = TemplateSwitchConfig::default();
    costs.base_cost.qrr = 4u64.into();
    tsa.set_costs(costs.clone());
    let _tsa_res = tsa
        .align(
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
        )
        .unwrap();
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
        aligner
            .align(
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
            .unwrap()
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

/// The costs of `test_files/config/chainalign/config.tsa`.
/// Template switches have base cost 2 and an inner length of at least 10.
const TS_END_TEST_COSTS: &str = r#"# Limits

left_flank_length = 0
right_flank_length = 0

# Base Cost

rrf_cost = 31
rqf_cost = 30
qrf_cost = 30
qqf_cost = 31
rrr_cost = 2
rqr_cost = 2
qrr_cost = 2
qqr_cost = 2

# Jump Costs

RQQROffset
 -inf -300 301
  inf    0 inf

RRQQOffset
 -inf -300 301
  inf    0 inf

Length
   0 10 200
 inf  0 inf

LengthDifference
 -inf -100 -20 -10 11 21 101
  inf    0   0   0  0  0 inf

ForwardAntiDescendantGap
 -inf   1
    0 inf

ReverseAntiDescendantGap
 -inf
    0

# Primary Edit Costs

SubstitutionCostTable
  |  A  C  G  T
--+------------
A |  0  2  2  2
C |  2  0  2  2
G |  2  2  0  2
T |  2  2  2  0

GapOpenCostVector
 A C G T
 3 3 3 3

GapExtendCostVector
 A C G T
 1 1 1 1

# Secondary Forward Edit Costs

SubstitutionCostTable
  |  A  C  G  T
--+------------
A |  0  4  4  4
C |  4  0  4  4
G |  4  4  0  4
T |  4  4  4  0

GapOpenCostVector
 A C G T
 3 3 3 3

GapExtendCostVector
 A C G T
 2 2 2 2

# Secondary Reverse Edit Costs

SubstitutionCostTable
  |  A  C  G  T
--+------------
A |  0  4  4  4
C |  4  0  4  4
G |  4  4  0  4
T |  4  4  4  0

GapOpenCostVector
 A C G T
 3 3 3 3

GapExtendCostVector
 A C G T
 2 2 2 2

# Left Flank Edit Costs

SubstitutionCostTable
  |  A  C  G  T
--+------------
A |  0  3  3  3
C |  3  0  3  3
G |  3  3  0  3
T |  3  3  3  0

GapOpenCostVector
 A C G T
 4 4 4 4

GapExtendCostVector
 A C G T
 1 1 1 1

# Right Flank Edit Costs

SubstitutionCostTable
  |  A  C  G  T
--+------------
A |  0  3  3  3
C |  3  0  3  3
G |  3  3  0  3
T |  3  3  3  0

GapOpenCostVector
 A C G T
 4 4 4 4

GapExtendCostVector
 A C G T
 1 1 1 1"#;

/// Align with [`TS_END_TEST_COSTS`] and return the cost and the number of template switches.
fn align_ts_end_test(reference: &[u8], query: &[u8]) -> (f64, f64) {
    let mut aligner = Aligner::<DnaAlphabet>::default();
    aligner.set_costs_parse(TS_END_TEST_COSTS).unwrap();
    let result = aligner
        .align(
            "reference",
            reference,
            "query",
            query,
            None,
            &[],
            None,
            None,
            false,
            TSMUncertaintyRangeExtensionMode::None,
        )
        .unwrap();
    println!("{result}");
    let statistics = result.statistics();
    (
        statistics.cost.raw(),
        statistics.template_switch_amount.raw(),
    )
}

/// A template switch may end at the end of the anti-descendant sequence.
/// The query is the reverse complement of the reference, so a single template switch covers both sequences completely.
#[test]
fn template_switch_ends_at_sequence_end_full_reverse_complement() {
    let (cost, template_switches) = align_ts_end_test(
        b"GAGGGTTTACCACTTCCAGTCATGGTGTAG",
        b"CTACACCATGACTGGAAGTGGTAAACCCTC",
    );
    assert_eq!((cost, template_switches), (2.0, 1.0));
}

/// A template switch may end at the end of the anti-descendant sequence.
/// The last 10 characters of the query are the reverse complement of the last 10 characters of the reference.
#[test]
fn template_switch_ends_at_sequence_end_inversion() {
    let (cost, template_switches) = align_ts_end_test(
        b"TTTGGTACGCTTCACCAACGGCGACAAGAAATACGTGGCG",
        b"TTTGGTACGCTTCACCAACGGCGACAAGAACGCCACGTAT",
    );
    assert_eq!((cost, template_switches), (2.0, 1.0));
}

/// The bounds of the anti-descendant gap at the template switch exit depend on the entrance index of the anti-descendant,
/// which differs from the entrance index of the descendant if the template switch starts off the main diagonal.
///
/// The query has two inserted characters, and then the last 11 characters of the reference are replaced by the
/// reverse complement of their last 10 characters.
/// The optimal alignment is a gap of length two (cost 4) and a template switch (cost 2)
/// whose anti-descendant part is one character longer than its descendant part and ends at the end of the reference.
#[test]
fn template_switch_exit_off_diagonal() {
    let (cost, template_switches) = align_ts_end_test(
        b"GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTT",
        b"GCTAAAGACAGGATTACATAACATACACGTCAAAGTTTCGTG",
    );
    assert_eq!((cost, template_switches), (6.0, 1.0));
}

/// Costs for testing the offset costs of forward template switches (issue #166).
///
/// Only forward template switches of the kinds selected by `qqf_cost` and `qrf_cost` are possible,
/// all other base costs are infinite.
/// The offset cost functions are given as two lines each (breakpoints and values).
/// The inner of a template switch has a length of at least 10 and costs zero if it matches exactly,
/// the anti-descendant must not advance during a forward template switch (at no cost),
/// and length differences up to 100 cost zero.
/// The edit costs are those of [`TS_END_TEST_COSTS`], i.e. a primary gap of length `l` costs `3 + (l - 1)`.
fn forward_ts_offset_test_costs(
    qqf_cost: &str,
    qrf_cost: &str,
    rq_qr_offset_costs: &str,
    rr_qq_offset_costs: &str,
) -> String {
    let (_, edit_costs) = TS_END_TEST_COSTS
        .split_once("# Primary Edit Costs")
        .unwrap();

    format!(
        r#"# Limits

left_flank_length = 0
right_flank_length = 0

# Base Cost

rrf_cost = inf
rqf_cost = inf
qrf_cost = {qrf_cost}
qqf_cost = {qqf_cost}
rrr_cost = inf
rqr_cost = inf
qrr_cost = inf
qqr_cost = inf

# Jump Costs

RQQROffset
{rq_qr_offset_costs}

RRQQOffset
{rr_qq_offset_costs}

Length
   0 10 200
 inf  0 inf

LengthDifference
 -inf -100 101
  inf    0 inf

ForwardAntiDescendantGap
 -inf   1
    0 inf

ReverseAntiDescendantGap
 -inf
    0

# Primary Edit Costs{edit_costs}"#
    )
}

/// Offset costs that are zero for offset zero, seven for offsets in `-20..0` and `1..21`, and infinite otherwise.
const FINITE_FORWARD_OFFSET_COSTS: &str = " -inf -20 0 1 21\n  inf   7 0 7 inf";
/// Offset costs that are zero for offset zero and infinite otherwise.
/// Since a forward template switch never has offset zero, this forbids forward template switches.
const ZERO_ONLY_FORWARD_OFFSET_COSTS: &str = " -inf 0 1\n  inf 0 inf";
/// Offset costs that are infinite everywhere.
const INFINITE_OFFSET_COSTS: &str = " -inf\n  inf";

/// Align with the given costs and return the cost, the number of template switches,
/// and the cost of the computed alignment recomputed independently of the A* search with
/// [`Alignment::compute_cost`](crate::alignment_result::alignment::Alignment::compute_cost).
fn align_forward_ts_offset_test(costs: &str, reference: &[u8], query: &[u8]) -> (f64, f64, f64) {
    align_and_recompute_cost(costs, Default::default(), reference, query)
}

/// Align with the given costs and min length strategy and return the cost, the number of template switches,
/// and the cost of the computed alignment recomputed independently of the A* search with
/// [`Alignment::compute_cost`](crate::alignment_result::alignment::Alignment::compute_cost).
fn align_and_recompute_cost(
    costs: &str,
    min_length_strategy: MinLengthStrategySelector,
    reference: &[u8],
    query: &[u8],
) -> (f64, f64, f64) {
    let mut aligner = Aligner::<DnaAlphabet>::default();
    aligner.set_costs_parse(costs).unwrap();
    aligner.set_min_length_strategy(min_length_strategy);
    let result = aligner
        .align(
            "reference",
            reference,
            "query",
            query,
            None,
            &[],
            None,
            None,
            false,
            TSMUncertaintyRangeExtensionMode::None,
        )
        .unwrap();
    println!("{result}");

    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = &result
    else {
        panic!("No alignment found");
    };
    let config = TemplateSwitchConfig::<DnaAlphabet, U64Cost>::from_str(costs).unwrap();
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(reference).unwrap();
    let query = VectorGenome::<DnaAlphabet>::from_slice_u8(query).unwrap();
    let recomputed_cost = alignment.compute_cost(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        0,
        0,
        &config,
        false,
    );

    (
        statistics.cost.raw(),
        statistics.template_switch_amount.raw(),
        recomputed_cost.as_primitive() as f64,
    )
}

/// The query contains a tandem duplication of a 10 character segment `Y` of the reference.
const TANDEM_DUPLICATION_REFERENCE: &[u8] = b"GCTAAAGACAATTACATAAC GTCAGCACGA AACTTGCATGGTCA";
const TANDEM_DUPLICATION_QUERY: &[u8] =
    b"GCTAAAGACAATTACATAAC GTCAGCACGA GTCAGCACGA AACTTGCATGGTCA";

/// The query contains 12 additional `A`s inside a run of 12 `A`s of the reference.
const HOMOPOLYMER_INSERTION_REFERENCE: &[u8] = b"GCTAAAGACAATTACATAAC AAAAAAAAAAAA AACTTGCATGGTCA";
const HOMOPOLYMER_INSERTION_QUERY: &[u8] =
    b"GCTAAAGACAATTACATAAC AAAAAAAAAAAA AAAAAAAAAAAA AACTTGCATGGTCA";

fn without_spaces(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().copied().filter(|c| *c != b' ').collect()
}

/// Issue #166: the initial entrance of a forward template switch has offset -1 or 1,
/// but its cost used to be computed for offset 0.
///
/// The query contains a tandem duplication of 10 characters, and only forward QQ template switches are allowed,
/// with base cost 2 and offset cost 7 for offsets in `-20..0` (zero only for offset zero, which a forward template switch never has).
/// Any alignment without template switch needs at least 10 insertions, costing at least `3 + 9 = 12`.
/// Any template switch costs at least `2 + 7 = 9`, and the QQ template switch with offset -10 that copies
/// the duplicated segment costs exactly 9, since its inner matches exactly and all other costs are zero.
/// Hence the optimal cost is 9.
/// With the bug, the offset cost was charged as `cost(0) + cost(-10) - cost(-1) = 0`, giving cost 2.
#[test]
fn issue_166_forward_qq_template_switch_pays_offset_costs() {
    let costs = forward_ts_offset_test_costs(
        "2",
        "inf",
        INFINITE_OFFSET_COSTS,
        FINITE_FORWARD_OFFSET_COSTS,
    );
    let result = align_forward_ts_offset_test(
        &costs,
        &without_spaces(TANDEM_DUPLICATION_REFERENCE),
        &without_spaces(TANDEM_DUPLICATION_QUERY),
    );
    assert_eq!(result, (9.0, 1.0, 9.0));
}

/// Same as [`issue_166_forward_qq_template_switch_pays_offset_costs`], but only forward QR template switches are allowed,
/// such that the RQ/QR offset costs are used.
/// The QR template switch with offset -10 copies the segment from the reference and costs exactly 9.
#[test]
fn issue_166_forward_qr_template_switch_pays_offset_costs() {
    let costs = forward_ts_offset_test_costs(
        "inf",
        "2",
        FINITE_FORWARD_OFFSET_COSTS,
        INFINITE_OFFSET_COSTS,
    );
    let result = align_forward_ts_offset_test(
        &costs,
        &without_spaces(TANDEM_DUPLICATION_REFERENCE),
        &without_spaces(TANDEM_DUPLICATION_QUERY),
    );
    assert_eq!(result, (9.0, 1.0, 9.0));
}

/// Issue #166: if the offset costs are finite only for offset zero, then forward template switches are impossible.
///
/// Only forward QQ template switches with base cost 1 are enabled, but their offset costs are infinite for all offsets except zero,
/// which a forward template switch never has. Hence the optimal alignment has no template switch,
/// and since the query contains 12 additional `A`s, its optimal cost is a single gap of length 12, costing `3 + 11 = 14`.
/// With the bug, the initial entrances at offsets -1 and 1 were charged the finite cost of offset zero,
/// such that a QQ template switch with offset -1 copying the `A`-run for cost 1 was found.
#[test]
fn issue_166_forward_qq_template_switch_with_infinite_offset_costs() {
    let costs = forward_ts_offset_test_costs(
        "1",
        "inf",
        INFINITE_OFFSET_COSTS,
        ZERO_ONLY_FORWARD_OFFSET_COSTS,
    );
    let result = align_forward_ts_offset_test(
        &costs,
        &without_spaces(HOMOPOLYMER_INSERTION_REFERENCE),
        &without_spaces(HOMOPOLYMER_INSERTION_QUERY),
    );
    assert_eq!(result, (14.0, 0.0, 14.0));
}

/// Same as [`issue_166_forward_qq_template_switch_with_infinite_offset_costs`], but only forward QR template switches are enabled,
/// such that the RQ/QR offset costs are used.
#[test]
fn issue_166_forward_qr_template_switch_with_infinite_offset_costs() {
    let costs = forward_ts_offset_test_costs(
        "inf",
        "1",
        ZERO_ONLY_FORWARD_OFFSET_COSTS,
        INFINITE_OFFSET_COSTS,
    );
    let result = align_forward_ts_offset_test(
        &costs,
        &without_spaces(HOMOPOLYMER_INSERTION_REFERENCE),
        &without_spaces(HOMOPOLYMER_INSERTION_QUERY),
    );
    assert_eq!(result, (14.0, 0.0, 14.0));
}

/// The min length strategies that keep the classic aligner exact.
/// [`MinLengthStrategySelector::PreprocessFilter`] is a heuristic that discards template switches whose inner
/// does not start with an exact match of minimum length, so it is excluded.
fn exact_min_length_strategies() -> [MinLengthStrategySelector; 4] {
    [
        MinLengthStrategySelector::None,
        MinLengthStrategySelector::Lookahead,
        MinLengthStrategySelector::PreprocessPrice,
        MinLengthStrategySelector::PreprocessLookahead,
    ]
}

/// Compute the cost of the given alignment independently of the A* search with
/// [`Alignment::compute_cost_and_lengths`](crate::alignment_result::alignment::Alignment::compute_cost_and_lengths),
/// including the geometry costs, and return it together with the number of consumed reference and query characters.
fn recompute_cost_and_lengths(
    costs: &str,
    reference: &[u8],
    query: &[u8],
    alignment: &[(usize, AlignmentType)],
) -> (f64, usize, usize) {
    let config = TemplateSwitchConfig::<DnaAlphabet, U64Cost>::from_str(costs).unwrap();
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(reference).unwrap();
    let query = VectorGenome::<DnaAlphabet>::from_slice_u8(query).unwrap();
    let alignment = Alignment::from(alignment.to_vec());
    let (cost, reference_length, query_length) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        0,
        0,
        &config,
        false,
    );
    (cost.as_primitive() as f64, reference_length, query_length)
}

/// Check that `witness` is a valid alignment of `reference` and `query` of cost `expected_cost`,
/// and that the classic aligner finds an alignment of cost exactly `expected_cost` with each exact min length strategy,
/// i.e. that the witness is optimal.
///
/// The witness is checked first, so a reported cost above `expected_cost` means that the classic aligner misses
/// an alignment that its own cost model accepts.
/// A reported cost below `expected_cost` (which is recomputed independently of the A* search)
/// means that the witness is not optimal.
fn assert_classic_aligner_finds_witness_cost(
    costs: &str,
    reference: &[u8],
    query: &[u8],
    witness: &[(usize, AlignmentType)],
    expected_cost: f64,
) {
    assert_eq!(
        recompute_cost_and_lengths(costs, reference, query, witness),
        (expected_cost, reference.len(), query.len()),
        "the witness alignment is invalid or has an unexpected cost"
    );

    let results: Vec<_> = exact_min_length_strategies()
        .into_iter()
        .map(|min_length_strategy| {
            let name = format!("{min_length_strategy:?}");
            (
                name,
                align_and_recompute_cost(costs, min_length_strategy, reference, query),
            )
        })
        .collect();
    for (name, (cost, _, recomputed_cost)) in &results {
        assert_eq!(
            cost, recomputed_cost,
            "{name}: reported cost differs from recomputed cost"
        );
    }
    assert!(
        results
            .iter()
            .all(|(_, (cost, _, _))| *cost == expected_cost),
        "expected cost {expected_cost} with all exact min length strategies, but got {results:?}"
    );
}

/// A reverse template switch with descendant `descendant`, ancestor `ancestor` and the given first offset,
/// whose inner consists of 9 secondary matches followed by one secondary insertion,
/// and whose anti-descendant gap is 10.
fn reverse_ts_with_trailing_insertion(
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
    first_offset: isize,
) -> [(usize, AlignmentType); 4] {
    [
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                descendant,
                ancestor,
                direction: TemplateSwitchDirection::Reverse,
                uncertainty_range: None,
                first_offset,
            },
        ),
        (9, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondaryInsertion),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 10,
            },
        ),
    ]
}

/// The tschainalign end-to-end test case `inversion_9_start`.
/// The first 9 characters of the query are the reverse complement of the first 9 characters of the reference,
/// followed by an extra `C`, and the last 30 characters are equal.
const INVERSION_9_START_REFERENCE: &[u8] = b"AATCGGAAC C TGGTACTAGGCATAAGACGATTCGTGCGAA";
const INVERSION_9_START_QUERY: &[u8] = b"GTTCCGATT C TGGTACTAGGCATAAGACGATTCGTGCGAA";

/// The alignment `[TSQRR:[-]:[-]:9:9=1I:10]30=` reported by tschainalign for `inversion_9_start`.
fn inversion_9_start_witness(
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
) -> Vec<(usize, AlignmentType)> {
    let mut witness = reverse_ts_with_trailing_insertion(descendant, ancestor, 9).to_vec();
    witness.push((30, AlignmentType::PrimaryMatch));
    witness
}

/// The alignment reported by tschainalign for `inversion_9_start` is valid under the cost model of the classic aligner
/// and costs 5.
///
/// It is a reverse QR template switch entering at `(0, 0)` with first offset 9 (base cost 2, offset cost 0).
/// Its ancestor runs backwards through `reference[0..9]`, whose reverse complement matches `query[0..9]`,
/// and then the ancestor is exhausted, so the 10th descendant character `query[9]` is a secondary insertion (cost 3).
/// The descendant length is 10, which is the minimum template switch length (length cost 0),
/// and the anti-descendant gap is 10 (length difference 0, reverse anti-descendant gap cost 0).
/// The remaining 30 characters match.
#[test]
fn ts_at_sequence_start_with_ancestor_shorter_than_min_length_is_valid() {
    let reference = without_spaces(INVERSION_9_START_REFERENCE);
    let query = without_spaces(INVERSION_9_START_QUERY);
    let witness = inversion_9_start_witness(
        TemplateSwitchDescendant::Query,
        TemplateSwitchAncestor::Reference,
    );
    assert_eq!(
        recompute_cost_and_lengths(TS_END_TEST_COSTS, &reference, &query, &witness),
        (5.0, 40, 40)
    );
}

/// The classic aligner must find an alignment of cost 5, which is optimal,
/// like the one of [`ts_at_sequence_start_with_ancestor_shorter_than_min_length_is_valid`].
///
/// It used to report cost 8 instead (`[TSRQR:[0,0]:[0,0]:10:1D9=1I:10]30=`), with every min length strategy:
/// the template switch entrance only generated the secondary root of a reverse template switch if the ancestor index
/// was at least the minimum template switch length (10), but here the ancestor index is `0 + 9 = 9`,
/// because the inner consumes only 9 ancestor characters (the 10th descendant character is a secondary insertion).
#[test]
fn ts_at_sequence_start_with_ancestor_shorter_than_min_length() {
    let reference = without_spaces(INVERSION_9_START_REFERENCE);
    let query = without_spaces(INVERSION_9_START_QUERY);
    assert_classic_aligner_finds_witness_cost(
        TS_END_TEST_COSTS,
        &reference,
        &query,
        &inversion_9_start_witness(
            TemplateSwitchDescendant::Query,
            TemplateSwitchAncestor::Reference,
        ),
        5.0,
    );
}

/// Same as [`ts_at_sequence_start_with_ancestor_shorter_than_min_length`], but with reference and query swapped,
/// such that the witness is a reverse RQ template switch.
#[test]
fn ts_at_sequence_start_with_ancestor_shorter_than_min_length_swapped() {
    let reference = without_spaces(INVERSION_9_START_QUERY);
    let query = without_spaces(INVERSION_9_START_REFERENCE);
    assert_classic_aligner_finds_witness_cost(
        TS_END_TEST_COSTS,
        &reference,
        &query,
        &inversion_9_start_witness(
            TemplateSwitchDescendant::Reference,
            TemplateSwitchAncestor::Query,
        ),
        5.0,
    );
}

/// Control for [`ts_at_sequence_start_with_ancestor_shorter_than_min_length`]:
/// with a common prefix of 5 characters, the witness template switch enters at `(5, 5)` with first offset 9,
/// so its ancestor index is `5 + 9 = 14`, which is at least the minimum template switch length.
/// Then the classic aligner found the optimal cost-5 alignment even before the ancestor index of template switch
/// entrances stopped being bounded by the minimum template switch length.
#[test]
fn ts_with_ancestor_shorter_than_min_length_after_prefix() {
    let prefix = b"CTGAC";
    let reference = [
        prefix.as_slice(),
        &without_spaces(INVERSION_9_START_REFERENCE),
    ]
    .concat();
    let query = [prefix.as_slice(), &without_spaces(INVERSION_9_START_QUERY)].concat();
    let mut witness = vec![(5, AlignmentType::PrimaryMatch)];
    witness.extend(inversion_9_start_witness(
        TemplateSwitchDescendant::Query,
        TemplateSwitchAncestor::Reference,
    ));
    assert_classic_aligner_finds_witness_cost(TS_END_TEST_COSTS, &reference, &query, &witness, 5.0);
}

/// Like [`ts_at_sequence_start_with_ancestor_shorter_than_min_length`], but the template switch does not start
/// at the start of the alignment: it enters at `(20, 20)` with first offset -11, such that the ancestor runs backwards
/// through `reference[0..9]` (ancestor index `20 - 11 = 9`).
///
/// The query copies the reverse complement of the first 9 characters of the reference plus an extra `C`
/// in place of the 10 characters `Z` of the reference.
/// The witness costs 5 (base cost 2, secondary insertion 3, all geometry costs 0), which is optimal.
/// The classic aligner used to report cost 7, since it neither walked the first offset below ancestor index 10
/// nor generated a secondary root with ancestor index below 10.
/// Hence the bug was not about template switches at the start of the alignment,
/// but about the ancestor index of reverse template switches being smaller than the minimum template switch length.
#[test]
fn ts_with_ancestor_shorter_than_min_length_at_ancestor_start() {
    //               [ A9  ]                    [  Z   ]
    let reference = b"AATCGGAAC CTAGTCAGGTA ACAGTAGCAG GCATAAGACGATTCGTGCGA";
    //               [ A9  ]                    [rc(A9)]C
    let query = b"AATCGGAAC CTAGTCAGGTA GTTCCGATT C GCATAAGACGATTCGTGCGA";
    let mut witness = vec![(20, AlignmentType::PrimaryMatch)];
    witness.extend(reverse_ts_with_trailing_insertion(
        TemplateSwitchDescendant::Query,
        TemplateSwitchAncestor::Reference,
        -11,
    ));
    witness.push((20, AlignmentType::PrimaryMatch));
    assert_classic_aligner_finds_witness_cost(
        TS_END_TEST_COSTS,
        &without_spaces(reference),
        &without_spaces(query),
        &witness,
        5.0,
    );
}

/// The forward analogue of [`ts_with_ancestor_shorter_than_min_length_at_ancestor_start`]:
/// a forward QR template switch whose ancestor ends at the end of the reference and consumes only 9 characters.
///
/// Only forward QR template switches with base cost 2 are allowed, and their offset costs are 1 for offsets in
/// `1..21` (see [`forward_ts_offset_test_costs`]).
/// The query contains an insertion of 10 characters `X`, which are the last 9 characters of the reference plus an extra `C`.
/// The witness is a forward QR template switch entering at `(20, 20)` with first offset 2
/// (ancestor `reference[22..31]`), whose inner has 9 secondary matches followed by a secondary insertion,
/// and whose anti-descendant gap is 0. It costs `2 + 1 + 3 = 6`, which is optimal.
/// The classic aligner used to report cost 7, since it required `ancestor_index + 10 <= reference.len()` for forward
/// template switches, but here `22 + 10 > 31`.
#[test]
fn forward_ts_with_ancestor_shorter_than_min_length_at_ancestor_end() {
    let costs = forward_ts_offset_test_costs(
        "inf",
        "2",
        " -inf -20 0 1 21\n  inf   1 0 1 inf",
        INFINITE_OFFSET_COSTS,
    );
    let reference = b"GCTAAAGACAATTACATAAC GT CAGCACGAT";
    //                                         [    X    ]
    let query = b"GCTAAAGACAATTACATAAC CAGCACGAT C GT CAGCACGAT";
    let witness = [
        (20, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                descendant: TemplateSwitchDescendant::Query,
                ancestor: TemplateSwitchAncestor::Reference,
                direction: TemplateSwitchDirection::Forward,
                uncertainty_range: None,
                first_offset: 2,
            },
        ),
        (9, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondaryInsertion),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 0,
            },
        ),
        (11, AlignmentType::PrimaryMatch),
    ];
    assert_classic_aligner_finds_witness_cost(
        &costs,
        &without_spaces(reference),
        &without_spaces(query),
        &witness,
        6.0,
    );
}

/// Align with every min length strategy, including the heuristic [`MinLengthStrategySelector::PreprocessFilter`],
/// check that each reported cost equals the independently recomputed cost, and return the costs by strategy name.
fn align_with_all_min_length_strategies(
    costs: &str,
    reference: &[u8],
    query: &[u8],
) -> Vec<(String, f64)> {
    exact_min_length_strategies()
        .into_iter()
        .chain([MinLengthStrategySelector::PreprocessFilter])
        .map(|min_length_strategy| {
            let name = format!("{min_length_strategy:?}");
            let (cost, _, recomputed_cost) =
                align_and_recompute_cost(costs, min_length_strategy, reference, query);
            assert_eq!(
                cost, recomputed_cost,
                "{name}: reported cost differs from recomputed cost"
            );
            (name, cost)
        })
        .collect()
}

/// Both sequences are shorter than the minimum template switch length (10), so no template switch is possible.
/// The preprocessed min length strategies used to panic with an integer underflow
/// (`len - template_switch_min_length`) on such inputs.
/// The optimal alignment has a single substitution (cost 2).
#[test]
fn min_length_strategies_with_both_sequences_shorter_than_min_length() {
    let results = align_with_all_min_length_strategies(TS_END_TEST_COSTS, b"CCCTAGTT", b"CCCAAGTT");
    assert!(
        results.iter().all(|(_, cost)| *cost == 2.0),
        "expected cost 2 with all min length strategies, but got {results:?}"
    );
}

/// `A + rc(A)` for a 10 character `A`, so it is its own reverse complement.
const PALINDROME_20: &[u8] = b"GATTACAGGC GCCTGTAATC";
/// A sequence shorter than the minimum template switch length (10).
const SHORT_5: &[u8] = b"TTGCA";

/// The reference is [`PALINDROME_20`] followed by the 5 characters of the query,
/// so the query is shorter than the minimum template switch length.
///
/// The optimal alignment is a reverse RR template switch at `(0, 0)` with first offset 20, whose inner matches all
/// 20 characters of the palindrome exactly (base cost 2, all geometry costs 0), followed by 5 primary matches.
/// Without a template switch, the 20 extra reference characters would cost at least `3 + 19 = 22`.
///
/// The preprocessed min length strategies used to panic on such inputs.
/// Since the long reference contains an exact match of minimum length with itself, the short query must not make them
/// treat the template switch as a non-match. Otherwise `PreprocessPrice` would use an inadmissible lower bound and
/// `PreprocessFilter` would drop the optimal template switch, even though its inner is an exact match.
#[test]
fn min_length_strategies_with_query_shorter_than_min_length() {
    let reference = [without_spaces(PALINDROME_20).as_slice(), SHORT_5].concat();
    let query = SHORT_5;
    let witness = [
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Reference,
                direction: TemplateSwitchDirection::Reverse,
                uncertainty_range: None,
                first_offset: 20,
            },
        ),
        (20, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 0,
            },
        ),
        (5, AlignmentType::PrimaryMatch),
    ];
    assert_classic_aligner_finds_witness_cost(TS_END_TEST_COSTS, &reference, query, &witness, 2.0);

    let results = align_with_all_min_length_strategies(TS_END_TEST_COSTS, &reference, query);
    assert!(
        results.iter().all(|(_, cost)| *cost == 2.0),
        "expected cost 2 with all min length strategies, but got {results:?}"
    );
}

/// Same as [`min_length_strategies_with_query_shorter_than_min_length`], but with reference and query swapped,
/// such that the optimal alignment uses a reverse QQ template switch.
#[test]
fn min_length_strategies_with_reference_shorter_than_min_length() {
    let reference = SHORT_5;
    let query = [without_spaces(PALINDROME_20).as_slice(), SHORT_5].concat();
    let witness = [
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                descendant: TemplateSwitchDescendant::Query,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
                uncertainty_range: None,
                first_offset: 20,
            },
        ),
        (20, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 0,
            },
        ),
        (5, AlignmentType::PrimaryMatch),
    ];
    assert_classic_aligner_finds_witness_cost(TS_END_TEST_COSTS, reference, &query, &witness, 2.0);

    let results = align_with_all_min_length_strategies(TS_END_TEST_COSTS, reference, &query);
    assert!(
        results.iter().all(|(_, cost)| *cost == 2.0),
        "expected cost 2 with all min length strategies, but got {results:?}"
    );
}

/// If all `Length` costs are infinite, then template switches are disabled and the minimum template switch length
/// is `usize::MAX`. The preprocessed min length strategies must not underflow on this either.
/// The optimal alignment is the same as without template switches: 5 deleted characters (`3 + 4 = 7`).
#[test]
fn min_length_strategies_with_template_switches_disabled() {
    let costs = TS_END_TEST_COSTS.replace("Length\n   0 10 200\n inf  0 inf", "Length\n   0\n inf");
    assert_ne!(costs, TS_END_TEST_COSTS);
    let reference = [without_spaces(PALINDROME_20).as_slice(), SHORT_5].concat();
    let query = without_spaces(PALINDROME_20);
    let results = align_with_all_min_length_strategies(&costs, &reference, &query);
    assert!(
        results.iter().all(|(_, cost)| *cost == 7.0),
        "expected cost 7 with all min length strategies, but got {results:?}"
    );
}

/// Issue #250: the cost limit is inclusive, so an alignment whose cost equals the limit is found.
///
/// The default aligner searches label-correcting (it maximises the total template switch length
/// as a secondary criterion), which used to report the cost limit as exceeded when the limit
/// pruned a node after the target had already been found.
#[test]
fn cost_limit_is_inclusive() {
    let aligner = Aligner::<DnaAlphabetOrN>::new();
    let align = |reference: &[u8], query: &[u8], cost_limit| {
        aligner
            .align(
                "reference",
                reference,
                "query",
                query,
                None,
                &[],
                cost_limit,
                None,
                true,
                TSMUncertaintyRangeExtensionMode::EqualCostIgnoreGeometry,
            )
            .unwrap()
    };

    for (reference, query) in [
        (&b"ACGTACGT"[..], &b"ACGACGT"[..]),
        (b"ACGTACGTAC", b"ACGTTCGTAC"),
        (b"AAAAAAAAAACCCCCCCCCCCC", b"AAAAAAAAAATTTTTTTTTTTC"),
    ] {
        let unlimited = align(reference, query, None);
        let AlignmentResult::WithTarget { .. } = &unlimited else {
            panic!("no alignment without cost limit: {unlimited:#?}");
        };
        let cost = unlimited.statistics().cost.raw();
        assert!(cost > 0.0 && cost.fract() == 0.0, "unexpected cost {cost}");
        let cost = cost as u64;

        for cost_limit in [cost, cost + 1] {
            let limited = align(reference, query, Some(cost_limit));
            let AlignmentResult::WithTarget { .. } = &limited else {
                panic!(
                    "cost limit {cost_limit} rejected an alignment of cost {cost}: {limited:#?}"
                );
            };
            assert_eq!(limited.statistics().cost.raw(), cost as f64);
        }

        let limited = align(reference, query, Some(cost - 1));
        assert!(
            matches!(limited, AlignmentResult::WithoutTarget { .. }),
            "cost limit {} accepted an alignment of cost {cost}: {limited:#?}",
            cost - 1,
        );
    }
}
