use num_traits::real::Real;

use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;

use crate::a_star_aligner::alignment_geometry::AlignmentCoordinates;

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
        );
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
    let result = aligner.align(
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
    );
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
