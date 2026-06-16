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
        true);
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
    let _tsa_res = tsa.align("ref", r, "qry", q, None, &[], None, None, false);
}
