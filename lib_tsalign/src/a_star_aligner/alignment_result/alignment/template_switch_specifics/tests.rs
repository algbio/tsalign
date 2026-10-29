use std::sync::LazyLock;

use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::{
        alphabet::Alphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
};
use generic_a_star::cost::U64Cost;

use crate::{
    a_star_aligner::{
        alignment_result::alignment::Alignment,
        template_switch_distance::{
            AlignmentType, TemplateSwitchAncestor, TemplateSwitchDescendant,
            TemplateSwitchDirection,
        },
    },
    config::{BaseCost, TemplateSwitchConfig},
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
};

static START_REFERENCE: &[u8] = b"AGAGAGCTCTAA";
static START_QUERY: &[u8] = b"AGAGAGCTTTAA";
static START_COSTS: LazyLock<Vec<U64Cost>> = LazyLock::new(|| {
    [
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&-6)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&2)
            + CONFIG.length_costs.evaluate(&2)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&-4)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&3)
            + CONFIG.length_costs.evaluate(&3)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&-2)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&4)
            + CONFIG.length_costs.evaluate(&4)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&0)
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'G').unwrap(),
                DnaAlphabet::ascii_to_character(b'T').unwrap(),
            )
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&5)
            + CONFIG.length_costs.evaluate(&5)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&2)
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'G').unwrap(),
                DnaAlphabet::ascii_to_character(b'T').unwrap(),
            )
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'A').unwrap(),
                DnaAlphabet::ascii_to_character(b'C').unwrap(),
            )
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&6)
            + CONFIG.length_costs.evaluate(&6)
            + CONFIG.length_difference_costs.evaluate(&0),
    ]
    .to_vec()
});
static START_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
    &[
        (6, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: -6,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (2, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 2,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
    &[
        (5, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: -4,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (3, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 3,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
    &[
        (4, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: -2,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (4, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 4,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
    &[
        (3, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 0,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (1, AlignmentType::SecondarySubstitution),
        (4, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 5,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
    &[
        (2, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 2,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (2, AlignmentType::SecondarySubstitution),
        (4, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 6,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
];

static END_REFERENCE: &[u8] = b"AACTCTAGAGAG";
static END_QUERY: &[u8] = b"AATTCTAGAGAG";
static END_COSTS: LazyLock<Vec<U64Cost>> = LazyLock::new(|| {
    [
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&10)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&2)
            + CONFIG.length_costs.evaluate(&2)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&10)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&3)
            + CONFIG.length_costs.evaluate(&3)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&10)
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&4)
            + CONFIG.length_costs.evaluate(&4)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&10)
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'A').unwrap(),
                DnaAlphabet::ascii_to_character(b'C').unwrap(),
            )
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&5)
            + CONFIG.length_costs.evaluate(&5)
            + CONFIG.length_difference_costs.evaluate(&0),
        CONFIG.base_cost.rqr
            + CONFIG.rq_qr_offset_costs.evaluate(&10)
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'A').unwrap(),
                DnaAlphabet::ascii_to_character(b'C').unwrap(),
            )
            + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                DnaAlphabet::ascii_to_character(b'G').unwrap(),
                DnaAlphabet::ascii_to_character(b'T').unwrap(),
            )
            + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&6)
            + CONFIG.length_costs.evaluate(&6)
            + CONFIG.length_difference_costs.evaluate(&0),
    ]
    .to_vec()
});
static END_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
    &[
        (1, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 10,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (2, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 2,
            },
        ),
        (6, AlignmentType::PrimaryMatch),
    ],
    &[
        (1, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 10,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (3, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 3,
            },
        ),
        (5, AlignmentType::PrimaryMatch),
    ],
    &[
        (1, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 10,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (4, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 4,
            },
        ),
        (4, AlignmentType::PrimaryMatch),
    ],
    &[
        (1, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 10,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (4, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondarySubstitution),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 5,
            },
        ),
        (3, AlignmentType::PrimaryMatch),
    ],
    &[
        (1, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: 10,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Reference,
                ancestor: TemplateSwitchAncestor::Query,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (4, AlignmentType::SecondaryMatch),
        (2, AlignmentType::SecondarySubstitution),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 6,
            },
        ),
        (2, AlignmentType::PrimaryMatch),
    ],
];

static CONFIG: LazyLock<TemplateSwitchConfig<DnaAlphabet, U64Cost>> =
    LazyLock::new(|| TemplateSwitchConfig {
        left_flank_length: 0,
        right_flank_length: 0,
        template_switch_min_length: 3,
        base_cost: BaseCost {
            rrf: 10u64.into(),
            rqf: 100u64.into(),
            qrf: 1000u64.into(),
            qqf: 10000u64.into(),
            rrr: 100000u64.into(),
            rqr: 1000000u64.into(),
            qrr: 10000000u64.into(),
            qqr: 100000000u64.into(),
        },
        primary_edit_costs: GapAffineAlignmentCostTable::new(
            "",
            [0u64, 2, 2, 2, 2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 0]
                .map(Into::into)
                .to_vec(),
            [3u64, 3, 3, 3].map(Into::into).to_vec(),
            [1u64, 1, 1, 1].map(Into::into).to_vec(),
        ),
        secondary_forward_edit_costs: GapAffineAlignmentCostTable::new(
            "",
            [0u64, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0]
                .map(Into::into)
                .to_vec(),
            [3u64, 6, 6, 6].map(Into::into).to_vec(),
            [1u64, 1, 1, 1].map(Into::into).to_vec(),
        ),
        secondary_reverse_edit_costs: GapAffineAlignmentCostTable::new(
            "",
            [0u64, 5, 5, 5, 5, 0, 5, 5, 5, 5, 0, 5, 5, 5, 5, 0]
                .map(Into::into)
                .to_vec(),
            [3u64, 7, 7, 7].map(Into::into).to_vec(),
            [1u64, 1, 1, 1].map(Into::into).to_vec(),
        ),
        left_flank_edit_costs: GapAffineAlignmentCostTable::new_zero(),
        right_flank_edit_costs: GapAffineAlignmentCostTable::new_zero(),
        rq_qr_offset_costs: CostFunction::try_from(
            (-20..=20)
                .map(|i| (i, U64Cost::from(17 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
        rr_qq_offset_costs: CostFunction::try_from(
            (-20..=20)
                .map(|i| (i, U64Cost::from(17 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
        length_costs: CostFunction::try_from(
            (0..=20)
                .map(|i| (i, U64Cost::from(19 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
        length_difference_costs: CostFunction::try_from(
            (-20..=20)
                .map(|i| (i, U64Cost::from(23 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
        forward_anti_descendant_gap_costs: CostFunction::try_from(
            (-20..=20)
                .map(|i| (i, U64Cost::from(29 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
        reverse_anti_descendant_gap_costs: CostFunction::try_from(
            (-20..=20)
                .map(|i| (i, U64Cost::from(31 * u64::try_from(i + 21).unwrap())))
                .collect::<Vec<_>>(),
        )
        .unwrap(),
    });

#[test]
fn move_template_switch_start_backwards() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
    let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
    let mut alignment = Alignment::from(START_ALIGNMENTS[0].to_vec());
    let (cost, len_a, len_b) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        2,
        2,
        &CONFIG,
        false,
    );
    assert_eq!(cost, START_COSTS[0]);

    for (expected_alignment, expected_cost) in START_ALIGNMENTS[1..].iter().zip(&START_COSTS[1..]) {
        assert!(alignment.move_template_switch_start_backwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            2,
            2,
            &mut 1,
        ));
        assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        assert_eq!(
            alignment.compute_cost_and_lengths(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &CONFIG,
                false,
            ),
            (*expected_cost, len_a, len_b),
        );
    }
}

#[test]
fn move_template_switch_start_forwards() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
    let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
    let mut alignment = Alignment::from(START_ALIGNMENTS.last().unwrap().to_vec());
    let (cost, len_a, len_b) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        2,
        2,
        &CONFIG,
        false,
    );
    assert_eq!(cost, *START_COSTS.last().unwrap());

    for (expected_alignment, expected_cost) in START_ALIGNMENTS
        .iter()
        .zip(START_COSTS.iter())
        .rev()
        .skip(1)
    {
        assert!(alignment.move_template_switch_start_forwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            2,
            2,
            &mut 1,
        ));
        assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        assert_eq!(
            alignment.compute_cost_and_lengths(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &CONFIG,
                false,
            ),
            (*expected_cost, len_a, len_b),
        );
    }
}

#[test]
fn move_template_switch_end_backwards() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
    let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
    let mut alignment = Alignment::from(END_ALIGNMENTS.last().unwrap().to_vec());
    let (cost, len_a, len_b) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        1,
        1,
        &CONFIG,
        false,
    );
    assert_eq!(cost, *END_COSTS.last().unwrap());

    for (expected_alignment, expected_cost) in
        END_ALIGNMENTS.iter().zip(END_COSTS.iter()).rev().skip(1)
    {
        assert!(alignment.move_template_switch_end_backwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            1,
            1,
            1,
        ));
        assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        assert_eq!(
            alignment.compute_cost_and_lengths(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                &CONFIG,
                false,
            ),
            (*expected_cost, len_a, len_b),
        );
    }
}

#[test]
fn move_template_switch_end_forwards() {
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
    let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
    let mut alignment = Alignment::from(END_ALIGNMENTS[0].to_vec());
    let (cost, len_a, len_b) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        1,
        1,
        &CONFIG,
        false,
    );
    assert_eq!(cost, *END_COSTS.first().unwrap());

    for (expected_alignment, expected_cost) in END_ALIGNMENTS[1..].iter().zip(&END_COSTS[1..]) {
        assert!(alignment.move_template_switch_end_forwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            1,
            1,
            1,
        ));
        assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        assert_eq!(
            alignment.compute_cost_and_lengths(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                &CONFIG,
                false,
            ),
            (*expected_cost, len_a, len_b),
        );
    }
}

#[test]
fn ari_41() {
    let reference = b"ACTCAAGTATTTAGAATTACTTCTGAGCACCTAAGACATATCAGTGAAACATCAATTTTGGCTTCCCCATGCTTTCTCTGGGTCTCCAAGTTCTCTGTATAAAGCTGAAGCTAGGAATGGAATGATTAGTGACTCTAAGGGAAAAATTTCTGCAGAGGTGGGCTATCCATCTGCCATCATAGGTGCTATCCTTGCAGTTGGCTCATGAGGTGTTTCAGATACTCTGCAGTGGTGGCTGTGCATTCTTAGCATGCCAGAGACCTGGAATACAGCCTCCTGCTTTCCCTCACTCAAGCCATGCCCTGAGAACAAGACTGGCCAGCCAGATTGTGAGCAGGTCTCTGTttcttacctattggtctactgaccccaaagtcacattttcgcctcctgtttctcaaaccttcctgcactgaaagataataaattgtgagcagagaagtcaaaattcctttggactctggttttCCAAAGATACCACCTAAAAGTCCAGTTTTCCCAGCAATATGCTGGGCCATGTGTACATTTCAAACTTCGGCTGAAAACACATTTCACAACAGCTATCATAAGGGAACAGAGAAAACACAACCATTGTTTAAAAGCCAGAAGTGGGTGGTGGTCCTAGTCAAGAATTGTTGGAAGACAAAACAATGTCGTCTTGTTGGTCTATGAACCACCATTGGGAAGCAGTGTTCTCCACcagtggttcccaaaatgggctgcacacactcaggggatgcttaagatcacccactgtaacatggatagaaagcattaagacctttcttcctatacgttcttaaaaaaaaaaaggaacaagaatgaagctttacttctttacttggtccacatactagatggccctgcctcaagtttagctgatgaaatatcaaggtgccaaaacagaaaggagttggcagtaatgccttacccttcattcacttgtgacacatcatgttgtgttgcaatttctatttgtctggttaatggatttatcacttatc".to_ascii_uppercase();
    let query = b"ACTCAAGTATTTAGAATTACTTCTGAGCACCTAAGACATATCAGTGAAACATCAATTTTGGCTTCCCCATGCTTTCTCTGGGTCTCCAAGTTCTCTGTATAAAGCTGAAGCTAGGAATGGAATGATTAGTGACTCTAAGGGAAAAATTTCTGCAGAGGTGGGCTATCCATCTGCCATCATAGGTGCTATCCTTGCAGTTGGCTCATGAGGTGTTTCAGATACTCTGCAGTGGTGGCTGTGCATTCTTAGCATGCCAGAGACCTGGAATACAGCCTCCTGCTTTCCCTCACTCAAGCCATGCCCTGAGAACAAGACTGGCCAGCCAGATTGTGAGCAGGTCTCTGTttcttacctattggtctactgaccccaaagtcacattttcgcctcctgtttctcaaaccttcctgcactgaaagataataaattgtgagcagagaagtcaaaattcctttggactctggttttCCAAAGATACCACCTAAAAGTCCAGTTTTCCCTGGACTTTTCAATATGCTGGGCCATGTGTACATTTCAAACTTCGGCTGAAAACACATTTCACAACAGCTATCATAAGGGAACAGAGAAAACACAACCATTGTTTAAAAGCCAGAAGTGGGTGGTGGTCCTAGTCAAGAATTGTTGGAAGACAAAACAATGTCGTCTTGTTGGTCTATGAACCACCATTGGGAAGCAGTGTTCTCCACcagtggttcccaaaatgggctgcacacactcaggggatgcttaagatcacccactgtaacatggatagaaagcattaagacctttcttcctatacgttcttaaaaaaaaaaaggaacaagaatgaagctttacttctttacttggtccacatactagatggccctgcctcaagtttagctgatgaaatatcaaggtgccaaaacagaaaggagttggcagtaatgccttacccttcattcacttgtgacacatcatgttgtgttgcaatttctatttgtctggttaatggatttatcacttatc".to_ascii_uppercase();
    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(&reference).unwrap();
    let query = VectorGenome::from_slice_u8(&query).unwrap();
    // 499=[TSQRR:[-]:[-]:-6:10=:3]502=
    let alignment = Alignment::from(vec![
        (499, AlignmentType::PrimaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchEntrance {
                first_offset: -6,
                uncertainty_range: None,
                descendant: TemplateSwitchDescendant::Query,
                ancestor: TemplateSwitchAncestor::Reference,
                direction: TemplateSwitchDirection::Reverse,
            },
        ),
        (10, AlignmentType::SecondaryMatch),
        (
            1,
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap: 3,
            },
        ),
        (502, AlignmentType::PrimaryMatch),
    ]);
    println!("Alignment: {alignment}");

    let (_, len_a, len_b) = alignment.compute_cost_and_lengths(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        0,
        0,
        &CONFIG,
        true,
    );
    println!("len_a: {len_a}, len_b: {len_b}");

    let mut tsm_entrance_compact_index = 1;
    let mut tsm_alignment = alignment.clone();
    for _ in 0..3 {
        assert!(tsm_alignment.move_template_switch_start_forwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &mut tsm_entrance_compact_index,
        ));
        println!("Alignment: {tsm_alignment}");
        let (_, len_a_new, len_b_new) = tsm_alignment.compute_cost_and_lengths(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &CONFIG,
            true,
        );
        assert_eq!((len_a, len_b), (len_a_new, len_b_new));
    }

    let mut tsm_entrance_compact_index = 1;
    let mut tsm_alignment = alignment.clone();
    for _ in 0..3 {
        assert!(tsm_alignment.move_template_switch_start_backwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &mut tsm_entrance_compact_index,
        ));
        println!("Alignment: {tsm_alignment}");
        let (_, len_a_new, len_b_new) = tsm_alignment.compute_cost_and_lengths(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &CONFIG,
            true,
        );
        assert_eq!((len_a, len_b), (len_a_new, len_b_new));
    }

    let tsm_entrance_compact_index = 1;
    let mut tsm_alignment = alignment.clone();
    for _ in 0..3 {
        assert!(tsm_alignment.move_template_switch_end_forwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            tsm_entrance_compact_index,
        ));
        println!("Alignment: {tsm_alignment}");
        let (_, len_a_new, len_b_new) = tsm_alignment.compute_cost_and_lengths(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &CONFIG,
            true,
        );
        assert_eq!((len_a, len_b), (len_a_new, len_b_new));
    }

    let tsm_entrance_compact_index = 1;
    let mut tsm_alignment = alignment.clone();
    for _ in 0..3 {
        assert!(tsm_alignment.move_template_switch_end_backwards(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            tsm_entrance_compact_index,
        ));
        println!("Alignment: {tsm_alignment}");
        let (_, len_a_new, len_b_new) = tsm_alignment.compute_cost_and_lengths(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            0,
            0,
            &CONFIG,
            false,
        );
        assert_eq!((len_a, len_b), (len_a_new, len_b_new));
    }
}
