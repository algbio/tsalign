use std::fs::File;

use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U64Cost;
use noisy_float::types::r64;

use crate::{
    a_star_aligner::{
        alignment_geometry::{AlignmentCoordinates, AlignmentRange},
        template_switch_distance::strategies::{
            AlignmentStrategySelection,
            chaining::NoChainingStrategy,
            node_ord::AntiDiagonalNodeOrdStrategy,
            primary_match::AllowPrimaryMatchStrategy,
            primary_range::NoPrunePrimaryRangeStrategy,
            secondary_deletion::AllowSecondaryDeletionStrategy,
            shortcut::NoShortcutStrategy,
            template_switch_count::NoTemplateSwitchCountStrategy,
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                PreprocessedLookaheadTemplateSwitchMinLengthStrategy,
                PreprocessedTemplateSwitchMinLengthStrategy,
            },
            template_switch_total_length::NoTemplateSwitchTotalLengthStrategy,
        },
        template_switch_distance_a_star_align,
    },
    config::TemplateSwitchConfig,
};

#[test]
fn test_tsnax_disc1_473() {
    let config = TemplateSwitchConfig::<DnaAlphabet, U64Cost>::read_plain(
        File::open("../test_files/config/bench/config.tsa").unwrap(),
    )
    .unwrap();

    let reference = VectorGenome::<DnaAlphabet>::from_slice_u8( b"TATTTGTTGTTAGCAGATTAAAGATTAGCTAGACGAAATTCTCTGGGAGGTACAGAAAGGAATATAGAATAAAAAGAGTAGAAACACAGTAGAATATCCTTTACCTCATGTGTTTTATGTATTGTATTATGGTGATTGAGACAAAAATTTTAATGTCATCTGATACTTAAGACAATAGTATTTTATTTTATTTTATTTTATTTTTGAGACAGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGTGATCTTGGTTCACCGCAACCTCCAATTCCCAGGTTCAAGTGATTCTCCTGTCTCAGCCTCCTGAGTAGCTGGAACTACAAGCATGTGCCACTGTGTCCAGCTAATTTTTGCATTTTTAGTAGAGACGGGGTTTCACCATATTGGTCAGGCTGGTCTCCAACTCCTG").unwrap();
    let query = VectorGenome::from_slice_u8( b"TATTTGTTGTTAGCAGATTAAAGATTAGCTAGATGAAATTCTCTGGGAGGTACAGAAAGGAATATAGAATAAAAAGAGTAGAAACACAGTAGAATATCCTTTACCTCATGTGTTTTATGTATTGTATTATGGTGATTGAGACAAAAATTTTAATGTCATCTGATACTTAAGACAATAGTATTTTATTTTATTTTATTTTTGAGATAGAGTCTCACTCTCACCCAGGCTGGAGTGCAGTGGCGTGATCTTGGTTCACCGCAACCTCCAATTCCCAGGTTCAAGTGATTCCCCTGTCTCAGCCTCCTGAGTAGCTGGAACTACAAGCATGTGCCACTGTGTCCAGCTAATTTTTGCATTTTTAGTAGAGATGGGGTTTCACCATATTGGTCAGGCTGGTCTCCAACTCCTG").unwrap();
    // R196..219Q196..212
    let range = AlignmentRange::new_offset_limit(
        AlignmentCoordinates::new(196, 196),
        AlignmentCoordinates::new(219, 212),
    );

    let result = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            DnaAlphabet,
            U64Cost,
            AntiDiagonalNodeOrdStrategy,
            NoTemplateSwitchMinLengthStrategy<U64Cost>,
            NoChainingStrategy<U64Cost>,
            NoTemplateSwitchCountStrategy,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            NoTemplateSwitchTotalLengthStrategy,
        >,
        _,
    >(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        "ref",
        "qry",
        range.clone(),
        &config,
        None,
        None,
        false,
        true,
        (),
    );
    assert_eq!(
        result.statistics().reference_offset,
        result.statistics().query_offset
    );
    // Sample alignment with offset 34: 165M[TSQQR:[0,0]:[0,0]:21:5M1D1M1I3M:17]2M1S1M
    let sample_alignment = format!(
        "Sample alignment with offset {}: {}",
        result.statistics().reference_offset,
        result.cigar()
    );
    println!("{sample_alignment}");
    assert_eq!(result.statistics().cost, r64(10.0));

    let result = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            DnaAlphabet,
            U64Cost,
            AntiDiagonalNodeOrdStrategy,
            LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            NoChainingStrategy<U64Cost>,
            NoTemplateSwitchCountStrategy,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            NoTemplateSwitchTotalLengthStrategy,
        >,
        _,
    >(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        "ref",
        "qry",
        range.clone(),
        &config,
        None,
        None,
        false,
        true,
        (),
    );
    println!("{sample_alignment}");
    assert_eq!(result.statistics().cost, r64(10.0));

    let result = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            DnaAlphabet,
            U64Cost,
            AntiDiagonalNodeOrdStrategy,
            PreprocessedTemplateSwitchMinLengthStrategy<false, U64Cost>,
            NoChainingStrategy<U64Cost>,
            NoTemplateSwitchCountStrategy,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            NoTemplateSwitchTotalLengthStrategy,
        >,
        _,
    >(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        "ref",
        "qry",
        range.clone(),
        &config,
        None,
        None,
        false,
        true,
        (),
    );
    println!("{sample_alignment}");
    assert_eq!(result.statistics().cost, r64(10.0));

    let result = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            DnaAlphabet,
            U64Cost,
            AntiDiagonalNodeOrdStrategy,
            PreprocessedLookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            NoChainingStrategy<U64Cost>,
            NoTemplateSwitchCountStrategy,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            NoTemplateSwitchTotalLengthStrategy,
        >,
        _,
    >(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        "ref",
        "qry",
        range.clone(),
        &config,
        None,
        None,
        false,
        true,
        (),
    );
    println!("{sample_alignment}");
    assert_eq!(result.statistics().cost, r64(10.0));
}
