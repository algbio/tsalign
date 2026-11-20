#![allow(rustdoc::redundant_explicit_links)]

use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::U32Cost;
use lib_tsalign::a_star_aligner::{
    alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    template_switch_distance::AlignmentType,
};
use log::{debug, info, trace};

use crate::{
    alignment::{coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::performance_parameters::AlignmentPerformanceParameters,
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
    costs::AlignmentCosts,
};

pub mod alignment;
pub mod anchors;
pub mod chain_align;
pub mod chaining_cost_function;
pub mod chaining_lower_bounds;
pub mod costs;
pub mod exact_chaining;
pub mod panic_on_extend;

/// A reverse complement function for DNA alphabets.
pub fn dna_rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => panic!("Unsupported character: {c}"),
    }
}

/// Perform preprocessing for tschainalign.
///
/// * `max_n` is the maximum sequence length that the lower bounds should support.
/// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
///   Set this to `k-1`, if the anchors are supposed to be `k`-mers.
/// * `alignment_costs` is the cost function for the alignment.
pub fn preprocess(
    max_n: usize,
    max_match_run: u32,
    alignment_costs: AlignmentCosts<U32Cost>,
) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(max_n, max_match_run, alignment_costs)
}

/// Align two sequences.
///
/// Note that `reference` and `query` are interchangeable, and the order has no meaning.
///
/// * `AlphabetType` must be a DNA alphabet.
///
/// * `reference` is the reference string in ASCII format. Only characters `A`, `C`, `G` and `T` are allowed.
/// * `query` is the query string in ASCII format. Only characters `A`, `C`, `G` and `T` are allowed.
/// * `range` is the range on which the alignment happens. Note that points 2 and 3 of a template switch may fall outside of this range.
/// * `performance_parameters` is a set of parameters for the aligner that only affect performance.
/// * `rc_fn` is a function that maps a character to its reverse complement.
/// * `reference_name` is the name of the reference string. It is irrelevant for the alignment, but will appear in e.g. the output of `tsalign show`.
/// * `query_name` is the name of the query string. It is irrelevant for the alignment, but will appear in e.g. the output of `tsalign show`.
/// * `chaining_lower_bounds` are the lower bounds computed with the function [`preprocess`].
#[expect(clippy::too_many_arguments)]
pub fn align<AlphabetType: Alphabet>(
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    performance_parameters: &AlignmentPerformanceParameters<U32Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    reference_name: &str,
    query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<AlignmentType, U32Cost> {
    debug!(
        "Reference sequence: {}",
        String::from_utf8_lossy(&reference)
    );
    debug!("Query sequence: {}", String::from_utf8_lossy(&query));
    info!("Aligning on subsequence {}", range);

    let sequences = AlignmentSequences::new_named(
        reference,
        query,
        reference_name.to_string(),
        query_name.to_string(),
        AlignmentCoordinates::new_primary(range.reference_offset(), range.query_offset()),
        AlignmentCoordinates::new_primary(range.reference_limit(), range.query_limit()),
    );
    let k = chaining_lower_bounds.max_match_run() + 1;

    let anchors = Anchors::new(&sequences, k, rc_fn);
    trace!("Anchors:\n{anchors}");
    let mut chaining_cost_function = ChainingCostFunction::new_from_lower_bounds(
        chaining_lower_bounds,
        &anchors,
        &sequences,
        performance_parameters.max_exact_cost_function_cost,
        rc_fn,
    );

    chain_align::align::<AlphabetType, _>(
        &sequences,
        performance_parameters,
        chaining_lower_bounds.alignment_costs(),
        rc_fn,
        chaining_lower_bounds.max_match_run(),
        &anchors,
        &mut chaining_cost_function,
    )
}

#[cfg(test)]
mod tests {
    use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;

    use crate::costs::{GapAffineCosts, TsLimits};

    use super::*;

    fn align_on_range(reference: &[u8], query: &[u8], range: AlignmentRange) -> (U32Cost, String) {
        let costs = AlignmentCosts {
            primary_costs: GapAffineCosts::new(2u32.into(), 3u32.into(), 1u32.into()),
            secondary_costs: GapAffineCosts::new(4u32.into(), 3u32.into(), 2u32.into()),
            ts_base_cost: U32Cost::from(2u32).into(),
            ts_limits: TsLimits {
                length_23: 8..401,
                ..TsLimits::new_unlimited()
            },
        };
        let k = 5;
        let max_n = reference.len().max(query.len()).next_power_of_two();
        let lower_bounds = preprocess(max_n, k - 1, costs);

        let result = align::<DnaAlphabet>(
            reference.to_vec(),
            query.to_vec(),
            range,
            &AlignmentPerformanceParameters::default(),
            &dna_rc_fn,
            "reference",
            "query",
            &lower_bounds,
        );
        (result.statistics().result.cost(), result.cigar())
    }

    /// The template switch starts and ends close to the range edges, matches in between.
    #[test]
    fn match_run_between_range_edge_and_template_switch() {
        let reference = b"CAAATCCCACACTTGTGTGAGGATGTACACGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGCATGGCCAG";
        let query = b"CAAATCCCACACTTGTGTGAGGATGTACACGTGTGTGTGTGTGTGTGTGTGCAAGGCCAG";

        let (cost, cigar) = align_on_range(reference, query, AlignmentRange::new(25, 25, 47, 35));
        assert_eq!(cost, U32Cost::from(2u32), "{cigar}");
    }

    /// The template of the template switch lies outside the range
    #[test]
    fn template_beyond_range_end() {
        let reference =
            b"GATTACAGCTTGACCTAGGCATCGTTAGCAATGCCGTAAGCTTCGAGGTACCTTAGGCTAACGTTGCAGTCCATGGAATC";
        // `reference` with 25..35 replaced by the reverse complement of 60..70.
        let query =
            b"GATTACAGCTTGACCTAGGCATCGTACTGCAACGTGTAAGCTTCGAGGTACCTTAGGCTAACGTTGCAGTCCATGGAATC";

        let (cost, cigar) = align_on_range(reference, query, AlignmentRange::new(20, 20, 40, 40));
        assert_eq!(cost, U32Cost::from(2u32), "{cigar}");
    }
}
