#![allow(rustdoc::redundant_explicit_links)]

use generic_a_star::cost::U32Cost;
use lib_tsalign::{alignment_geometry::AlignmentRange, alignment_result::AlignmentResult};
use log::{debug, info, trace};

use crate::{
    alignment::{coordinates::PrimaryAlignmentCoordinates, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::{
        limits::AlignmentLimits, performance_parameters::AlignmentPerformanceParameters,
    },
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::AlignmentCosts,
};

pub(crate) mod alignment;
pub(crate) mod alignment_history;
pub(crate) mod anchors;
pub mod chain_align;
pub(crate) mod chaining_cost_function;
pub mod chaining_lower_bounds;
pub mod config;
pub mod costs;
pub(crate) mod gap_fill;
mod panic_on_extend;
#[cfg(test)]
mod tests;

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
/// * `anchor_k` is the size of the anchors.
/// * `max_anchor_mutations` is the maximum number of mutations allowed in an anchor.
/// * `alignment_costs` is the cost function for the alignment.
///
/// Panics if the `config` is not supported, see [`ChainingLowerBoundConfig::validate`].
/// Validate it beforehand to report unsupported parameters as an error.
pub fn preprocess(
    max_n: usize,
    config: ChainingLowerBoundConfig,
    alignment_costs: AlignmentCosts<U32Cost>,
) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(max_n, config, alignment_costs)
}

/// Align two sequences without limits.
///
/// See [`align_with_limits`] for the parameters.
#[expect(clippy::too_many_arguments)]
pub fn align(
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    performance_parameters: &AlignmentPerformanceParameters<U32Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    reference_name: &str,
    query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<U32Cost> {
    align_with_limits(
        reference,
        query,
        range,
        performance_parameters,
        &AlignmentLimits::none(),
        rc_fn,
        reference_name,
        query_name,
        chaining_lower_bounds,
    )
}

/// Align two sequences.
///
/// Note that `reference` and `query` are interchangeable, and the order has no meaning.
///
/// * `reference` is the reference string in ASCII format. Only characters `A`, `C`, `G` and `T` are allowed.
/// * `query` is the query string in ASCII format. Only characters `A`, `C`, `G` and `T` are allowed.
/// * `range` is the range on which the alignment happens. Note that points 2 and 3 of a template switch may fall outside of this range.
/// * `performance_parameters` is a set of parameters for the aligner that only affect performance.
///   They do not affect the cost of the optimal chain.
/// * `limits` are a cost limit and an iteration limit, see [`AlignmentLimits`].
///   If a limit is reached, then the result contains no alignment, but reports the limit.
/// * `rc_fn` is a function that maps a character to its reverse complement.
/// * `reference_name` is the name of the reference string. It is irrelevant for the alignment, but will appear in e.g. the output of `tsalign show`.
/// * `query_name` is the name of the query string. It is irrelevant for the alignment, but will appear in e.g. the output of `tsalign show`.
/// * `chaining_lower_bounds` are the lower bounds computed with the function [`preprocess`].
///   Their config also determines the anchors, i.e. the anchor size `k` and the maximum number of mutations allowed in an anchor.
#[expect(clippy::too_many_arguments)]
pub fn align_with_limits(
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    performance_parameters: &AlignmentPerformanceParameters<U32Cost>,
    limits: &AlignmentLimits<U32Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    reference_name: &str,
    query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<U32Cost> {
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
        PrimaryAlignmentCoordinates::new(range.reference_offset(), range.query_offset()),
        PrimaryAlignmentCoordinates::new(range.reference_limit(), range.query_limit()),
    );
    // The anchors, the gap alignments and the lower bounds must agree on these parameters,
    // so they are taken from the config of the lower bounds only.
    let k = chaining_lower_bounds.config().anchor_k;
    let max_mutations = chaining_lower_bounds.config().max_anchor_mutations;

    let anchors = if max_mutations == 0 {
        Anchors::new_exact(&sequences, k, rc_fn)
    } else {
        let mut anchors = Anchors::new_inexact(
            &sequences,
            k,
            max_mutations,
            chaining_lower_bounds.alignment_costs(),
            rc_fn,
        );
        // Pruning dominated anchors is only proven to be sound for the max-match-run chaining.
        if performance_parameters.prune_dominated_anchors
            && matches!(
                chaining_lower_bounds.config().inexact_lower_bound,
                InexactLowerBoundKind::MaxMatchRun
            )
        {
            anchors.prune_dominated_primary_anchors(
                &sequences,
                k,
                max_mutations,
                chaining_lower_bounds.alignment_costs(),
            );
        }
        anchors
    };
    trace!("Anchors:\n{anchors}");
    let mut chaining_cost_function = ChainingCostFunction::new_from_lower_bounds(
        chaining_lower_bounds,
        &anchors,
        &sequences,
        performance_parameters.max_exact_cost_function_cost,
        rc_fn,
    );

    chain_align::align(
        &sequences,
        performance_parameters,
        limits,
        chaining_lower_bounds.alignment_costs(),
        rc_fn,
        chaining_lower_bounds.config(),
        &anchors,
        &mut chaining_cost_function,
    )
}
