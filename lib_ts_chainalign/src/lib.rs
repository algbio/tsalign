#![allow(rustdoc::redundant_explicit_links)]

use generic_a_star::{arena::ArenaIndex, cost::U32Cost};
use std::time::Instant;

use clap::ValueEnum;
use lib_tsalign::{
    alignment_geometry::{AlignmentRange, InvalidAlignmentRange},
    alignment_result::{AlignmentResult, ChainStatistics},
    arena_index::ArenaIndexType,
    error::SequenceRole,
};
use log::{debug, info, trace};
use noisy_float::types::R64;
use thiserror::Error;

use crate::{
    alignment::{coordinates::PrimaryAlignmentCoordinates, sequences::AlignmentSequences},
    anchors::{Anchors, inexact_pruning::AnchorType},
    chain_align::{
        limits::AlignmentLimits,
        performance_parameters::{AlignmentPerformanceParameters, AnchorPruning},
    },
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind, InvalidChainingLowerBoundConfig},
    costs::AlignmentCosts,
};

pub(crate) mod alignment;
pub(crate) mod alignment_history;
pub(crate) mod anchors;
pub mod cache;
pub mod chain_align;
pub(crate) mod chaining_cost_function;
pub mod chaining_lower_bounds;
pub mod config;
pub mod costs;
pub(crate) mod gap_fill;
mod panic_on_extend;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod validation_tests;

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

/// A reason why [`align`] fails: the inputs are invalid, or the chaining search failed ([`Self::Search`]).
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ChainAlignError {
    /// The config of the chaining lower bounds is not supported, see [`ChainingLowerBoundConfig::validate`].
    ///
    /// This can only happen if the lower bounds were read from a file.
    #[error("Unsupported chaining lower bound config: {0}")]
    InvalidConfig(#[from] InvalidChainingLowerBoundConfig),

    /// The alignment range does not fit the sequences.
    #[error("Invalid alignment range: {0}")]
    InvalidRange(#[from] InvalidAlignmentRange),

    /// A sequence is longer than the maximum sequence length of the chaining lower bounds.
    #[error(
        "The {sequence} has length {length}, but the chaining lower bounds were preprocessed only up to length {max_n}"
    )]
    SequenceTooLong {
        sequence: SequenceRole,
        length: usize,
        max_n: usize,
    },

    /// A sequence contains a character other than `A`, `C`, `G` and `T`.
    #[error(
        "{} contains character {character:?} at position {position}, but tschainalign supports only the characters A, C, G and T",
        .sequence.capitalised()
    )]
    NonDnaCharacter {
        sequence: SequenceRole,
        position: usize,
        character: char,
    },

    /// The chaining search failed.
    ///
    /// This happens if it needs more nodes at the same time than the arena index type allows
    /// ([`AStarError::ArenaFull`](generic_a_star::AStarError::ArenaFull)).
    /// Then the alignment can be repeated with a larger index type, see
    /// [`AlignmentPerformanceParameters::arena_index`].
    #[error(transparent)]
    Search(#[from] generic_a_star::AStarError),
}

/// Perform preprocessing for tschainalign.
///
/// * `max_n` is the maximum sequence length that the lower bounds should support.
/// * `config` contains the size of the anchors `anchor_k`, the maximum number of mutations allowed in an anchor
///   `max_anchor_mutations` and the kind of lower bound for inexact anchors.
/// * `alignment_costs` is the cost function for the alignment.
///
/// # Errors
///
/// Returns an error if the `config` is not supported, see [`ChainingLowerBoundConfig::validate`].
pub fn preprocess(
    max_n: usize,
    config: ChainingLowerBoundConfig,
    alignment_costs: AlignmentCosts<U32Cost>,
) -> Result<ChainingLowerBounds<U32Cost>, InvalidChainingLowerBoundConfig> {
    ChainingLowerBounds::try_new(max_n, config, alignment_costs)
}

/// Check that tschainalign can align the given sequences on the given range.
///
/// The sequences must consist of the characters `A`, `C`, `G` and `T` only,
/// and the range must fit the sequences, see [`AlignmentRange::validate`].
///
/// [`align`] calls this as well, so this is only needed for reporting invalid inputs before preprocessing.
///
/// # Errors
///
/// Returns [`ChainAlignError::NonDnaCharacter`] or [`ChainAlignError::InvalidRange`].
pub fn validate_sequences(
    reference: &[u8],
    query: &[u8],
    range: &AlignmentRange,
) -> Result<(), ChainAlignError> {
    for (sequence, role) in [
        (reference, SequenceRole::Reference),
        (query, SequenceRole::Query),
    ] {
        if let Some(position) = sequence
            .iter()
            .position(|c| !matches!(c, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(ChainAlignError::NonDnaCharacter {
                sequence: role,
                position,
                character: char::from(sequence[position]),
            });
        }
    }

    range.validate(reference.len(), query.len())?;
    Ok(())
}

/// Align two sequences without limits.
///
/// See [`align_with_limits`] for the parameters.
///
/// # Errors
///
/// Returns an error if the inputs are invalid, see [`ChainAlignError`] and [`validate_sequences`],
/// or [`ChainAlignError::Search`] if the chaining search needs more nodes than the arena index type allows.
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
) -> Result<AlignmentResult<U32Cost>, ChainAlignError> {
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
///   They must have been preprocessed for at least the length of the longer sequence.
///
/// # Errors
///
/// Returns an error if the inputs are invalid, see [`ChainAlignError`] and [`validate_sequences`],
/// or [`ChainAlignError::Search`] if the chaining search needs more nodes than the arena index type allows.
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
) -> Result<AlignmentResult<U32Cost>, ChainAlignError> {
    match performance_parameters.arena_index {
        ArenaIndexType::U32 => align_with_index::<u32>(
            reference,
            query,
            range,
            performance_parameters,
            limits,
            rc_fn,
            reference_name,
            query_name,
            chaining_lower_bounds,
        ),
        ArenaIndexType::U64 => align_with_index::<u64>(
            reference,
            query,
            range,
            performance_parameters,
            limits,
            rc_fn,
            reference_name,
            query_name,
            chaining_lower_bounds,
        ),
    }
}

/// [`align_with_limits`] with the node index type `I` of the chaining search,
/// ignoring [`AlignmentPerformanceParameters::arena_index`].
#[expect(clippy::too_many_arguments)]
pub(crate) fn align_with_index<I: ArenaIndex>(
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    performance_parameters: &AlignmentPerformanceParameters<U32Cost>,
    limits: &AlignmentLimits<U32Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    reference_name: &str,
    query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> Result<AlignmentResult<U32Cost>, ChainAlignError> {
    validate_sequences(&reference, &query, &range)?;
    chaining_lower_bounds.config().validate()?;
    let max_n = chaining_lower_bounds.max_n();
    for (sequence, role) in [
        (&reference, SequenceRole::Reference),
        (&query, SequenceRole::Query),
    ] {
        if sequence.len() > max_n {
            return Err(ChainAlignError::SequenceTooLong {
                sequence: role,
                length: sequence.len(),
                max_n,
            });
        }
    }

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

    let mut statistics = ChainStatistics::default();
    let anchors = generate_and_prune_anchors(
        &sequences,
        performance_parameters,
        rc_fn,
        chaining_lower_bounds,
        &mut statistics,
    );
    trace!("Anchors:\n{anchors}");

    let cost_function_start_time = Instant::now();
    let mut chaining_cost_function = ChainingCostFunction::new_from_lower_bounds(
        chaining_lower_bounds,
        &anchors,
        &sequences,
        performance_parameters.max_exact_cost_function_cost,
        rc_fn,
    );
    statistics.cost_function_seconds = seconds_since(cost_function_start_time);

    let mut result = chain_align::align::<_, I>(
        &sequences,
        performance_parameters,
        limits,
        &mut statistics,
        chaining_lower_bounds.alignment_costs(),
        rc_fn,
        chaining_lower_bounds.config(),
        &anchors,
        &mut chaining_cost_function,
    )?;
    result.statistics_mut().chain_statistics = Some(statistics);
    Ok(result)
}

fn seconds_since(start_time: Instant) -> R64 {
    start_time.elapsed().as_secs_f64().try_into().unwrap()
}

/// Generates the anchors and prunes them according to the performance parameters.
///
/// Fills the anchor counts, the pruning statistics and the times in `statistics`.
fn generate_and_prune_anchors(
    sequences: &AlignmentSequences,
    performance_parameters: &AlignmentPerformanceParameters<U32Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
    statistics: &mut ChainStatistics,
) -> Anchors<U32Cost> {
    // The anchors, the gap alignments and the lower bounds must agree on these parameters,
    // so they are taken from the config of the lower bounds only.
    let k = chaining_lower_bounds.config().anchor_k;
    let max_mutations = chaining_lower_bounds.config().max_anchor_mutations;
    let max_match_run_lower_bound = matches!(
        chaining_lower_bounds.config().inexact_lower_bound,
        InexactLowerBoundKind::MaxMatchRun
    );
    let costs = chaining_lower_bounds.alignment_costs();

    let anchor_pruning = performance_parameters.anchor_pruning;
    statistics.anchor_pruning = anchor_pruning
        .to_possible_value()
        .unwrap()
        .get_name()
        .to_string();
    if anchor_pruning == AnchorPruning::Inexact {
        statistics.inexact_pruning_rules = performance_parameters.inexact_pruning.rules.to_string();
    }

    let anchor_start_time = Instant::now();
    let mut anchors = if max_mutations == 0 {
        Anchors::new_exact(sequences, k, rc_fn)
    } else {
        Anchors::new_inexact(sequences, k, max_mutations, costs, rc_fn)
    };
    statistics.anchor_seconds = seconds_since(anchor_start_time);
    statistics.primary_anchors = anchors.primary_anchor_count() as u64;
    statistics.secondary_anchors = anchors.secondary_anchor_count() as u64;

    // Exact anchors are never pruned.
    let pruning_start_time = Instant::now();
    if max_mutations > 0 {
        match anchor_pruning {
            AnchorPruning::None => {}
            AnchorPruning::Exact => {
                // Pruning dominated anchors is only proven to be sound for the max-match-run chaining.
                if max_match_run_lower_bound {
                    let pruned =
                        anchors.prune_dominated_primary_anchors(sequences, k, max_mutations, costs);
                    statistics
                        .pruned_anchors
                        .insert("primary-exact-domination".to_string(), pruned as u64);
                    statistics
                        .matched_anchors
                        .insert("primary-exact-domination".to_string(), pruned as u64);
                }
            }
            AnchorPruning::Inexact => {
                let (primary_counts, secondary_counts) = anchors.prune_inexact(
                    sequences,
                    k,
                    max_mutations,
                    costs,
                    max_match_run_lower_bound,
                    &performance_parameters.inexact_pruning,
                    rc_fn,
                );
                primary_counts.insert_into(
                    AnchorType::Primary,
                    &mut statistics.matched_anchors,
                    &mut statistics.pruned_anchors,
                );
                secondary_counts.insert_into(
                    AnchorType::Secondary,
                    &mut statistics.matched_anchors,
                    &mut statistics.pruned_anchors,
                );
            }
        }
    }
    statistics.pruning_seconds = seconds_since(pruning_start_time);
    statistics.primary_anchors_after_pruning = anchors.primary_anchor_count() as u64;
    statistics.secondary_anchors_after_pruning = anchors.secondary_anchor_count() as u64;

    anchors
}
