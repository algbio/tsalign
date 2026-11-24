/// Limits for tschainalign.
///
/// Unlike the [`AlignmentPerformanceParameters`](super::performance_parameters::AlignmentPerformanceParameters),
/// the limits may change the result: if a limit is reached, then no alignment is returned.
/// Use the `default()` method for no limits.
///
/// tschainalign computes the optimal chain in rounds.
/// Each round runs an A* search for the cheapest chain with regard to the current chaining cost function,
/// whose costs are lower bounds of the exact costs of the gaps between anchors.
/// Then the chain is evaluated: its gaps are aligned exactly, and their exact costs replace the lower bounds
/// in the chaining cost function.
/// If the evaluated cost equals the cost found by the search, then the chain is optimal and the alignment is built from it.
/// Otherwise, the next round starts.
/// The cost of the chain found in a round is a lower bound of the cost of the optimal chain,
/// and it does not decrease from round to round.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignmentLimits<Cost> {
    /// An inclusive limit on the chain cost.
    ///
    /// The chain cost is the sum of the costs of the gaps and anchors of the chain, which tschainalign minimises.
    /// If the optimal chain costs at most this limit, then the alignment is returned as without a limit.
    /// Otherwise, the result is [`AStarResult::ExceededCostLimit`](generic_a_star::AStarResult::ExceededCostLimit)
    /// without an alignment.
    ///
    /// The limit is applied to the chaining searches of all rounds: they do not explore partial chains whose cost
    /// plus lower bound exceeds the limit, and stop as soon as the cheapest remaining partial chain exceeds it.
    /// The gap alignments of the chain evaluation are not limited.
    ///
    /// For inexact anchors, the cost of the final alignment may be lower than the chain cost
    /// (see [issue #221](https://github.com/algbio/tsalign/issues/221)).
    /// The limit applies to the chain cost, so an alignment whose cost is at most the limit
    /// is not returned if its chain cost exceeds the limit.
    pub cost_limit: Option<Cost>,

    /// A limit on the number of chaining rounds.
    ///
    /// If the optimal chain is not confirmed within this many rounds, then the alignment is aborted,
    /// and [`AlignmentStatistics::exceeded_iteration_limit`](lib_tsalign::alignment_result::AlignmentStatistics::exceeded_iteration_limit)
    /// is set.
    /// A limit of zero aborts before the first round.
    pub iteration_limit: Option<u64>,
}

impl<Cost> AlignmentLimits<Cost> {
    /// No limits.
    pub fn none() -> Self {
        Self {
            cost_limit: None,
            iteration_limit: None,
        }
    }
}

#[cfg(test)]
mod tests;
