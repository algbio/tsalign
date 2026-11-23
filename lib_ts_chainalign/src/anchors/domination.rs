//! Pruning of dominated inexact primary anchors.
//!
//! Inexact anchors create many chains of equal cost, which slows down the chaining.
//! Some of the anchors are redundant: every chain that uses them can be replaced by a chain that does not use them
//! and does not cost more.
//! Removing such anchors before the chaining does not change the cost of the optimal chain.
//!
//! # The pruning rule
//!
//! A primary anchor `A2` is *dominated* if
//!
//! * the anchors are generated with exactly one allowed mutation (`max_anchor_mutations == 1`),
//! * one of the sides of `A2` has length `k + 1`,
//! * `A2` has a positive cost, and
//! * there is a primary anchor `A1` with the same start as `A2`, such that `A2` is `A1` extended by a single alignment
//!   column `c` (a match, a substitution or a gap of length one) at its end, and `cost(A1) + cost(c) <= cost(A2)`.
//!
//! In other words, `A2` is `A1` plus one column that the gap alignment after `A2` could align just as cheaply.
//! With one allowed mutation, an anchor with a side of length `k + 1` consists of `k` matches and one gap character,
//! and removing its last column results in an exact anchor (if the last column is the gap)
//! or in an anchor with one side of length `k - 1` (if the last column is a match).
//! Hence, if the anchors are complete, all primary anchors with a side of length `k + 1` are dominated.
//!
//! # Soundness
//!
//! The rule is restricted to the case in which we can prove that pruning does not increase the cost of the optimal chain:
//!
//! * Only for `max_anchor_mutations == 1`.
//!   For more mutations, the alignment after the shortened anchor may violate the limit on the length of match runs,
//!   and there is not always an anchor that allows to repair the chain.
//! * Only for the max-match-run chaining, since the proof uses its rules for the alignments between anchors.
//! * Only if extending a gap is not more expensive than opening a new gap,
//!   since the removed column may be merged with a gap at the start of the following gap alignment.
//! * Only primary anchors.
//!   The length of a template switch is measured between the start of its first and the end of its last secondary anchor.
//!   Replacing a secondary anchor by a shorter one, or adding a further secondary anchor to a template switch, changes
//!   this length, and may move it out of the range allowed by the configuration.
//!   Without any restriction, pruning dominated secondary anchors increases the cost of
//!   `test_files/twin_in_place_inversion.fa` with `k = 9` from 5 to 18.
//!
//! Dominators are never dominated themselves: they have no side of length `k + 1`.
//! Anchors with zero cost are never dominated.
//!
//! The proof is an exchange argument: in a chain, `A2` is replaced by `A1`, and the removed column is prepended to the
//! alignment after the anchor.
//! If this alignment then violates the limit on the length of match runs, it contains a further anchor that is not
//! dominated, and the chain is routed through that anchor.

use generic_a_star::cost::AStarCost;
use rustc_hash::FxHashMap;

use crate::{
    alignment::coordinates::range::PrimaryAlignmentRange, anchors::primary::PrimaryAnchor,
    costs::GapAffineCosts,
};

#[cfg(test)]
mod tests;

/// Returns true if pruning dominated primary anchors is known to not change the cost of the optimal chain.
///
/// See the [module documentation](self) for the conditions.
/// The caller must additionally ensure that the max-match-run chaining is used.
pub fn is_pruning_sound<Cost: AStarCost>(
    max_mutations: usize,
    primary_costs: &GapAffineCosts<Cost>,
) -> bool {
    max_mutations == 1 && primary_costs.gap_extend <= primary_costs.gap_open
}

/// Removes all dominated primary anchors from the given anchors and returns the number of removed anchors.
///
/// See the [module documentation](self) for the definition of dominated anchors.
/// This must only be called if [`is_pruning_sound`] returns true.
///
/// The anchors may be in any order, and their order is preserved.
/// `seq1` and `seq2` are the complete sequences that the coordinates of the anchors refer to.
pub fn prune_dominated_primary_anchors<Cost: AStarCost>(
    anchors: &mut Vec<PrimaryAnchor<Cost>>,
    seq1: &[u8],
    seq2: &[u8],
    k: usize,
    primary_costs: &GapAffineCosts<Cost>,
) -> usize {
    debug_assert!(primary_costs.gap_extend <= primary_costs.gap_open);

    let mut cost_by_range: FxHashMap<PrimaryAlignmentRange, Cost> = FxHashMap::default();
    for anchor in anchors.iter() {
        // The anchors contain only one anchor per range, but in case of duplicates, the cheapest anchor counts.
        cost_by_range
            .entry(anchor.range())
            .and_modify(|cost| *cost = (*cost).min(anchor.cost()))
            .or_insert_with(|| anchor.cost());
    }

    let len_before = anchors.len();
    anchors.retain(|anchor| !is_dominated(anchor, &cost_by_range, seq1, seq2, k, primary_costs));
    len_before - anchors.len()
}

/// Returns true if the given anchor is dominated by an anchor in `cost_by_range`.
fn is_dominated<Cost: AStarCost>(
    anchor: &PrimaryAnchor<Cost>,
    cost_by_range: &FxHashMap<PrimaryAlignmentRange, Cost>,
    seq1: &[u8],
    seq2: &[u8],
    k: usize,
    primary_costs: &GapAffineCosts<Cost>,
) -> bool {
    let (start_a, start_b) = (anchor.start().a(), anchor.start().b());
    let (end_a, end_b) = (anchor.end().a(), anchor.end().b());
    let (len_a, len_b) = (end_a - start_a, end_b - start_b);

    // Only anchors with a side of length k + 1 are candidates.
    if len_a != k + 1 && len_b != k + 1 {
        return false;
    }
    // Anchors with zero cost are never pruned, since they may be chained for free with an overlap,
    // which the proof does not cover (see `PrimaryAnchor::is_direct_free_predecessor_of`).
    // Such an anchor is only a candidate if gaps are free.
    if anchor.cost().is_zero() {
        return false;
    }

    let last_both_cost = if seq1[end_a - 1] == seq2[end_b - 1] {
        Cost::zero()
    } else {
        primary_costs.substitution
    };

    // The dominator is the anchor without its last column, which is one of:
    // a gap in the second sequence, a gap in the first sequence, or a match or substitution.
    [
        (1, 0, primary_costs.gap_open),
        (0, 1, primary_costs.gap_open),
        (1, 1, last_both_cost),
    ]
    .into_iter()
    .any(|(delta_a, delta_b, column_cost)| {
        let (dominator_len_a, dominator_len_b) = (len_a - delta_a, len_b - delta_b);
        // A dominator must not be a candidate itself, such that it is never pruned.
        if dominator_len_a > k || dominator_len_b > k {
            return false;
        }

        let dominator_range = PrimaryAlignmentRange::new_from_ranges(
            start_a..end_a - delta_a,
            start_b..end_b - delta_b,
        );
        cost_by_range
            .get(&dominator_range)
            .and_then(|dominator_cost| dominator_cost.checked_add(&column_cost))
            .is_some_and(|dominated_cost| dominated_cost <= anchor.cost())
    })
}
