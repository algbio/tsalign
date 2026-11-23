use std::hash::Hash;

use crate::alignment_history::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;
pub mod unsigned_int_alignment_history_vec;

/// Keeps track of the most recent alignments required to reach a DP node.
///
/// Only keeps track of enough state to prove that a DP node does not end on a possibly valid anchor,
/// i.e. that when backtracking the node, at least `max_anchor_mutations + 1` mutations are required to create an anchor of lenght `anchor_k`.
///
/// # Properties
///
/// The alignment history can be described by a tuple of two flags.
/// First, `is_anchor_length` is true if and only if the history covers at least `anchor_k` bases for one of the two sequences.
/// Second, `is_invalid` is true if and only if the history contains more than `max_anchor_mutations` mutations.
///
/// When extending the alignment, the history develops in two phases: the initial phase and the general phase.
/// There is one transition from initial phase to general phase at some point, and no other phase transition.
/// The initial phase is when the history is too short to describe a full anchor, but there are not enough mutations to invalidate an anchor.
/// Formally, in the initial phase, `is_anchor_length` and `is_invalid` are both false.
/// As soon as an extension causes one of them to become true, the anchor has transitioned to the general phase.
/// Once in the general phase, `is_invalid` must always be true, otherwise an extension would result in a valid anchor.
/// However, it suffices to store as little history as necessary to keep `is_invalid` true, hence in the general phase `is_anchor_length` may be true or false.
pub trait AlignmentHistory: Default + Eq + Hash + Ord + Copy {
    /// Returns true if the alignment history has enough capacity for the given parameters.
    fn has_enough_capacity(anchor_k: u8, max_anchor_mutations: u8) -> bool;

    /// Returns the length of the alignment history vector.
    fn len(&self) -> usize;

    /// Returns true if the alignment history vector is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the maximum length that the history can have given the parameters.
    fn required_capacity(anchor_k: u8, max_anchor_mutations: u8) -> usize {
        let anchor_k = <usize as From<u8>>::from(anchor_k);
        let max_anchor_mutations = <usize as From<u8>>::from(max_anchor_mutations);

        // The required length of the initial phase is lower than that of the general phase, so we can ignore it.
        // Specifically, an initial history transitions to a general history by adding an alignment operation without removing any at the same time.
        // If any would be removed, then the initial phase history would have had enough mismatches to be general phase, or enough alignments to form a valid anchor.

        // Required len of the general phase.
        // In the general phase, there need to be `max_anchor_mutations + 1` mismatches, and each sequence has at most `anchor_k` characters in the history.
        // So, the longest history we can construct starts with `max_anchor_mutations + 1` gaps which take a length of at least `(max_anchor_mutations + 1).div_ceil(2)` for one of the sequences.
        // Then, it is followed by the maximum possible amount of matches, which is `anchor_k - (max_anchor_mutations + 1).div_ceil(2)`.
        // Summing up, we get a history of length `max_anchor_mutations + 1 + anchor_k - (max_anchor_mutations + 1).div_ceil(2)`.
        // We can rewrite that as follows:
        anchor_k + max_anchor_mutations.div_ceil(2)
    }

    /// If possible, returns the extension of this the alignment history vector with another alignment history operation.
    ///
    /// If the extension would result in a valid anchor, then None is returned.
    fn try_extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Option<Self>;
}
