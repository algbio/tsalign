//! A compact summary of an alignment history that determines which continuations of an alignment are valid.
//!
//! # Formalisation
//!
//! An alignment is a word over the [`AlignmentHistoryOperation`]s.
//! A subalignment (i.e. a contiguous subword) is *anchor-like* if it contains at most `max_anchor_mutations` mutations
//! and covers at least `anchor_k` characters of at least one of the two sequences.
//! An alignment between two chained anchors is *valid* if it contains no anchor-like subalignment.
//!
//! For an alignment history `H` and `0 <= j <= max_anchor_mutations`, let `L_j` be the longest suffix of `H` with at most `j` mutations.
//! The profile of `H` consists of the two vectors `len_a(j) = |L_j|_a` and `len_b(j) = |L_j|_b`,
//! where `|L_j|_a` and `|L_j|_b` are the number of characters of the respective sequence covered by `L_j`.
//!
//! 1. **Validity:** a valid history `H` extended by an operation stays valid if and only if `len_a(max_anchor_mutations) < anchor_k`
//!    and `len_b(max_anchor_mutations) < anchor_k` hold for the profile of the extension.
//!    This is because any new anchor-like subalignment is a suffix, and `L_max_anchor_mutations` is the longest suffix with few enough mutations.
//! 2. **Extension:** extending by a match increments all entries of the profile.
//!    Extending by a mutation that covers `e_a` characters of sequence a and `e_b` characters of sequence b
//!    sets `len_a'(0) = 0` and `len_a'(j) = len_a(j - 1) + e_a` for `j >= 1`, and analogously for `len_b`.
//! 3. **Domination:** both extension maps are monotone with respect to the componentwise order.
//!    Hence, if the profile of `A` is componentwise at most the profile of `B` and `B` extended by some continuation `C` is valid,
//!    then `A` extended by `C` is valid as well, and its profile is again at most the profile of `B` extended by `C`.
//!    In this case we say that `A` gives at least as much freedom to future alignments as `B`.
//!
//! The profile only depends on `L_max_anchor_mutations`, and there is a bijection between `L_max_anchor_mutations` and the profile.
//! Hence, it suffices to store `L_max_anchor_mutations` as the alignment history.

use crate::alignment_history::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;

/// The profile of an alignment history, see the [module documentation](self).
///
/// The entries `len_a(j)` and `len_b(j)` are stored in byte `j` of `len_a` and `len_b` respectively.
/// Unused bytes are zero.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct HistoryProfile {
    len_a: u64,
    len_b: u64,
}

/// The highest bit of each byte.
const HIGH_BITS: u64 = 0x8080_8080_8080_8080;

impl HistoryProfile {
    /// The maximum value of `max_anchor_mutations` that is supported by the packed representation.
    pub const MAX_SUPPORTED_ANCHOR_MUTATIONS: u8 = 7;

    /// The maximum length of a suffix that is supported by the packed representation.
    pub const MAX_SUPPORTED_LEN: usize = 127;

    /// Computes the profile of the given alignment history, which is given from newest to oldest operation.
    ///
    /// Operations after the `max_anchor_mutations + 1`th mutation are ignored.
    ///
    /// Returns `None` if `max_anchor_mutations` is larger than [`Self::MAX_SUPPORTED_ANCHOR_MUTATIONS`].
    ///
    /// # Panics
    ///
    /// Panics if a suffix covers more than [`Self::MAX_SUPPORTED_LEN`] characters of a sequence.
    pub fn new(
        history_newest_first: impl IntoIterator<Item = AlignmentHistoryOperation>,
        max_anchor_mutations: u8,
    ) -> Option<Self> {
        if max_anchor_mutations > Self::MAX_SUPPORTED_ANCHOR_MUTATIONS {
            return None;
        }
        let max_anchor_mutations = usize::from(max_anchor_mutations);

        let mut values_a = [0u8; 8];
        let mut values_b = [0u8; 8];
        let mut len_a = 0;
        let mut len_b = 0;
        let mut mutations = 0;

        for alignment in history_newest_first {
            if alignment.mutations() > 0 {
                // `L_mutations` ends right before this operation.
                values_a[mutations] = len_a;
                values_b[mutations] = len_b;
                mutations += 1;

                if mutations > max_anchor_mutations {
                    break;
                }
            }

            len_a += u8::try_from(alignment.len_a()).unwrap();
            len_b += u8::try_from(alignment.len_b()).unwrap();
            assert!(usize::from(len_a.max(len_b)) <= Self::MAX_SUPPORTED_LEN);
        }

        // The remaining suffixes contain the whole history.
        for mutations in mutations..=max_anchor_mutations {
            values_a[mutations] = len_a;
            values_b[mutations] = len_b;
        }

        Some(Self {
            len_a: u64::from_le_bytes(values_a),
            len_b: u64::from_le_bytes(values_b),
        })
    }

    /// The number of characters of sequence a covered by the longest suffix with at most `mutations` mutations.
    pub fn len_a(&self, mutations: u8) -> u8 {
        self.len_a.to_le_bytes()[usize::from(mutations)]
    }

    /// The number of characters of sequence b covered by the longest suffix with at most `mutations` mutations.
    pub fn len_b(&self, mutations: u8) -> u8 {
        self.len_b.to_le_bytes()[usize::from(mutations)]
    }

    /// Returns true if this profile is componentwise at most `other`.
    ///
    /// This means that an alignment with this profile gives at least as much freedom to future alignments as an alignment with the `other` profile.
    pub fn is_at_most(&self, other: &Self) -> bool {
        bytewise_at_most(self.len_a, other.len_a) && bytewise_at_most(self.len_b, other.len_b)
    }
}

/// Returns true if each byte of `x` is at most the corresponding byte of `y`.
///
/// All bytes must be smaller than 128.
fn bytewise_at_most(x: u64, y: u64) -> bool {
    debug_assert_eq!(x & HIGH_BITS, 0);
    debug_assert_eq!(y & HIGH_BITS, 0);

    // Each byte of `y | HIGH_BITS` is at least 128 and each byte of `x` is at most 127, so the subtraction never borrows across bytes.
    // The highest bit of a byte of the difference stays set if and only if the byte of `y` is at least the byte of `x`.
    ((y | HIGH_BITS) - x) & HIGH_BITS == HIGH_BITS
}
