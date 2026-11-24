//! A compact summary of an alignment history that determines which continuations of an alignment are anchor-free.
//!
//! # Formalisation
//!
//! An alignment is a word over the [`AlignmentHistoryOperation`]s.
//! A subalignment (i.e. a contiguous subword) is *anchor-like* if it contains at most `max_anchor_mutations` mutations
//! and covers at least `anchor_k` characters of at least one of the two sequences.
//! Gap alignments between inexact anchors must not contain any anchor-like subalignment.
//!
//! For an alignment history `H` and `0 <= j <= max_anchor_mutations`, let `L_j` be the longest suffix of `H` with at most `j` mutations.
//! The *profile* of `H` consists of the two vectors `len_a(j) = |L_j|_a` and `len_b(j) = |L_j|_b`,
//! where `|L_j|_a` and `|L_j|_b` are the numbers of characters of the respective sequence covered by `L_j`.
//!
//! 1. **Validity:** an anchor-free history `H` extended by an operation stays anchor-free if and only if
//!    `len_a(max_anchor_mutations) < anchor_k` and `len_b(max_anchor_mutations) < anchor_k` hold for the profile of the extension.
//!    This is because any new anchor-like subalignment is a suffix,
//!    and `L_max_anchor_mutations` is the longest suffix with few enough mutations.
//! 2. **Extension:** extending by a match increments all entries of the profile.
//!    Extending by a mutation that covers `e_a` characters of sequence a and `e_b` characters of sequence b
//!    sets `len_a'(0) = 0` and `len_a'(j) = len_a(j - 1) + e_a` for `j >= 1`, and analogously for `len_b`.
//!
//! The profile is in bijection with `L_max_anchor_mutations`, so two histories with the same profile
//! behave identically for every continuation.
//! Hence, the profile can replace the alignment history in search nodes.
//!
//! # Truncation
//!
//! Every mutation covers at least one character, so an anchor-free `L_j` contains at most `2 * anchor_k - 2` mutations.
//! If `max_anchor_mutations >= 2 * anchor_k - 1`, then an anchor-free history therefore has fewer than `max_anchor_mutations` mutations,
//! so `L_j` is the whole history for all `j >= 2 * anchor_k - 2`.
//! Hence, only the entries `0 <= j <= M` with `M = min(max_anchor_mutations, 2 * anchor_k - 1)` are stored,
//! and the extension and validity rules above hold for `M` in place of `max_anchor_mutations`.
//! For supported configurations (see [`ChainingLowerBoundConfig::validate`](crate::config::ChainingLowerBoundConfig::validate)),
//! `max_anchor_mutations < anchor_k`, so `M = max_anchor_mutations`.
//!
//! # Encoding
//!
//! Each of the two vectors `len(0) <= len(1) <= ... <= len(M)` is stored in unary in one `u128`:
//! bit `len(j) + j` is set for each `j`, and all other bits are zero.
//! In other words, reading from the least significant bit, each entry `j` is encoded by `len(j) - len(j - 1)` zeros
//! (with `len(-1) = 0`) followed by a one.
//! With this encoding, the extension rules become shifts:
//!
//! * a match inserts a zero in front, i.e. shifts left by one;
//! * a mutation drops the highest set bit (entry `M`), shifts left by `1 + e` and sets the lowest bit (the new entry `0` is zero).
//!
//! The highest set bit is at position `len(M) + M`, so the history is anchor-free if and only if both words are below `2^(anchor_k + M)`.
//! Since `len(M)` is at most `anchor_k` directly after an extension, the encoding needs `anchor_k + M + 1` bits.
//! Hence, it supports all parameters with `anchor_k + M <= 127`, see [`HistoryProfile::is_supported`].

use std::fmt::{Debug, Display};

use crate::alignment_history::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;

/// The profile of an alignment history, see the [module documentation](self).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HistoryProfile {
    /// The unary encoding of `len_a(0..=M)`.
    len_a: u128,
    /// The unary encoding of `len_b(0..=M)`.
    len_b: u128,
}

impl HistoryProfile {
    /// The number of profile entries per sequence minus one, i.e. `M` in the [module documentation](self).
    fn max_index(anchor_k: u8, max_anchor_mutations: u8) -> u32 {
        u32::from(max_anchor_mutations).min(2 * u32::from(anchor_k).saturating_sub(1) + 1)
    }

    /// The bit index that must not be reached by the encoding of an anchor-free history.
    #[inline]
    fn limit_bit(anchor_k: u8, max_anchor_mutations: u8) -> u32 {
        u32::from(anchor_k) + Self::max_index(anchor_k, max_anchor_mutations)
    }

    /// Returns true if profiles can be stored for the given parameters.
    ///
    /// This is the case if `anchor_k + min(max_anchor_mutations, 2 * anchor_k - 1) <= 127` and `anchor_k >= 1`.
    pub fn is_supported(anchor_k: u8, max_anchor_mutations: u8) -> bool {
        anchor_k >= 1 && Self::limit_bit(anchor_k, max_anchor_mutations) <= 127
    }

    /// The profile of the empty alignment history.
    ///
    /// # Panics
    ///
    /// Panics if the parameters are not supported, see [`Self::is_supported`].
    pub fn new_empty(anchor_k: u8, max_anchor_mutations: u8) -> Self {
        assert!(
            Self::is_supported(anchor_k, max_anchor_mutations),
            "Alignment history profiles do not support anchor_k = {anchor_k} and max_anchor_mutations = {max_anchor_mutations}"
        );
        let ones = (1u128 << (Self::max_index(anchor_k, max_anchor_mutations) + 1)) - 1;
        Self {
            len_a: ones,
            len_b: ones,
        }
    }

    /// If possible, returns the profile of the extension of this history with another alignment history operation.
    ///
    /// If the extension would end in an anchor-like subalignment, then `None` is returned.
    ///
    /// The parameters must be the same as the ones this profile was created with in [`Self::new_empty`].
    #[inline]
    pub fn try_extend(
        self,
        operation: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Option<Self> {
        let (len_a, len_b) = match operation {
            AlignmentHistoryOperation::Match => (self.len_a << 1, self.len_b << 1),
            AlignmentHistoryOperation::Substitution => (
                extend_mutation(self.len_a, 1),
                extend_mutation(self.len_b, 1),
            ),
            AlignmentHistoryOperation::GapInA => (
                extend_mutation(self.len_a, 0),
                extend_mutation(self.len_b, 1),
            ),
            AlignmentHistoryOperation::GapInB => (
                extend_mutation(self.len_a, 1),
                extend_mutation(self.len_b, 0),
            ),
        };

        let limit_bit = Self::limit_bit(anchor_k, max_anchor_mutations);
        debug_assert!(limit_bit <= 127);
        ((len_a | len_b) >> limit_bit == 0).then_some(Self { len_a, len_b })
    }

    /// Returns the entries `len_a(0..=M)` of the profile.
    pub fn len_a(&self) -> impl Iterator<Item = usize> {
        decode(self.len_a)
    }

    /// Returns the entries `len_b(0..=M)` of the profile.
    pub fn len_b(&self) -> impl Iterator<Item = usize> {
        decode(self.len_b)
    }
}

/// Extend a unary encoded profile vector by a mutation that covers `extension` characters of the sequence.
#[inline]
fn extend_mutation(encoding: u128, extension: u32) -> u128 {
    debug_assert_ne!(encoding, 0);
    // Drop entry `M`, which is the highest set bit.
    let highest_bit = 1u128 << (127 - encoding.leading_zeros());
    ((encoding ^ highest_bit) << (1 + extension)) | 1
}

/// Decode a unary encoded profile vector.
fn decode(mut encoding: u128) -> impl Iterator<Item = usize> {
    let mut index = 0;
    std::iter::from_fn(move || {
        if encoding == 0 {
            None
        } else {
            let position = encoding.trailing_zeros();
            encoding &= encoding - 1;
            let result = usize::try_from(position).unwrap() - index;
            index += 1;
            Some(result)
        }
    })
}

impl Display for HistoryProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "a[")?;
        for (index, len) in self.len_a().enumerate() {
            if index > 0 {
                write!(f, ",")?;
            }
            write!(f, "{len}")?;
        }
        write!(f, "]b[")?;
        for (index, len) in self.len_b().enumerate() {
            if index > 0 {
                write!(f, ",")?;
            }
            write!(f, "{len}")?;
        }
        write!(f, "]")
    }
}

impl Debug for HistoryProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HistoryProfile({self})")
    }
}
