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
//! Both rules are monotone in the entries of the profile:
//! if every entry of a profile `P` is at most the corresponding entry of a profile `Q`,
//! then this stays true after extending both by the same operation.
//! Hence, every continuation that is anchor-free after `Q` is also anchor-free after `P`
//! (see [`HistoryProfile::is_componentwise_le`]).
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
//! The profile is stored in a single `u64` with one byte per entry:
//! byte `j` holds `len_a(j)` and byte `4 + j` holds `len_b(j)` for `0 <= j <= M`, and all other bytes are zero.
//! Hence, `M <= 3` is supported.
//! Directly after an extension, every entry is at most `anchor_k` (before, every entry is at most `anchor_k - 1`),
//! so with `anchor_k <= 127` every byte stays below `128`, see [`HistoryProfile::is_supported`].
//! The extension rules work on both sequences at once:
//!
//! * a match adds one to every used byte;
//! * a mutation shifts the word left by one byte and clears the bytes that are unused or hold an entry `0`
//!   (this drops the old entries `M`), and then adds `e_a` to the bytes `1..=M` and `e_b` to the bytes `5..=4 + M`.
//!
//! Validity is checked by adding `128 - anchor_k` to the bytes of the two entries `M` and testing their highest bits.
//! The constants for these operations are precomputed in [`HistoryProfileParameters`].
//! Since all bytes are below `128`, the componentwise comparison of two profiles is a single SWAR subtraction,
//! see [`HistoryProfile::is_componentwise_le`].

use std::fmt::{Debug, Display};

use crate::alignment_history::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;

/// The highest bit of each byte.
const HIGH_BITS: u64 = 0x8080_8080_8080_8080;

/// The byte offset of the entries of sequence b.
const B_OFFSET: u32 = 4;

/// The largest supported `M`, i.e. the largest supported number of profile entries per sequence minus one.
pub const MAX_PROFILE_INDEX: u8 = 3;

/// The largest supported `anchor_k`.
pub const MAX_PROFILE_ANCHOR_K: u8 = 127;

/// The profile of an alignment history, see the [module documentation](self).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HistoryProfile {
    /// The byte encoding of `len_a(0..=M)` and `len_b(0..=M)`.
    lens: u64,
}

/// Precomputed constants for the [`HistoryProfile`]s of fixed `anchor_k` and `max_anchor_mutations`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HistoryProfileParameters {
    /// `M` in the [module documentation](self).
    max_index: u8,
    /// One in every used byte.
    match_increment: u64,
    /// `0xff` in every used byte except the bytes of the entries `0`.
    mutation_mask: u64,
    /// One in the bytes `1..=M`.
    mutation_increment_a: u64,
    /// One in the bytes `5..=4 + M`.
    mutation_increment_b: u64,
    /// `128 - anchor_k` in the bytes of the entries `M`.
    limit_increment: u64,
    /// The highest bit of the bytes of the entries `M`.
    limit_bits: u64,
}

impl HistoryProfileParameters {
    /// Precompute the constants for the given parameters.
    ///
    /// # Panics
    ///
    /// Panics if the parameters are not supported, see [`HistoryProfile::is_supported`].
    pub fn new(anchor_k: u8, max_anchor_mutations: u8) -> Self {
        assert!(
            HistoryProfile::is_supported(anchor_k, max_anchor_mutations),
            "Alignment history profiles do not support anchor_k = {anchor_k} and max_anchor_mutations = {max_anchor_mutations}; \
             they need 1 <= anchor_k <= {MAX_PROFILE_ANCHOR_K} and max_anchor_mutations <= {MAX_PROFILE_INDEX}"
        );
        let max_index = HistoryProfile::max_index(anchor_k, max_anchor_mutations);
        let m = u32::from(max_index);

        // One in the bytes `0..=M`.
        let ones_a = (1u64 << (8 * (m + 1))) / 0xff;
        let ones = ones_a | (ones_a << (8 * B_OFFSET));
        let entries_zero = 1 | (1 << (8 * B_OFFSET));
        let entries_max = (1 << (8 * m)) | (1 << (8 * (B_OFFSET + m)));
        let mutation_increment_a = ones_a & !1;

        Self {
            max_index,
            match_increment: ones,
            mutation_mask: (ones & !entries_zero) * 0xff,
            mutation_increment_a,
            mutation_increment_b: mutation_increment_a << (8 * B_OFFSET),
            limit_increment: entries_max * u64::from(128 - anchor_k),
            limit_bits: entries_max * 0x80,
        }
    }
}

impl HistoryProfile {
    /// The number of profile entries per sequence minus one, i.e. `M` in the [module documentation](self).
    fn max_index(anchor_k: u8, max_anchor_mutations: u8) -> u8 {
        let cap = (2 * u16::from(anchor_k)).saturating_sub(1);
        u8::try_from(u16::from(max_anchor_mutations).min(cap)).unwrap()
    }

    /// Returns true if profiles can be stored for the given parameters.
    ///
    /// This is the case if `1 <= anchor_k <= 127` and `min(max_anchor_mutations, 2 * anchor_k - 1) <= 3`.
    pub fn is_supported(anchor_k: u8, max_anchor_mutations: u8) -> bool {
        (1..=MAX_PROFILE_ANCHOR_K).contains(&anchor_k)
            && Self::max_index(anchor_k, max_anchor_mutations) <= MAX_PROFILE_INDEX
    }

    /// The profile of the empty alignment history.
    pub fn empty() -> Self {
        Self { lens: 0 }
    }

    /// If possible, returns the profile of the extension of this history with another alignment history operation.
    ///
    /// If the extension would end in an anchor-like subalignment, then `None` is returned.
    ///
    /// All profiles that are compared or combined must use the same `parameters`.
    #[inline]
    pub fn try_extend(
        self,
        operation: AlignmentHistoryOperation,
        parameters: HistoryProfileParameters,
    ) -> Option<Self> {
        let shifted = (self.lens << 8) & parameters.mutation_mask;
        let lens = match operation {
            AlignmentHistoryOperation::Match => self.lens + parameters.match_increment,
            AlignmentHistoryOperation::Substitution => {
                shifted + parameters.mutation_increment_a + parameters.mutation_increment_b
            }
            AlignmentHistoryOperation::GapInA => shifted + parameters.mutation_increment_b,
            AlignmentHistoryOperation::GapInB => shifted + parameters.mutation_increment_a,
        };
        debug_assert_eq!(lens & HIGH_BITS, 0, "{self} {operation}");

        ((lens + parameters.limit_increment) & parameters.limit_bits == 0).then_some(Self { lens })
    }

    /// Returns true if every entry of this profile is at most the corresponding entry of `other`.
    ///
    /// Then every continuation that is anchor-free after `other` is also anchor-free after this profile,
    /// see the [module documentation](self).
    /// Both profiles must use the same parameters.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "domination test for the planned Pareto prune list (issue #224)"
        )
    )]
    pub fn is_componentwise_le(&self, other: &Self) -> bool {
        // All bytes are below 128, so each byte of `other | HIGH_BITS` is larger than the corresponding byte of `self`,
        // and the subtraction does not borrow across bytes.
        // The highest bit of a byte of the difference is set if and only if the byte of `other` is at least the byte of `self`.
        ((other.lens | HIGH_BITS) - self.lens) & HIGH_BITS == HIGH_BITS
    }

    /// Returns the entries `len_a(0..=M)` of the profile.
    #[cfg(test)]
    pub fn len_a(&self, parameters: HistoryProfileParameters) -> impl Iterator<Item = usize> {
        self.entries(0, parameters.max_index)
    }

    /// Returns the entries `len_b(0..=M)` of the profile.
    #[cfg(test)]
    pub fn len_b(&self, parameters: HistoryProfileParameters) -> impl Iterator<Item = usize> {
        self.entries(B_OFFSET, parameters.max_index)
    }

    /// Returns the bytes `byte_offset..=byte_offset + max_index`.
    fn entries(&self, byte_offset: u32, max_index: u8) -> impl Iterator<Item = usize> {
        let lens = self.lens;
        (0..=u32::from(max_index))
            .map(move |j| usize::from((lens >> (8 * (byte_offset + j))) as u8))
    }
}

impl Display for HistoryProfile {
    /// Displays the entries of the profile.
    ///
    /// Since the profile does not know `M`, all four entries per sequence are shown, including unused ones.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (name, byte_offset) in [("a", 0), ("b", B_OFFSET)] {
            write!(f, "{name}[")?;
            for (index, len) in self.entries(byte_offset, MAX_PROFILE_INDEX).enumerate() {
                if index > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{len}")?;
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

impl Debug for HistoryProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HistoryProfile({self})")
    }
}
