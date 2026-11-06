use std::{fmt::Display, hash::Hash, iter};

use num_traits::{NumCast, PrimInt, Unsigned};

use crate::alignment_history::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;

/// Keeps track of the most recent alignments required to reach a DP node.
///
/// Only keeps track of enough state to prove that a DP node does not end on a possibly valid anchor,
/// i.e. that when backtracking the node, at least `max_anchor_mutations + 1` mutations are required to create an anchor of length `anchor_k`.
///
/// # Properties
///
/// The history stores the longest suffix of the alignment that contains at most `max_anchor_mutations` mutations, from newest to oldest operation.
/// The alignment ending in this history is valid (i.e. contains no subalignment that could be an anchor)
/// if and only if this suffix covers fewer than `anchor_k` characters of each of the two sequences.
///
/// Operations before the suffix are never required.
/// Extending by a match extends the suffix, and extending by a mutation extends the suffix and drops everything up to and including its oldest mutation, if it already had `max_anchor_mutations` mutations.
/// In both cases, the new suffix is a suffix of the old suffix extended by the new operation.
///
/// See [`HistoryProfile`](crate::alignment_history::profile::HistoryProfile) for a formalisation.
pub trait AlignmentHistory: Default + Eq + Hash + Ord + Copy {
    /// Returns true if the alignment history has enough capacity for the given parameters.
    fn has_enough_capacity(anchor_k: u8, max_anchor_mutations: u8) -> bool;

    /// Returns the length of the alignment history vector.
    fn len(&self) -> usize;

    /// Returns true if the alignment history vector is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the operations of the alignment history, from newest to oldest.
    fn operations(&self) -> impl Iterator<Item = AlignmentHistoryOperation>;

    /// Returns the maximum length that the history can have given the parameters.
    fn required_capacity(anchor_k: u8, max_anchor_mutations: u8) -> usize {
        let anchor_k = <usize as From<u8>>::from(anchor_k);
        let max_anchor_mutations = <usize as From<u8>>::from(max_anchor_mutations);

        // The history contains at most `max_anchor_mutations` mutations and covers at most `anchor_k - 1` characters of each sequence.
        // Matches and substitutions cover a character of both sequences, while gaps cover a character of only one sequence.
        // Hence, the longest history consists of `d` gaps in b, `i` gaps in a and `anchor_k - 1 - max(d, i)` matches,
        // with `d + i <= max_anchor_mutations` and `max(d, i) <= anchor_k - 1`.
        // Its length is `anchor_k - 1 + min(d, i)`, which is maximised by `min(d, i) = min(floor(max_anchor_mutations / 2), anchor_k - 1)`.
        let max_len_per_sequence = anchor_k.saturating_sub(1);
        max_len_per_sequence + (max_anchor_mutations / 2).min(max_len_per_sequence)
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

pub trait HistoryInt: PrimInt + Unsigned + NumCast + Hash {
    /// The number of bits used for storing the length of the history.
    ///
    /// As the history is always an even number of bits, this is rounded up correspondingly.
    const LEN_BITS: usize;
    /// The number of bits used for storing the history itself.
    ///
    /// This is always even.
    const HISTORY_BITS: usize;
    /// A mask for the lowest [`Self::LEN_BITS`] bits of the integer type.
    const LEN_MASK: Self;
    /// A mask for the highest [`Self::HISTORY_BITS`] bits of the integer type.
    const HISTORY_MASK: Self;

    /// Returns the maximum length of the alignment history vector that can be stored in this integer type.
    fn max_len() -> usize {
        Self::HISTORY_BITS / 2
    }
}

/// A packed representation of an alignment history.
///
/// The representation is just one unsigned integer.
/// The lowest [`HistoryInt::LEN_BITS`] bits are used to store the length of the history in number of operations (each operation takes two bits).
/// The remaining [`HistoryInt::HISTORY_BITS`] bits are used to store the history itself.
/// The first history item is stored in the lowest bits of the history part.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    data: UnsignedInt,
}

impl<UnsignedInt: HistoryInt> UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    /// Create a new empty alignment history vector.
    pub fn new() -> Self {
        Self {
            data: UnsignedInt::zero(),
        }
    }

    /// Adds an element at the end of the alignment history vector.
    ///
    /// If the vector is full, then the first element will be deleted and all previous elements will be shifted to the left by one index.
    #[allow(dead_code)]
    fn push_back(self, alignment: AlignmentHistoryOperation) -> Self {
        let len = self.len();
        let alignment_bits = UnsignedInt::from(alignment.to_bits()).unwrap();

        if len == UnsignedInt::max_len() {
            // Shift the history to the left by one index and add the new alignment operation at the end.
            let len = UnsignedInt::from(len).unwrap();
            let data = (self.data >> 2)
                | (alignment_bits << (UnsignedInt::HISTORY_BITS + UnsignedInt::LEN_BITS - 2));
            Self {
                data: (data & UnsignedInt::HISTORY_MASK) | (len & UnsignedInt::LEN_MASK),
            }
        } else {
            // Add the new alignment operation at the end.
            let alignment = UnsignedInt::from(alignment.to_bits()).unwrap();
            let data = (self.data) | (alignment) << (UnsignedInt::LEN_BITS + 2 * len);
            let len = UnsignedInt::from(len + 1).unwrap();

            Self {
                data: (data & UnsignedInt::HISTORY_MASK) | (len & UnsignedInt::LEN_MASK),
            }
        }
    }

    /// Removes the element at the front of the alignment history vector and returns it.
    fn pop_front(self) -> Option<(Self, AlignmentHistoryOperation)> {
        if self.is_empty() {
            None
        } else {
            let len = self.len();
            let alignment_bits =
                (self.data >> UnsignedInt::LEN_BITS) & UnsignedInt::from(0b11).unwrap();
            let alignment =
                AlignmentHistoryOperation::from_bits(alignment_bits.to_usize().unwrap()).unwrap();

            let len = len - 1;
            let data =
                (self.data >> 2) & UnsignedInt::HISTORY_MASK | UnsignedInt::from(len).unwrap();

            Some((Self { data }, alignment))
        }
    }
}

impl<UnsignedInt: HistoryInt> AlignmentHistory for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    fn has_enough_capacity(anchor_k: u8, max_anchor_mutations: u8) -> bool {
        Self::required_capacity(anchor_k, max_anchor_mutations) <= UnsignedInt::max_len()
    }

    fn len(&self) -> usize {
        (self.data & UnsignedInt::LEN_MASK).to_usize().unwrap()
    }

    fn operations(&self) -> impl Iterator<Item = AlignmentHistoryOperation> {
        *self
    }

    fn try_extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Option<Self> {
        let anchor_k = <usize as From<u8>>::from(anchor_k);
        let max_anchor_mutations = <usize as From<u8>>::from(max_anchor_mutations);

        // Find the longest suffix of the extended alignment with at most `max_anchor_mutations` mutations.
        // It is a suffix of the extended history, since the history is the longest such suffix before the extension.
        let mut len_a = 0;
        let mut len_b = 0;
        let mut mutations = 0;
        let mut len = 0;
        for alignment in iter::once(alignment).chain(*self) {
            mutations += alignment.mutations();
            if mutations > max_anchor_mutations {
                break;
            }

            len_a += alignment.len_a();
            len_b += alignment.len_b();
            len += 1;

            if len_a >= anchor_k || len_b >= anchor_k {
                // The suffix could be a valid anchor.
                return None;
            }
        }

        Some(iter::once(alignment).chain(*self).take(len).collect())
    }
}

impl<UnsignedInt: HistoryInt> Default for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    fn default() -> Self {
        Self::new()
    }
}

impl<UnsignedInt: HistoryInt> Iterator for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    type Item = AlignmentHistoryOperation;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_empty() {
            None
        } else {
            let (vec, alignment) = self.pop_front().unwrap();
            *self = vec;
            Some(alignment)
        }
    }
}

impl<UnsignedInt: HistoryInt> FromIterator<AlignmentHistoryOperation>
    for UnsignedIntAlignmentHistoryVec<UnsignedInt>
{
    fn from_iter<T: IntoIterator<Item = AlignmentHistoryOperation>>(iter: T) -> Self {
        let mut iter = iter.into_iter();
        let mut len = 0;
        let mut data = UnsignedInt::zero();
        for alignment in iter.by_ref().take(UnsignedInt::max_len()) {
            let alignment_bits = UnsignedInt::from(alignment.to_bits()).unwrap();
            data = data | (alignment_bits << (UnsignedInt::LEN_BITS + 2 * len));
            len += 1;
        }

        assert!(
            iter.next().is_none(),
            "Alignment history vector is too long",
        );

        let data = data | UnsignedInt::from(len).unwrap();
        Self { data }
    }
}

impl<UnsignedInt: HistoryInt> Display for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for alignment in *self {
            write!(f, "{alignment}")?;
        }
        Ok(())
    }
}

impl HistoryInt for u8 {
    const LEN_BITS: usize = 2;
    const HISTORY_BITS: usize = 6;
    const LEN_MASK: Self = (1 << Self::LEN_BITS) - 1;
    const HISTORY_MASK: Self = ((1 << Self::HISTORY_BITS) - 1) << Self::LEN_BITS;
}

impl HistoryInt for u16 {
    const LEN_BITS: usize = 4;
    const HISTORY_BITS: usize = 12;
    const LEN_MASK: Self = (1 << Self::LEN_BITS) - 1;
    const HISTORY_MASK: Self = ((1 << Self::HISTORY_BITS) - 1) << Self::LEN_BITS;
}

impl HistoryInt for u32 {
    const LEN_BITS: usize = 4;
    const HISTORY_BITS: usize = 28;
    const LEN_MASK: Self = (1 << Self::LEN_BITS) - 1;
    const HISTORY_MASK: Self = ((1 << Self::HISTORY_BITS) - 1) << Self::LEN_BITS;
}

impl HistoryInt for u64 {
    const LEN_BITS: usize = 6;
    const HISTORY_BITS: usize = 58;
    const LEN_MASK: Self = (1 << Self::LEN_BITS) - 1;
    const HISTORY_MASK: Self = ((1 << Self::HISTORY_BITS) - 1) << Self::LEN_BITS;
}

impl HistoryInt for u128 {
    const LEN_BITS: usize = 6;
    const HISTORY_BITS: usize = 122;
    const LEN_MASK: Self = (1 << Self::LEN_BITS) - 1;
    const HISTORY_MASK: Self = ((1 << Self::HISTORY_BITS) - 1) << Self::LEN_BITS;
}
