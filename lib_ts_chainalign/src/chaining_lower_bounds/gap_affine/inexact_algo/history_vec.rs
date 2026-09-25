use std::{fmt::Display, hash::Hash, iter};

use num_traits::{NumCast, PrimInt, Unsigned};

use crate::chaining_lower_bounds::gap_affine::inexact_algo::history_alignment_operations::AlignmentHistoryOperation;

#[cfg(test)]
mod tests;

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

pub trait HistoryInt: PrimInt + Unsigned + NumCast + Hash {
    /// The number of bits used for storing the length of the history.
    ///
    /// As the history is always an even number of bits, this is rounded up correspondingly.
    const LEN_BITS: usize;
    /// The number of bits used for storing the history itself.
    ///
    /// This is always even.
    const HISTORY_BITS: usize;
    /// A mask for the lowest [`Self::len_bits()`] bits of the integer type.
    const LEN_MASK: Self;
    /// A mask for the highest [`Self::history_bits()`] bits of the integer type.
    const HISTORY_MASK: Self;

    /// Returns the maximum length of the alignment history vector that can be stored in this integer type.
    fn max_len() -> usize {
        Self::HISTORY_BITS / 2
    }
}

/// A packed representation of an alignment history.
///
/// The representation is just one unsigned integer.
/// The lowest [`Self::len_bits()`] bits are used to store the length of the history in number of operations (each operation takes two bits).
/// The remaining [`Self::history_bits()`] bits are used to store the history itself.
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

    fn try_extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Option<Self> {
        let anchor_k = <usize as From<u8>>::from(anchor_k);
        let max_anchor_mutations = <usize as From<u8>>::from(max_anchor_mutations);

        // We probably can simplify how many numbers we track here, but we keep it like this because it is easier to reason about the correctness of the algorithm.
        // Hopefully the optimiser will take good care of this.
        let mut len_a = 0;
        let mut len_b = 0;
        let mut mutations = 0;
        let mut previous_len_a = 0;
        let mut previous_len_b = 0;
        let mut previous_mutations = 0;
        let extension: Self = iter::once(alignment)
            .chain(*self)
            .filter(|alignment| {
                // Extend until at least one of `is_anchor_length` or `is_invalid` becomes true.
                let result =
                    if len_a < anchor_k && len_b < anchor_k && mutations <= max_anchor_mutations {
                        len_a += alignment.len_a();
                        len_b += alignment.len_b();
                        mutations += alignment.mutations();
                        true
                    } else {
                        false
                    };

                previous_len_a += alignment.len_a();
                previous_len_b += alignment.len_b();
                previous_mutations += alignment.mutations();

                result
            })
            .collect();

        let len_a = len_a;
        let len_b = len_b;
        let mutations = mutations;
        let previous_len_a = previous_len_a - alignment.len_a();
        let previous_len_b = previous_len_b - alignment.len_b();
        let previous_mutations = previous_mutations - alignment.mutations();

        let previous_is_anchor_length = previous_len_a.max(previous_len_b) >= anchor_k;
        let previous_is_invalid = previous_mutations > max_anchor_mutations;
        let previous_general_phase = previous_is_anchor_length || previous_is_invalid;

        let current_is_anchor_length = len_a.max(len_b) >= anchor_k;
        let current_is_invalid = mutations > max_anchor_mutations;
        let current_is_general_phase = current_is_anchor_length || current_is_invalid;

        match (previous_general_phase, current_is_general_phase) {
            (_, true) => current_is_invalid.then_some(extension),
            (true, false) => unreachable!(
                "It should not be possible to transition from the general phase back to the initial phase."
            ),
            (false, false) => {
                // In the initial phase it is possible that some histories cannot be extended in to an invalid anchor, because the contain too many matches.
                // However, with inexact anchors, we chain without overlap, so even sequences of `anchor_k - 1` matches are required to cover all possible alignments.
                Some(extension)
            }
        }
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
