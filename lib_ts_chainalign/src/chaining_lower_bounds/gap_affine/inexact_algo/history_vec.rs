use std::hash::Hash;

use num_traits::{NumCast, PrimInt, Unsigned};

use crate::chaining_lower_bounds::gap_affine::inexact_algo::history_alignment_operations::AlignmentHistoryOperation;

pub trait AlignmentHistory: Default + Eq + Hash + Ord + Copy {
    /// Returns the length of the alignment history vector.
    fn len(&self) -> usize;

    /// Returns true if the alignment history vector is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns true if the alignment history vector can be extended with another alignment history operation without becoming a valid anchor.
    fn can_extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> bool;

    /// Return the extension of this the alignment history vector with another alignment history operation.
    ///
    /// This method panics in debug mode if the alignment history cannot be extended with the given alignment history operation.
    /// In release mode, this error is silently ignored.
    fn extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self;
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
        (Self::HISTORY_BITS / Self::from(2u8).unwrap())
            .to_usize()
            .unwrap()
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

    /// Sets the length of the alignment history vector.
    ///
    /// Note that this does not alter the history, so the vector may become invalid.
    fn set_len(&mut self, len: usize) {
        debug_assert!(len <= UnsignedInt::max_len());

        let len = UnsignedInt::from(len).unwrap();
        self.data = (self.data & UnsignedInt::HISTORY_MASK) | (len & UnsignedInt::LEN_MASK);
    }

    /// Adds an element at the end of the alignment history vector.
    ///
    /// If the vector is full, then the first element will be deleted and all previous elements will be shifted to the left by one index.
    fn push_back(self, alignment: AlignmentHistoryOperation) -> Self {
        let len = self.len();
        let alignment_bits = UnsignedInt::from(alignment.to_bits()).unwrap();

        if len == UnsignedInt::max_len() {
            // Shift the history to the left by one index and add the new alignment operation at the end.
            let len = UnsignedInt::from(len).unwrap();
            let data = (self.data >> 2) | (alignment_bits << (UnsignedInt::HISTORY_BITS - 2));
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
        if self.is_empty() { None } else { todo!() }
    }
}

impl<UnsignedInt: HistoryInt> AlignmentHistory for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    fn len(&self) -> usize {
        (self.data & UnsignedInt::LEN_MASK).try_into().ok().unwrap()
    }

    fn can_extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> bool {
        todo!()
    }

    fn extend(
        &self,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self {
        debug_assert!(self.can_extend(alignment, anchor_k, max_anchor_mutations));

        todo!()
    }
}

impl<UnsignedInt: HistoryInt> Default for UnsignedIntAlignmentHistoryVec<UnsignedInt> {
    fn default() -> Self {
        Self::new()
    }
}
