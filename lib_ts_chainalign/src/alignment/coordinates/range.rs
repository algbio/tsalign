use std::ops::Range;

use crate::alignment::coordinates::PrimaryAlignmentCoordinates;

/// A range of alignment coordinates in the primary sequence space.
/// The offset is inclusive, the limit is exclusive.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct PrimaryAlignmentRange {
    /// Inclusive start of the range.
    offset: PrimaryAlignmentCoordinates,
    /// Exclusive end of the range.
    limit: PrimaryAlignmentCoordinates,
}

impl PrimaryAlignmentRange {
    pub fn new(offset: PrimaryAlignmentCoordinates, limit: PrimaryAlignmentCoordinates) -> Self {
        assert!(offset <= limit);
        Self { offset, limit }
    }

    pub fn new_from_ranges(seq1: Range<usize>, seq2: Range<usize>) -> Self {
        let offset = PrimaryAlignmentCoordinates::new(seq1.start, seq2.start);
        let limit = PrimaryAlignmentCoordinates::new(seq1.end, seq2.end);
        Self::new(offset, limit)
    }

    pub fn offset(&self) -> PrimaryAlignmentCoordinates {
        self.offset
    }

    pub fn limit(&self) -> PrimaryAlignmentCoordinates {
        self.limit
    }
}

impl std::fmt::Display for PrimaryAlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {})", self.offset, self.limit)
    }
}