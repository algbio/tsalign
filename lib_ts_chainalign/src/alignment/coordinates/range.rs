use std::ops::Range;

use crate::alignment::{
    coordinates::{
        AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
        SpecificSecondaryAlignmentCoordinates,
    },
    ts_kind::TsKind,
};

/// A range of alignment coordinates in the primary sequence space.
/// The offset is inclusive, the limit is exclusive.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct PrimaryAlignmentRange {
    /// Inclusive start of the range.
    offset: PrimaryAlignmentCoordinates,
    /// Exclusive end of the range.
    limit: PrimaryAlignmentCoordinates,
}

/// A range of alignment coordinates in the secondary sequence space, without specifying which secondary sequence space.
/// The offset is inclusive, the limit is exclusive.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct AnySecondaryAlignmentRange {
    /// Inclusive start of the range.
    offset: AnySecondaryAlignmentCoordinates,
    /// Exclusive end of the range.
    limit: AnySecondaryAlignmentCoordinates,
}

/// A range of alignment coordinates in a specified secondary sequence space.
/// The offset is inclusive, the limit is exclusive.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct SpecificSecondaryAlignmentRange {
    range: AnySecondaryAlignmentRange,
    ts_kind: TsKind,
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

    pub fn new_equal_length(offset: PrimaryAlignmentCoordinates, length: usize) -> Self {
        let limit = PrimaryAlignmentCoordinates::new(offset.a() + length, offset.b() + length);
        Self::new(offset, limit)
    }

    pub fn offset(&self) -> PrimaryAlignmentCoordinates {
        self.offset
    }

    pub fn limit(&self) -> PrimaryAlignmentCoordinates {
        self.limit
    }

    pub fn len_a(&self) -> usize {
        self.limit.a() - self.offset.a()
    }

    pub fn len_b(&self) -> usize {
        self.limit.b() - self.offset.b()
    }
}

impl AnySecondaryAlignmentRange {
    pub fn new(
        offset: AnySecondaryAlignmentCoordinates,
        limit: AnySecondaryAlignmentCoordinates,
    ) -> Self {
        assert!(offset <= limit);
        Self { offset, limit }
    }

    pub fn new_from_ranges(
        ancestor_offset: usize,
        ancestor_limit: usize,
        descendant: Range<usize>,
    ) -> Self {
        debug_assert!(ancestor_offset >= ancestor_limit);
        let offset = AnySecondaryAlignmentCoordinates::new(ancestor_offset, descendant.start);
        let limit = AnySecondaryAlignmentCoordinates::new(ancestor_limit, descendant.end);
        Self::new(offset, limit)
    }

    pub fn new_equal_length(offset: AnySecondaryAlignmentCoordinates, length: usize) -> Self {
        let limit = AnySecondaryAlignmentCoordinates::new(
            offset.ancestor().checked_sub(length).unwrap(),
            offset.descendant() + length,
        );
        Self::new(offset, limit)
    }

    pub fn offset(&self) -> AnySecondaryAlignmentCoordinates {
        self.offset
    }

    pub fn limit(&self) -> AnySecondaryAlignmentCoordinates {
        self.limit
    }

    pub fn into_specific(self, ts_kind: TsKind) -> SpecificSecondaryAlignmentRange {
        SpecificSecondaryAlignmentRange::new(self.offset, self.limit, ts_kind)
    }

    pub fn len_ancestor(&self) -> usize {
        self.offset
            .ancestor()
            .checked_sub(self.limit.ancestor())
            .unwrap()
    }

    pub fn len_descendant(&self) -> usize {
        self.limit
            .descendant()
            .checked_sub(self.offset.descendant())
            .unwrap()
    }
}

impl SpecificSecondaryAlignmentRange {
    pub fn new(
        offset: AnySecondaryAlignmentCoordinates,
        limit: AnySecondaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> Self {
        assert!(offset <= limit);
        Self {
            range: AnySecondaryAlignmentRange::new(offset, limit),
            ts_kind,
        }
    }

    pub fn new_from_ranges(
        ancestor: Range<usize>,
        descendant: Range<usize>,
        ts_kind: TsKind,
    ) -> Self {
        let offset = AnySecondaryAlignmentCoordinates::new(ancestor.start, descendant.start);
        let limit = AnySecondaryAlignmentCoordinates::new(ancestor.end, descendant.end);
        Self::new(offset, limit, ts_kind)
    }

    pub fn new_equal_length(offset: SpecificSecondaryAlignmentCoordinates, length: usize) -> Self {
        let limit = AnySecondaryAlignmentCoordinates::new(
            offset.ancestor().checked_sub(length).unwrap(),
            offset.descendant() + length,
        );
        Self::new(offset.into(), limit, offset.ts_kind())
    }

    pub fn offset(&self) -> AnySecondaryAlignmentCoordinates {
        self.range.offset
    }

    pub fn limit(&self) -> AnySecondaryAlignmentCoordinates {
        self.range.limit
    }

    pub fn ts_kind(&self) -> TsKind {
        self.ts_kind
    }

    pub fn into_any(self) -> AnySecondaryAlignmentRange {
        self.range
    }
}

impl std::fmt::Display for PrimaryAlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "R[{}, {})", self.offset, self.limit)
    }
}

impl std::fmt::Display for AnySecondaryAlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "R[{}, {})", self.offset, self.limit)
    }
}

impl std::fmt::Display for SpecificSecondaryAlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "R{}[{}, {})",
            self.ts_kind.digits(),
            self.range.offset,
            self.range.limit,
        )
    }
}

impl From<SpecificSecondaryAlignmentRange> for AnySecondaryAlignmentRange {
    fn from(value: SpecificSecondaryAlignmentRange) -> Self {
        value.range
    }
}
