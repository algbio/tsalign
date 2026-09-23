use std::{fmt::Display, ops::Range};

use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::{
            AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            range::AnySecondaryAlignmentRange,
        },
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::primary::PrimaryAnchor,
};

/// A secondary anchor.
///
/// This is an anchor between the ancestor in reverse direction and the descendant in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its ancestor ordinate and finally by its descendant ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SecondaryAnchor<Cost> {
    range: AnySecondaryAlignmentRange,
    pub(super) cost: Cost,
}

impl<Cost> SecondaryAnchor<Cost> {
    pub fn new(range: AnySecondaryAlignmentRange, cost: Cost) -> Self {
        Self { range, cost }
    }

    pub fn new_from_ranges(
        ancestor_offset: usize,
        ancestor_limit: usize,
        descendant: Range<usize>,
        cost: Cost,
    ) -> Self {
        Self::new(
            AnySecondaryAlignmentRange::new_from_ranges(
                ancestor_offset,
                ancestor_limit,
                descendant,
            ),
            cost,
        )
    }

    pub fn new_exact(offset: AnySecondaryAlignmentCoordinates, length: usize) -> Self
    where
        Cost: Zero,
    {
        Self::new(
            AnySecondaryAlignmentRange::new_equal_length(offset, length),
            Cost::zero(),
        )
    }

    pub fn new_from_inexact_alignment_anchor(
        inexact_alignment_anchors::anchor::Anchor {
            offset_a,
            limit_a,
            offset_b,
            limit_b,
            cost,
        }: inexact_alignment_anchors::anchor::Anchor<Cost>,
        ancestor_len: usize,
    ) -> Self {
        Self::new(
            AnySecondaryAlignmentRange::new_from_ranges(
                ancestor_len - offset_a,
                ancestor_len - limit_a,
                offset_b..limit_b,
            ),
            cost,
        )
    }

    pub fn range(&self) -> AnySecondaryAlignmentRange {
        self.range
    }

    pub fn start(&self) -> AnySecondaryAlignmentCoordinates {
        self.range.offset()
    }

    pub fn end(&self) -> AnySecondaryAlignmentCoordinates {
        self.range.limit()
    }

    pub fn cost(&self) -> Cost
    where
        Cost: Copy,
    {
        self.cost
    }

    /// Returns true if the anchor is at the given range.
    ///
    /// Does not check the `ts_kind`, and will produce false positives if the range given has the wrong `ts_kind`.
    pub fn is_at(&self, range: AnySecondaryAlignmentRange) -> bool {
        self.range == range
    }

    pub fn chaining_gaps(&self, second: &Self) -> Option<(usize, usize)> {
        let gap_start = self.end();
        let gap_end = second.start();

        let gap1 = gap_start.ancestor().checked_sub(gap_end.ancestor())?;
        let gap2 = gap_end.descendant().checked_sub(gap_start.descendant())?;

        Some((gap1, gap2))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &PrimaryAnchor<Cost>,
        ts_kind: TsKind,
    ) -> Option<usize> {
        let gap_start = self.end();
        let gap_end = second.start();

        let gap_start = gap_start.descendant();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_end.a(),
            TsDescendant::Seq2 => gap_end.b(),
        };

        gap_end.checked_sub(gap_start)
    }

    pub fn chaining_jump_gap_from_start(
        &self,
        start: PrimaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> usize
    where
        Cost: Display,
    {
        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => start.a(),
            TsDescendant::Seq2 => start.b(),
        };
        let gap_end = self.start().descendant();

        gap_end
            .checked_sub(gap_start)
            .unwrap_or_else(|| panic!("self: {self}, start: {start}, ts_kind: {ts_kind}"))
    }

    pub fn chaining_jump_gap_to_end(
        &self,
        end: PrimaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> usize {
        let gap_start = self.end().descendant();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => end.a(),
            TsDescendant::Seq2 => end.b(),
        };

        gap_end.checked_sub(gap_start).unwrap()
    }

    /// Returns true if this anchor preceedes the other anchor with an overlap of k-1 characters in both sequences.
    /// Also returns true only if both anchors have zero cost.
    #[deprecated(
        note = "We need to redefine this based on how we want to chain things in the future."
    )]
    pub fn is_direct_free_predecessor_of(&self, successor: &Self) -> bool
    where
        Cost: Zero,
    {
        // TODO how do we account for predecessors that are inexact anchors? do we need to?
        // Maybe trim the exact matches at the start and end of inexact anchors?
        self.range.offset().increment_both(1) == successor.range.offset()
            && self.range.limit().increment_both(1) == successor.range.limit()
            && self.cost.is_zero()
            && successor.cost.is_zero()
    }

    /// Returns the length of the 2-3 alignment of a TS that starts in `self` and ends in `until`.
    ///
    /// The length is the maximum of the difference of the two sequences.
    pub fn ts_length_until(&self, until: &Self) -> usize {
        let start = self.start();
        let end = until.end();

        (start.ancestor().checked_sub(end.ancestor()).unwrap())
            .max(end.descendant().checked_sub(start.descendant()).unwrap())
    }
}

impl<Cost: Display> Display for SecondaryAnchor<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SA(({}, {}], [{}, {}), {})",
            self.range.offset().ancestor(),
            self.range.limit().ancestor(),
            self.range.offset().descendant(),
            self.range.limit().descendant(),
            self.cost
        )
    }
}

impl<Cost> From<(AnySecondaryAlignmentRange, Cost)> for SecondaryAnchor<Cost> {
    fn from(value: (AnySecondaryAlignmentRange, Cost)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl<Cost: Ord> Ord for SecondaryAnchor<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.range
            .offset()
            .ancestor()
            .min(self.range.offset().descendant())
            .cmp(
                &other
                    .range
                    .offset()
                    .ancestor()
                    .min(other.range.offset().descendant()),
            )
            .then_with(|| {
                self.range
                    .offset()
                    .ancestor()
                    .cmp(&other.range.offset().ancestor())
            })
            .then_with(|| {
                self.range
                    .offset()
                    .descendant()
                    .cmp(&other.range.offset().descendant())
            })
            .then_with(|| self.range.limit().cmp(&other.range.limit()))
            .then_with(|| self.cost.cmp(&other.cost))
    }
}

impl<Cost: Ord> PartialOrd for SecondaryAnchor<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
