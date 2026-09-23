use std::{fmt::Display, ops::Range};

use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::{PrimaryAlignmentCoordinates, range::PrimaryAlignmentRange},
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::secondary::SecondaryAnchor,
};

/// A primary anchor.
///
/// This is an anchor between the two sequences in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its first ordinate and finally by its second ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PrimaryAnchor<Cost> {
    range: PrimaryAlignmentRange,
    cost: Cost,
}

impl<Cost> PrimaryAnchor<Cost> {
    pub fn new(range: PrimaryAlignmentRange, cost: Cost) -> Self {
        Self { range, cost }
    }

    pub fn new_from_ranges(seq1: Range<usize>, seq2: Range<usize>, cost: Cost) -> Self {
        Self::new(PrimaryAlignmentRange::new_from_ranges(seq1, seq2), cost)
    }

    pub fn new_exact(offset: PrimaryAlignmentCoordinates, length: usize) -> Self
    where
        Cost: Zero,
    {
        Self::new(
            PrimaryAlignmentRange::new_equal_length(offset, length),
            Cost::zero(),
        )
    }

    pub fn range(&self) -> PrimaryAlignmentRange {
        self.range
    }

    pub fn start(&self) -> PrimaryAlignmentCoordinates {
        self.range.offset()
    }

    pub fn end(&self) -> PrimaryAlignmentCoordinates {
        self.range.limit()
    }

    pub fn cost(&self) -> Cost
    where
        Cost: Copy,
    {
        self.cost
    }

    pub fn is_at(&self, range: PrimaryAlignmentRange) -> bool {
        self.range == range
    }

    pub fn chaining_gaps(&self, second: &Self) -> Option<(usize, usize)> {
        let gap_start = self.end();
        let gap_end = second.start();
        primary_chaining_gaps(gap_start, gap_end)
    }

    /// Returns the gap between given start coordinates and the start of this anchor.
    pub fn chaining_gaps_from_start(&self, start: PrimaryAlignmentCoordinates) -> (usize, usize)
    where
        Cost: Display,
    {
        let gap_end = self.start();
        primary_chaining_gaps(start, gap_end)
            .unwrap_or_else(|| panic!("self: {self}, start: {start}"))
    }

    /// Returns the gap between the end of this anchor and given end coordinates.
    pub fn chaining_gaps_to_end(&self, end: PrimaryAlignmentCoordinates) -> (usize, usize)
    where
        Cost: Display,
    {
        let gap_start = self.end();
        primary_chaining_gaps(gap_start, end).unwrap_or_else(|| panic!("self: {self}, end: {end}"))
    }

    /// Returns the gap in the descendant for the 12-jump from this anchor to the given anchor.
    pub fn chaining_jump_gap(
        &self,
        second: &SecondaryAnchor<Cost>,
        ts_kind: TsKind,
    ) -> Option<usize> {
        let gap_start = self.end();
        let gap_end = second.start();

        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_start.a(),
            TsDescendant::Seq2 => gap_start.b(),
        };
        let gap_end = gap_end.descendant();

        gap_end.checked_sub(gap_start)
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
}

fn primary_chaining_gaps(
    gap_start: PrimaryAlignmentCoordinates,
    gap_end: PrimaryAlignmentCoordinates,
) -> Option<(usize, usize)> {
    let gap1 = gap_end.a().checked_sub(gap_start.a())?;
    let gap2 = gap_end.b().checked_sub(gap_start.b())?;

    Some((gap1, gap2))
}

impl<Cost: Display> Display for PrimaryAnchor<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PA([{}, {}), [{}, {}), {})",
            self.range.offset().a(),
            self.range.limit().a(),
            self.range.offset().b(),
            self.range.limit().b(),
            self.cost
        )
    }
}

impl<Cost> From<(PrimaryAlignmentRange, Cost)> for PrimaryAnchor<Cost> {
    fn from(value: (PrimaryAlignmentRange, Cost)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl<Cost> From<inexact_alignment_anchors::anchor::Anchor<Cost>> for PrimaryAnchor<Cost> {
    fn from(
        inexact_alignment_anchors::anchor::Anchor {
            offset_a,
            limit_a,
            offset_b,
            limit_b,
            cost,
        }: inexact_alignment_anchors::anchor::Anchor<Cost>,
    ) -> Self {
        Self::new(
            PrimaryAlignmentRange::new_from_ranges(offset_a..limit_a, offset_b..limit_b),
            cost,
        )
    }
}

impl<Cost: Ord> Ord for PrimaryAnchor<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.range
            .offset()
            .a()
            .min(self.range.offset().b())
            .cmp(&other.range.offset().a().min(other.range.offset().b()))
            .then_with(|| self.range.offset().a().cmp(&other.range.offset().a()))
            .then_with(|| self.range.offset().b().cmp(&other.range.offset().b()))
            .then_with(|| self.range.limit().cmp(&other.range.limit()))
            .then_with(|| self.cost.cmp(&other.cost))
    }
}

impl<Cost: Ord> PartialOrd for PrimaryAnchor<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
