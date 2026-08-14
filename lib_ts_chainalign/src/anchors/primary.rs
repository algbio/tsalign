use std::fmt::Display;

use crate::{
    alignment::{
        coordinates::{AlignmentCoordinates, PrimaryAlignmentCoordinates},
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::secondary::SecondaryAnchor,
};

/// A primary anchor.
///
/// This is an anchor between the two sequences in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its first ordinate and finally by its second ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimaryAnchor<Cost> {
    coordinates: PrimaryAlignmentCoordinates,
    cost: Cost,
}

impl<Cost> PrimaryAnchor<Cost> {
    pub fn new(seq1: usize, seq2: usize, cost: Cost) -> Self {
        Self::new_from_start(PrimaryAlignmentCoordinates::new(seq1, seq2), cost)
    }

    pub fn new_from_start(alignment_coordinates: PrimaryAlignmentCoordinates, cost: Cost) -> Self {
        Self {
            coordinates: alignment_coordinates,
            cost,
        }
    }

    pub fn new_from_end(
        alignment_coordinates: &AlignmentCoordinates,
        k: usize,
        cost: Cost,
    ) -> Self {
        Self::new(
            alignment_coordinates
                .primary_ordinate_a()
                .unwrap()
                .checked_sub(k)
                .unwrap(),
            alignment_coordinates
                .primary_ordinate_b()
                .unwrap()
                .checked_sub(k)
                .unwrap(),
            cost,
        )
    }

    pub fn start(&self) -> PrimaryAlignmentCoordinates {
        self.coordinates
    }

    pub fn end(&self, k: usize) -> PrimaryAlignmentCoordinates {
        self.coordinates.increment_both(k)
    }

    pub fn cost(&self) -> Cost
    where
        Cost: Copy,
    {
        self.cost
    }

    pub fn is_at(&self, coordinates: PrimaryAlignmentCoordinates) -> bool {
        self.coordinates == coordinates
    }

    pub fn chaining_gaps(&self, second: &Self, k: usize) -> Option<(usize, usize)> {
        let gap_start = self.end(k);
        let gap_end = second.start();
        primary_chaining_gaps(gap_start, gap_end)
    }

    pub fn chaining_gaps_from_start(&self, start: PrimaryAlignmentCoordinates) -> (usize, usize)
    where
        Cost: Display,
    {
        let gap_end = self.start();
        primary_chaining_gaps(start, gap_end)
            .unwrap_or_else(|| panic!("self: {self}, start: {start}"))
    }

    pub fn chaining_gaps_to_end(&self, end: PrimaryAlignmentCoordinates, k: usize) -> (usize, usize)
    where
        Cost: Display,
    {
        let gap_start = self.end(k);
        primary_chaining_gaps(gap_start, end)
            .unwrap_or_else(|| panic!("self: {self}, end: {end}, k: {k}"))
    }

    /// Returns the gap in the descendant for the 12-jump from this anchor to the given anchor.
    pub fn chaining_jump_gap(
        &self,
        second: &SecondaryAnchor<Cost>,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(k);
        let gap_end = second.start();

        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_start.a(),
            TsDescendant::Seq2 => gap_start.b(),
        };
        let gap_end = gap_end.descendant();

        gap_end.checked_sub(gap_start)
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.coordinates.increment_both(1) == successor.coordinates
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
        write!(f, "PA({}, {})", self.coordinates, self.cost)
    }
}

impl<Cost> From<(usize, usize, Cost)> for PrimaryAnchor<Cost> {
    fn from(value: (usize, usize, Cost)) -> Self {
        Self::new(value.0, value.1, value.2)
    }
}

impl<Cost: Ord> Ord for PrimaryAnchor<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.coordinates
            .a()
            .min(self.coordinates.b())
            .cmp(&other.coordinates.a().min(other.coordinates.b()))
            .then_with(|| self.coordinates.a().cmp(&other.coordinates.a()))
            .then_with(|| self.coordinates.b().cmp(&other.coordinates.b()))
            .then_with(|| self.cost.cmp(&other.cost))
    }
}

impl<Cost: Ord> PartialOrd for PrimaryAnchor<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
