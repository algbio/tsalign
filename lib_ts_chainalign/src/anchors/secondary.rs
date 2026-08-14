use std::fmt::Display;

use crate::{
    alignment::{
        coordinates::{AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates},
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::primary::PrimaryAnchor,
};

/// A secondary anchor.
///
/// This is an anchor between the ancestor in reverse direction and the descendant in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its ancestor ordinate and finally by its descendant ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SecondaryAnchor<Cost> {
    coordinates: AnySecondaryAlignmentCoordinates,
    pub(super) cost: Cost,
}

impl<Cost> SecondaryAnchor<Cost> {
    pub fn new(ancestor: usize, descendant: usize, cost: Cost) -> Self {
        Self::new_from_start(
            AnySecondaryAlignmentCoordinates::new(ancestor, descendant),
            cost,
        )
    }

    pub fn new_from_start(
        alignment_coordinates: AnySecondaryAlignmentCoordinates,
        cost: Cost,
    ) -> Self {
        Self {
            coordinates: alignment_coordinates,
            cost,
        }
    }

    pub fn start(&self) -> AnySecondaryAlignmentCoordinates {
        self.coordinates
    }

    pub fn end(&self, k: usize) -> AnySecondaryAlignmentCoordinates {
        self.coordinates.increment_both(k)
    }

    pub fn cost(&self) -> Cost
    where
        Cost: Copy,
    {
        self.cost
    }

    /// Returns true if the anchor is at the given coordinates.
    ///
    /// Does not check the `ts_kind`, and will produce false positives if the coordinates given have the wrong `ts_kind`.
    pub fn is_at(&self, coordinates: AnySecondaryAlignmentCoordinates) -> bool {
        self.coordinates == coordinates
    }

    pub fn chaining_gaps(&self, second: &Self, k: usize) -> Option<(usize, usize)> {
        let gap_start = self.end(k);
        let gap_end = second.start();

        let gap1 = gap_start.ancestor().checked_sub(gap_end.ancestor())?;
        let gap2 = gap_end.descendant().checked_sub(gap_start.descendant())?;

        Some((gap1, gap2))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &PrimaryAnchor<Cost>,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(k);
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
    ) -> usize {
        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => start.a(),
            TsDescendant::Seq2 => start.b(),
        };
        let gap_end = self.start().descendant();

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn chaining_jump_gap_to_end(
        &self,
        end: PrimaryAlignmentCoordinates,
        ts_kind: TsKind,
        k: usize,
    ) -> usize {
        let gap_start = self.end(k).descendant();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => end.a(),
            TsDescendant::Seq2 => end.b(),
        };

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.coordinates.increment_both(1) == successor.coordinates
    }

    /// Returns the length of the 2-3 alignment of a TS that starts in `self` and ends in `until`.
    ///
    /// The length is the maximum of the difference of the two sequences.
    pub fn ts_length_until(&self, until: &Self, k: usize) -> usize {
        let start = self.start();
        let end = until.end(k);

        (start.ancestor().checked_sub(end.ancestor()).unwrap())
            .max(end.descendant().checked_sub(start.descendant()).unwrap())
    }
}

impl<Cost: Display> Display for SecondaryAnchor<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SA({}, {})", self.coordinates, self.cost)
    }
}

impl<Cost> From<(usize, usize, Cost)> for SecondaryAnchor<Cost> {
    fn from(value: (usize, usize, Cost)) -> Self {
        Self::new(value.0, value.1, value.2)
    }
}

impl<Cost: Ord> Ord for SecondaryAnchor<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.coordinates
            .ancestor()
            .min(self.coordinates.descendant())
            .cmp(
                &other
                    .coordinates
                    .ancestor()
                    .min(other.coordinates.descendant()),
            )
            .then_with(|| {
                self.coordinates
                    .ancestor()
                    .cmp(&other.coordinates.ancestor())
            })
            .then_with(|| {
                self.coordinates
                    .descendant()
                    .cmp(&other.coordinates.descendant())
            })
            .then_with(|| self.cost.cmp(&other.cost))
    }
}

impl<Cost: Ord> PartialOrd for SecondaryAnchor<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
