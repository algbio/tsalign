use std::{cmp::Ordering, iter::Peekable};

use crate::{
    alignment::{
        coordinates::{AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates},
        ts_kind::TsKind,
    },
    anchors::{Anchors, index::AnchorIndex, primary::PrimaryAnchor, secondary::SecondaryAnchor},
};

struct PrimaryAnchorToIndexIter<CoordinateIter: Iterator, AnchorIter: Iterator> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

struct SecondaryAnchorToIndexIter<CoordinateIter: Iterator, AnchorIter: Iterator> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

pub trait PartialIntoAnchorIndex {
    type IntoPartSource;
    type IntoTarget;

    fn source_part(&self) -> &Self::IntoPartSource;

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget;
}

impl<
    Cost: Ord,
    Coordinates: PartialIntoAnchorIndex<IntoPartSource = PrimaryAlignmentCoordinates> + Clone,
    CoordinateIter: Iterator<Item = Coordinates>,
    AnchorIter: Iterator<Item = (AnchorIndex, PrimaryAnchor<Cost>)>,
> Iterator for PrimaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Coordinates::IntoTarget;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.source_part().cmp(&anchor.start()) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                }
                Ordering::Equal => {
                    // Do not advance the coordinate iterator, since further anchors may start at the same coordinates.
                    // Duplicates of the coordinates (with higher costs) are skipped by the `Less` branch
                    // once the anchors have advanced beyond them.
                    let result = Some(coordinate_anchor.clone().partial_into(*anchor_index));
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        None
    }
}

impl<
    Cost: Ord,
    Coordinates: PartialIntoAnchorIndex<IntoPartSource = AnySecondaryAlignmentCoordinates> + Clone,
    CoordinateIter: Iterator<Item = Coordinates>,
    AnchorIter: Iterator<Item = (AnchorIndex, SecondaryAnchor<Cost>)>,
> Iterator for SecondaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Coordinates::IntoTarget;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.source_part().cmp(&anchor.start()) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                }
                Ordering::Equal => {
                    // Do not advance the coordinate iterator, since further anchors may start at the same coordinates.
                    // Duplicates of the coordinates (with higher costs) are skipped by the `Less` branch
                    // once the anchors have advanced beyond them.
                    let result = Some(coordinate_anchor.clone().partial_into(*anchor_index));
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        None
    }
}

impl<Cost> Anchors<Cost> {
    /// Returns an iterator over the primary anchors that start at the given primary alignment coordinates.
    ///
    /// For each primary anchor whose start is among the given coordinates, the iterator returns
    /// the first item with these coordinates, converted with [`PartialIntoAnchorIndex::partial_into`].
    /// If several anchors start at the same coordinates, then the same item is returned for each of them,
    /// and if the coordinates occur several times, then only the first occurrence is used.
    /// Anchors whose start is not among the given coordinates and coordinates that are not the start of an anchor are skipped.
    ///
    /// The items must be sorted by their coordinates.
    /// If the items are pairs of coordinates and costs (e.g. the targets reached by a gap aligner, which may
    /// report the same coordinates several times with different costs), then they must be sorted by coordinates and then by cost,
    /// such that each anchor gets the minimum cost of its start coordinates.
    pub fn primary_anchor_to_index_iter<
        Coordinates: PartialIntoAnchorIndex<IntoPartSource = PrimaryAlignmentCoordinates> + Clone,
    >(
        &self,
        iter: impl IntoIterator<Item = Coordinates>,
    ) -> impl Iterator<Item = Coordinates::IntoTarget>
    where
        Cost: Ord + Copy,
    {
        PrimaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_primaries().peekable(),
        }
    }

    /// Returns an iterator over the secondary anchors of the given kind that start at the given secondary alignment coordinates.
    ///
    /// Works like [`Self::primary_anchor_to_index_iter`], see there for the requirements on the order of the items.
    pub fn secondary_anchor_to_index_iter<
        Coordinates: PartialIntoAnchorIndex<IntoPartSource = AnySecondaryAlignmentCoordinates> + Clone,
    >(
        &self,
        iter: impl IntoIterator<Item = Coordinates>,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = Coordinates::IntoTarget>
    where
        Cost: Ord + Copy,
    {
        SecondaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_secondaries(ts_kind).peekable(),
        }
    }
}

impl PartialIntoAnchorIndex for PrimaryAlignmentCoordinates {
    type IntoPartSource = Self;

    type IntoTarget = AnchorIndex;

    fn source_part(&self) -> &Self::IntoPartSource {
        self
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        target
    }
}

impl PartialIntoAnchorIndex for AnySecondaryAlignmentCoordinates {
    type IntoPartSource = Self;

    type IntoTarget = AnchorIndex;

    fn source_part(&self) -> &Self::IntoPartSource {
        self
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        target
    }
}

impl<T> PartialIntoAnchorIndex for (PrimaryAlignmentCoordinates, T) {
    type IntoPartSource = PrimaryAlignmentCoordinates;

    type IntoTarget = (AnchorIndex, T);

    fn source_part(&self) -> &Self::IntoPartSource {
        &self.0
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        (target, self.1)
    }
}

impl<T> PartialIntoAnchorIndex for (AnySecondaryAlignmentCoordinates, T) {
    type IntoPartSource = AnySecondaryAlignmentCoordinates;

    type IntoTarget = (AnchorIndex, T);

    fn source_part(&self) -> &Self::IntoPartSource {
        &self.0
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        (target, self.1)
    }
}
