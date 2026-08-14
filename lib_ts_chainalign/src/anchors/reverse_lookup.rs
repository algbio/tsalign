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
    Coordinates: PartialIntoAnchorIndex<IntoPartSource = PrimaryAlignmentCoordinates>,
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
                    let result = Some(
                        self.coordinate_iter
                            .next()
                            .unwrap()
                            .partial_into(*anchor_index),
                    );
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
    Coordinates: PartialIntoAnchorIndex<IntoPartSource = AnySecondaryAlignmentCoordinates>,
    CoordinateIter: Iterator<Item = Coordinates>,
    AnchorIter: Iterator<Item = (AnchorIndex, SecondaryAnchor<Cost>)>,
> Iterator for SecondaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Coordinates::IntoTarget;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.source_part().cmp(anchor) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                }
                Ordering::Equal => {
                    let result = Some(
                        self.coordinate_iter
                            .next()
                            .unwrap()
                            .partial_into(*anchor_index),
                    );
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
    /// Returns an iterator over the primary anchor indices that correspond to the given primary alignment coordinates.
    ///
    /// If a pair of primary alignment coordinates does not correspond to a primary anchor, then the iterator returns `Some(None)`.
    pub fn primary_anchor_to_index_iter<
        Coordinates: PartialIntoAnchorIndex<IntoPartSource = PrimaryAlignmentCoordinates>,
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

    /// Returns an iterator over the secondary anchor indices that correspond to the given secondary alignment coordinates.
    ///
    /// If a pair of secondary alignment coordinates does not correspond to a secondary anchor, then the iterator returns `Some(None)`.
    pub fn secondary_anchor_to_index_iter<
        Coordinates: PartialIntoAnchorIndex<IntoPartSource = AnySecondaryAlignmentCoordinates>,
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
