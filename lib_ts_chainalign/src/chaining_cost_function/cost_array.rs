use std::ops::{Index, IndexMut};

use bitvec::{bitvec, order::LocalBits, vec::BitVec};
use num_traits::Zero;

use crate::anchors::index::AnchorIndex;

pub struct ChainingCostArray<Cost> {
    len: [AnchorIndex; 2],
    data: Vec<Cost>,
    /// True if the anchor indices of the two dimensions map to the same set of anchors.
    maps_between_same_anchors: bool,
    /// For each ordinate 0, store the indices of the corresponding ordinates 1, ordered by their cost.
    cost_order_permutation: Vec<Vec<AnchorIndex>>,
    /// Used to mark chaining costs as exact, i.e. once the cost between two anchors is exactly computed via alignment, it is marked as exact and will not be recomputed.
    is_exact: BitVec<usize, LocalBits>,
}

impl<Cost> ChainingCostArray<Cost> {
    pub fn new_from_cost(
        dim: [impl Into<AnchorIndex>; 2],
        cost: Cost,
        maps_between_same_anchors: bool,
    ) -> Self
    where
        Cost: Clone,
    {
        let dim = dim.map(Into::into);

        Self {
            len: dim,
            data: vec![cost; dim[0].as_usize() * dim[1].as_usize()],
            maps_between_same_anchors,
            cost_order_permutation: vec![Vec::new(); dim[0].as_usize()],
            is_exact: bitvec![usize, LocalBits; 0; dim[0].as_usize() * dim[1].as_usize()],
        }
    }

    pub fn is_exact(&self, c1: AnchorIndex, c2: AnchorIndex) -> bool {
        self.is_exact[coordinates_to_index(c1, c2, self.len)]
    }

    pub fn set_exact(&mut self, c1: AnchorIndex, c2: AnchorIndex) {
        self.is_exact
            .set(coordinates_to_index(c1, c2, self.len), true);
    }

    pub fn dim(&self) -> (AnchorIndex, AnchorIndex) {
        (self.len[0], self.len[1])
    }

    fn compute_cost_order_permutation_if_missing(&mut self, c1: AnchorIndex)
    where
        Cost: Copy + Ord,
    {
        if self.cost_order_permutation[c1.as_usize()].is_empty() {
            if self.maps_between_same_anchors {
                self.cost_order_permutation[c1.as_usize()].extend(
                    ((0..c1.as_usize()).chain(c1.as_usize() + 1..))
                        .take(self.len[1].as_usize().saturating_sub(1))
                        .map(AnchorIndex::from),
                );
            } else {
                self.cost_order_permutation[c1.as_usize()]
                    .extend((0..self.len[1].as_usize()).map(AnchorIndex::from));
            }
            self.cost_order_permutation[c1.as_usize()]
                .sort_unstable_by_key(|c2| self.data[coordinates_to_index(c1, *c2, self.len)]);
        }
    }

    /// Iterate over all ordinates `c2` that belong to this `c1` in order of their cost.
    #[expect(dead_code)]
    pub fn iter_in_cost_order(
        &mut self,
        c1: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.iter_in_cost_order_from(c1, AnchorIndex::zero())
    }

    pub fn iter_in_cost_order_from(
        &mut self,
        c1: AnchorIndex,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.compute_cost_order_permutation_if_missing(c1);
        let data = &self.data;
        let len = self.len;
        self.cost_order_permutation[c1.as_usize()]
            .iter()
            .copied()
            .skip(offset.as_usize())
            .map(move |c2| (c2, data[coordinates_to_index(c1, c2, len)]))
    }
}

impl<Cost> Index<[AnchorIndex; 2]> for ChainingCostArray<Cost> {
    type Output = Cost;

    fn index(&self, index: [AnchorIndex; 2]) -> &Self::Output {
        &self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

impl<Cost> IndexMut<[AnchorIndex; 2]> for ChainingCostArray<Cost> {
    fn index_mut(&mut self, index: [AnchorIndex; 2]) -> &mut Self::Output {
        self.cost_order_permutation[index[0].as_usize()].clear();
        &mut self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

#[inline]
fn coordinates_to_index(c1: AnchorIndex, c2: AnchorIndex, len: [AnchorIndex; 2]) -> usize {
    c1.as_usize() * len[1].as_usize() + c2.as_usize()
}
