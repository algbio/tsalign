use std::{cmp::Ordering, iter};

use generic_a_star::{
    AStarNode,
    arena::{ArenaIndex, NodeArena, NodeIndex},
    closed_lists::AStarClosedList,
    cost::AStarCost,
    reset::Reset,
};
use itertools::Itertools;
use num_traits::CheckedSub;
use rustc_hash::{FxHashMapSeed, FxSeededState};

use crate::{
    anchors::index::AnchorIndex,
    chain_align::chainer::{Identifier, Node},
};

/// A closed list for the chainer that keeps one node per identifier without its offset.
///
/// The nodes are small, so they are moved out of the [`NodeArena`] of the search and stored here directly.
/// This avoids looking them up in the arena when checking successors.
pub struct ChainerClosedList<Cost> {
    start: Option<Node<Cost>>,
    start_to_primary: Option<Node<Cost>>,
    start_to_secondary: [Option<Node<Cost>>; 4],

    primary_to_primary: Vec<Option<Node<Cost>>>,
    primary_to_secondary: [Vec<Option<Node<Cost>>>; 4],

    secondary_to_secondary: [FxHashMapSeed<(AnchorIndex, AnchorIndex), Node<Cost>>; 4],
    secondary_to_primary: [FxHashMapSeed<(AnchorIndex, AnchorIndex), Node<Cost>>; 4],

    end: Option<Node<Cost>>,

    len: usize,
}

impl<Cost: AStarCost> ChainerClosedList<Cost> {
    fn new_list() -> Self {
        Self {
            start: None,
            start_to_primary: None,
            start_to_secondary: [None; 4],

            primary_to_primary: Vec::new(),
            primary_to_secondary: iter::repeat_n(Vec::new(), 4)
                .collect_vec()
                .try_into()
                .unwrap(),

            secondary_to_secondary: iter::repeat_n(
                FxHashMapSeed::with_hasher(FxSeededState::with_seed(0)),
                4,
            )
            .collect_vec()
            .try_into()
            .unwrap(),
            secondary_to_primary: iter::repeat_n(
                FxHashMapSeed::with_hasher(FxSeededState::with_seed(0)),
                4,
            )
            .collect_vec()
            .try_into()
            .unwrap(),

            end: None,

            len: 0,
        }
    }

    /// Inserts the node if it is better than the node with the same identifier, if any.
    ///
    /// Returns the node that is not kept, if there was a node with the same identifier.
    fn insert_node(&mut self, identifier: Identifier, node: Node<Cost>) -> Option<Node<Cost>> {
        match identifier {
            Identifier::Start => {
                debug_assert!(node.cost.is_zero());
                debug_assert!(self.start.is_none());
                self.start = Some(node);
                self.len += 1;
                None
            }
            Identifier::StartToPrimary { offset } => {
                let current = &mut self.start_to_primary;
                if let Some(current) = current.as_mut() {
                    if Ordering::Less
                        == node
                            .offset_zero_cost
                            .cmp(&current.offset_zero_cost)
                            .then(current.identifier.offset().cmp(&offset))
                    {
                        // Insert better node.
                        let previous = *current;
                        *current = node;
                        Some(previous)
                    } else {
                        // Keep current node.
                        Some(node)
                    }
                } else {
                    *current = Some(node);
                    self.len += 1;
                    None
                }
            }
            Identifier::StartToSecondary { ts_kind, offset } => {
                let current = &mut self.start_to_secondary[ts_kind.index()];
                if let Some(current) = current.as_mut() {
                    if Ordering::Less
                        == node
                            .offset_zero_cost
                            .cmp(&current.offset_zero_cost)
                            .then(current.identifier.offset().cmp(&offset))
                    {
                        // Insert better node.
                        let previous = *current;
                        *current = node;
                        Some(previous)
                    } else {
                        // Keep current node.
                        Some(node)
                    }
                } else {
                    *current = Some(node);
                    self.len += 1;
                    None
                }
            }
            Identifier::PrimaryToPrimary { index, offset } => {
                if let Some(missing_length) =
                    (index + 1).checked_sub(&self.primary_to_primary.len().into())
                {
                    self.primary_to_primary
                        .extend(iter::repeat_n(None, missing_length.as_usize()));
                }

                let current = &mut self.primary_to_primary[index.as_usize()];
                if let Some(current) = current.as_mut() {
                    if Ordering::Less
                        == node
                            .offset_zero_cost
                            .cmp(&current.offset_zero_cost)
                            .then(current.identifier.offset().cmp(&offset))
                    {
                        // Insert better node.
                        let previous = *current;
                        *current = node;
                        Some(previous)
                    } else {
                        // Keep current node.
                        Some(node)
                    }
                } else {
                    *current = Some(node);
                    self.len += 1;
                    None
                }
            }
            Identifier::PrimaryToSecondary {
                index,
                ts_kind,
                offset,
            } => {
                if let Some(missing_length) = (index + 1)
                    .checked_sub(&self.primary_to_secondary[ts_kind.index()].len().into())
                {
                    self.primary_to_secondary[ts_kind.index()]
                        .extend(iter::repeat_n(None, missing_length.as_usize()));
                }

                let current = &mut self.primary_to_secondary[ts_kind.index()][index.as_usize()];
                if let Some(current) = current.as_mut() {
                    if Ordering::Less
                        == node
                            .offset_zero_cost
                            .cmp(&current.offset_zero_cost)
                            .then(current.identifier.offset().cmp(&offset))
                    {
                        // Insert better node.
                        let previous = *current;
                        *current = node;
                        Some(previous)
                    } else {
                        // Keep current node.
                        Some(node)
                    }
                } else {
                    *current = Some(node);
                    self.len += 1;
                    None
                }
            }
            Identifier::SecondaryToSecondary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => {
                let mut result = None;
                self.secondary_to_secondary[ts_kind.index()]
                    .entry((index, first_secondary_index))
                    .and_modify(|current| {
                        if Ordering::Less
                            == node
                                .offset_zero_cost
                                .cmp(&current.offset_zero_cost)
                                .then(current.identifier.offset().cmp(&offset))
                        {
                            // Insert better node.
                            result = Some(*current);
                            *current = node;
                        } else {
                            // Keep current node.
                            result = Some(node);
                        }
                    })
                    .or_insert_with(|| {
                        self.len += 1;
                        node
                    });

                result
            }
            Identifier::SecondaryToPrimary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => {
                let mut result = None;
                self.secondary_to_primary[ts_kind.index()]
                    .entry((index, first_secondary_index))
                    .and_modify(|current| {
                        if Ordering::Less
                            == node
                                .offset_zero_cost
                                .cmp(&current.offset_zero_cost)
                                .then(current.identifier.offset().cmp(&offset))
                        {
                            // Insert better node.
                            result = Some(*current);
                            *current = node;
                        } else {
                            // Keep current node.
                            result = Some(node);
                        }
                    })
                    .or_insert_with(|| {
                        self.len += 1;
                        node
                    });

                result
            }
            Identifier::End => {
                let current = &mut self.end;
                if let Some(current) = current.as_mut() {
                    if node.offset_zero_cost < current.offset_zero_cost {
                        // Insert better node.
                        let previous = *current;
                        *current = node;
                        Some(previous)
                    } else {
                        // Keep current node.
                        Some(node)
                    }
                } else {
                    *current = Some(node);
                    self.len += 1;
                    None
                }
            }
        }
    }

    fn get_node(&self, identifier: &Identifier) -> Option<&Node<Cost>> {
        match identifier {
            Identifier::Start => self.start.as_ref(),
            Identifier::StartToPrimary { offset } => {
                (identifier.offset() <= *offset).then_some(self.start_to_primary.as_ref())?
            }
            Identifier::StartToSecondary { ts_kind, offset } => (identifier.offset() <= *offset)
                .then_some(self.start_to_secondary[ts_kind.index()].as_ref())?,
            Identifier::PrimaryToPrimary { index, offset } => (identifier.offset() <= *offset)
                .then_some(self.primary_to_primary.get(index.as_usize())?.as_ref())?,
            Identifier::PrimaryToSecondary {
                index,
                ts_kind,
                offset,
            } => (identifier.offset() <= *offset).then_some(
                self.primary_to_secondary[ts_kind.index()]
                    .get(index.as_usize())?
                    .as_ref(),
            )?,
            Identifier::SecondaryToSecondary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => (identifier.offset() <= *offset).then_some(
                self.secondary_to_secondary[ts_kind.index()].get(&(*index, *first_secondary_index)),
            )?,
            Identifier::SecondaryToPrimary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => (identifier.offset() <= *offset).then_some(
                self.secondary_to_primary[ts_kind.index()].get(&(*index, *first_secondary_index)),
            )?,
            Identifier::End => self.end.as_ref(),
        }
    }

    fn iter_nodes(&self) -> impl Iterator<Item = &Node<Cost>> {
        self.start
            .iter()
            .chain(&self.start_to_primary)
            .chain(
                self.start_to_secondary
                    .iter()
                    .flat_map(IntoIterator::into_iter),
            )
            .chain(self.primary_to_primary.iter().flatten())
            .chain(
                self.primary_to_secondary
                    .iter()
                    .flat_map(IntoIterator::into_iter)
                    .flatten(),
            )
            .chain(
                self.secondary_to_secondary
                    .iter()
                    .flat_map(FxHashMapSeed::values),
            )
            .chain(
                self.secondary_to_primary
                    .iter()
                    .flat_map(FxHashMapSeed::values),
            )
            .chain(&self.end)
    }

    fn can_skip(&self, node: &Node<Cost>) -> bool {
        if let Some(current) = self.get_node(node.identifier()) {
            Ordering::Less
                != node
                    .offset_zero_cost
                    .cmp(&current.offset_zero_cost)
                    .then(current.identifier.offset().cmp(&node.identifier.offset()))
        } else {
            false
        }
    }
}

impl<Cost: AStarCost, I: ArenaIndex> AStarClosedList<Node<Cost>, I> for ChainerClosedList<Cost> {
    fn new() -> Self {
        Self::new_list()
    }

    fn len(&self) -> usize {
        self.len
    }

    fn insert(&mut self, index: NodeIndex<I>, nodes: &mut NodeArena<Node<Cost>, I>) -> bool {
        let node = nodes.remove(index);
        self.insert_node(node.identifier, node).is_some()
    }

    fn get<'nodes>(
        &'nodes self,
        identifier: &Identifier,
        _nodes: &'nodes NodeArena<Node<Cost>, I>,
    ) -> Option<&'nodes Node<Cost>> {
        self.get_node(identifier)
    }

    fn iter<'nodes>(
        &'nodes self,
        _nodes: &'nodes NodeArena<Node<Cost>, I>,
    ) -> impl Iterator<Item = &'nodes Node<Cost>>
    where
        Node<Cost>: 'nodes,
    {
        self.iter_nodes()
    }

    fn can_skip_node(
        &self,
        node: &Node<Cost>,
        _nodes: &NodeArena<Node<Cost>, I>,
        _is_label_setting: bool,
    ) -> bool {
        self.can_skip(node)
    }
}

impl<Cost> Reset for ChainerClosedList<Cost> {
    fn reset(&mut self) {
        self.start.reset();
        self.start_to_primary.reset();
        self.start_to_secondary.reset();

        self.primary_to_primary.clear();
        self.primary_to_secondary.iter_mut().for_each(Vec::clear);

        self.secondary_to_secondary.reset();
        self.secondary_to_primary.reset();

        self.end.reset();

        self.len = 0;
    }
}
