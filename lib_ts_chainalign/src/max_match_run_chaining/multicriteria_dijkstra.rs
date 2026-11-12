//! A prune list for the multicriteria Dijkstra specialised to two criteria.
//!
//! With two criteria, the set of pairwise non-dominated costs forms a "staircase":
//! sorted by increasing first criterion, the second criterion is strictly decreasing.
//! This allows to check for domination with a single binary search,
//! and to find all entries dominated by a new cost as a contiguous range.

use std::{cmp::Reverse, collections::BinaryHeap, hash::Hash};

use multicriteria_dijkstra::{
    node::{MulticriteriaDijkstraCost, MulticriteriaDijkstraNode},
    prune_list::{
        IndexedNode, MulticriteriaDijkstraPruneList, NodeWithIsTarget, NodeWithOptionalPredecessor,
        NodeWithOptionalPredecessorAndIsTarget,
    },
};
use rustc_hash::FxHashMap;

#[cfg(test)]
mod tests;

/// A cost with exactly two criteria, where smaller is better in both.
///
/// The implementation of [`MulticriteriaDijkstraCost::dominates`] must be consistent with
/// the component-wise comparison, i.e. `a.dominates(b)` iff `a.first() <= b.first() && a.second() <= b.second()`.
/// The [`PruneList`] only uses [`first`](Self::first) and [`second`](Self::second).
pub trait TwoCriteriaCost: MulticriteriaDijkstraCost {
    type First: Ord + Copy;
    type Second: Ord + Copy;

    fn first(&self) -> Self::First;

    fn second(&self) -> Self::Second;
}

type Index = usize;

/// A prune list for two-criteria optimisation.
///
/// * Each identifier stores its Pareto front of non-dominated nodes as a staircase,
///   such that domination checks are `O(log n)` and pruning is `O(log n + k)` for `k` pruned nodes.
/// * Target nodes are kept in a global staircase.
///   Pruning of nodes dominated by targets is done lazily:
///   a node is considered pruned if it is dominated by a target node other than itself.
/// * Open nodes are kept in a binary heap ordered lexicographically by cost, with lazy deletion of pruned nodes.
///   This makes the search label-setting in the first criterion, which reduces the number of nodes that get expanded
///   but are later pruned.
/// * All inserted nodes are kept for backtracking, including pruned ones,
///   since they may be predecessors of non-pruned nodes.
pub struct PruneList<Node: MulticriteriaDijkstraNode>
where
    Node::Cost: TwoCriteriaCost,
{
    nodes: Vec<StoredNode<Node>>,
    fronts: FxHashMap<Node::Identifier, ParetoFront<Node::Cost>>,
    targets: ParetoFront<Node::Cost>,
    open: BinaryHeap<Reverse<OpenEntry<Node::Cost>>>,
}

struct StoredNode<Node> {
    node: Node,
    predecessor_index: Option<Index>,
    /// True if the node was removed from the Pareto front of its identifier.
    ///
    /// Domination by targets is not tracked by this flag, but evaluated lazily.
    pruned: bool,
}

type OpenEntry<Cost> = (
    <Cost as TwoCriteriaCost>::First,
    <Cost as TwoCriteriaCost>::Second,
    Index,
);

/// A set of pairwise non-dominated costs.
///
/// Entries are sorted by strictly increasing `first` and strictly decreasing `second`.
struct ParetoFront<Cost: TwoCriteriaCost> {
    entries: Vec<FrontEntry<Cost>>,
}

struct FrontEntry<Cost: TwoCriteriaCost> {
    first: Cost::First,
    second: Cost::Second,
    index: Index,
}

impl<Cost: TwoCriteriaCost> ParetoFront<Cost> {
    /// Returns the entry that dominates the given cost, if any.
    ///
    /// If there are multiple, the one with the largest `first` is returned.
    fn dominator(&self, first: Cost::First, second: Cost::Second) -> Option<&FrontEntry<Cost>> {
        // The entry with the largest `first <= first` has the smallest `second` among all entries with `first <= first`.
        let end = self.entries.partition_point(|entry| entry.first <= first);
        end.checked_sub(1)
            .map(|candidate| &self.entries[candidate])
            .filter(|candidate| candidate.second <= second)
    }

    fn is_dominated(&self, first: Cost::First, second: Cost::Second) -> bool {
        self.dominator(first, second).is_some()
    }

    /// Returns true if the given cost is dominated by an entry with a different index.
    fn is_dominated_by_other(
        &self,
        first: Cost::First,
        second: Cost::Second,
        index: Index,
    ) -> bool {
        // If the dominator is the entry itself, then all entries with smaller `first` have a larger `second`,
        // so no other entry dominates it.
        self.dominator(first, second)
            .is_some_and(|dominator| dominator.index != index)
    }

    /// Inserts the given cost and removes all entries dominated by it.
    /// Calls `on_removed` with the index of each removed entry.
    ///
    /// The cost must not be dominated by any entry.
    fn insert(
        &mut self,
        first: Cost::First,
        second: Cost::Second,
        index: Index,
        mut on_removed: impl FnMut(Index),
    ) {
        debug_assert!(!self.is_dominated(first, second));

        // Dominated entries have `first >= first` and `second >= second`.
        // Since `second` is decreasing, they form a prefix of the entries with `first >= first`.
        let start = self.entries.partition_point(|entry| entry.first < first);
        let end = start + self.entries[start..].partition_point(|entry| entry.second >= second);

        if start == end {
            self.entries.insert(
                start,
                FrontEntry {
                    first,
                    second,
                    index,
                },
            );
        } else {
            for entry in &self.entries[start..end] {
                on_removed(entry.index);
            }
            self.entries[start] = FrontEntry {
                first,
                second,
                index,
            };
            self.entries.drain(start + 1..end);
        }
    }

    fn iter(&self) -> impl Iterator<Item = &FrontEntry<Cost>> {
        self.entries.iter()
    }

    fn clear(&mut self) {
        self.entries.clear();
    }
}

impl<Node: MulticriteriaDijkstraNode> PruneList<Node>
where
    Node::Cost: TwoCriteriaCost,
{
    /// Returns true if the node at the given index is dominated by a target node other than itself.
    fn is_dominated_by_target(&self, index: Index) -> bool {
        let cost = self.nodes[index].node.cost();
        self.targets
            .is_dominated_by_other(cost.first(), cost.second(), index)
    }

    /// Returns true if the given front entry is pruned by a target node other than itself.
    fn is_entry_dominated_by_target(&self, entry: &FrontEntry<Node::Cost>) -> bool {
        self.targets
            .is_dominated_by_other(entry.first, entry.second, entry.index)
    }

    /// Iterate over the Pareto fronts of all identifiers.
    ///
    /// The nodes of each front are ordered by increasing first criterion and decreasing second criterion.
    /// Nodes dominated by a target node other than themselves are skipped.
    pub fn iter_fronts(
        &self,
    ) -> impl Iterator<Item = (&Node::Identifier, impl Iterator<Item = &Node>)> {
        self.fronts.iter().map(|(identifier, front)| {
            (
                identifier,
                front
                    .iter()
                    .filter(|entry| !self.is_entry_dominated_by_target(entry))
                    .map(|entry| &self.nodes[entry.index].node),
            )
        })
    }

    /// Returns the node with the smallest first criterion at the given identifier.
    ///
    /// Nodes dominated by a target node other than themselves are skipped.
    pub fn min_first_node_at(&self, identifier: &Node::Identifier) -> Option<&Node>
    where
        Node::Identifier: Hash,
    {
        self.fronts
            .get(identifier)?
            .iter()
            .find(|entry| !self.is_entry_dominated_by_target(entry))
            .map(|entry| &self.nodes[entry.index].node)
    }
}

impl<Node: MulticriteriaDijkstraNode> MulticriteriaDijkstraPruneList for PruneList<Node>
where
    Node::Cost: TwoCriteriaCost,
    Node::Identifier: Hash + Clone,
{
    type Node = Node;

    type Index = Index;

    fn insert(
        &mut self,
        NodeWithOptionalPredecessorAndIsTarget {
            node,
            predecessor_index,
            is_target,
        }: NodeWithOptionalPredecessorAndIsTarget<Self::Node, Self::Index>,
    ) -> Option<Self::Index> {
        let first = node.cost().first();
        let second = node.cost().second();

        // Check if node is dominated.
        if self.targets.is_dominated(first, second)
            || self
                .fronts
                .get(node.identifier())
                .is_some_and(|front| front.is_dominated(first, second))
        {
            return None;
        }

        let index = self.nodes.len();

        // Insert into the front of the identifier, pruning what it dominates.
        let nodes = &mut self.nodes;
        self.fronts
            .entry(node.identifier().clone())
            .or_insert_with(|| ParetoFront {
                entries: Vec::new(),
            })
            .insert(first, second, index, |removed| {
                nodes[removed].pruned = true;
            });

        // Open nodes at other identifiers that are dominated by a target are pruned lazily.
        if is_target {
            self.targets.insert(first, second, index, |_| {});
        }

        self.nodes.push(StoredNode {
            node,
            predecessor_index,
            pruned: false,
        });
        self.open.push(Reverse((first, second, index)));

        Some(index)
    }

    fn next(&mut self) -> Option<IndexedNode<&Self::Node, Self::Index>> {
        while let Some(Reverse((_, _, index))) = self.open.pop() {
            if self.nodes[index].pruned || self.is_dominated_by_target(index) {
                continue;
            }

            return Some(IndexedNode {
                node: &self.nodes[index].node,
                index,
            });
        }

        None
    }

    fn get(
        &self,
        index: Self::Index,
    ) -> Option<NodeWithOptionalPredecessor<&Self::Node, Self::Index>> {
        // Pruned nodes are returned as well, since they may be predecessors of non-pruned nodes.
        self.nodes
            .get(index)
            .map(|node| NodeWithOptionalPredecessor {
                node: &node.node,
                predecessor_index: node.predecessor_index,
            })
    }

    fn iter_nodes_with_predecessor_at_identifier<'this, 'identifier>(
        &'this self,
        identifier: &'identifier <Self::Node as MulticriteriaDijkstraNode>::Identifier,
    ) -> impl use<'this, Node>
    + Iterator<Item = NodeWithOptionalPredecessor<&'this Self::Node, Self::Index>>
    where
        Self::Node: 'this,
    {
        self.fronts
            .get(identifier)
            .into_iter()
            .flat_map(ParetoFront::iter)
            .filter(|entry| {
                !self
                    .targets
                    .is_dominated_by_other(entry.first, entry.second, entry.index)
            })
            .map(|entry| {
                let node = &self.nodes[entry.index];
                NodeWithOptionalPredecessor {
                    node: &node.node,
                    predecessor_index: node.predecessor_index,
                }
            })
    }

    fn clear(&mut self) {
        self.nodes.clear();
        self.fronts.clear();
        self.targets.clear();
        self.open.clear();
    }

    fn supports_backtracking(&self) -> bool {
        true
    }
}

impl<Node: MulticriteriaDijkstraNode> Default for PruneList<Node>
where
    Node::Cost: TwoCriteriaCost,
{
    fn default() -> Self {
        Self {
            nodes: Vec::new(),
            fronts: FxHashMap::default(),
            targets: ParetoFront {
                entries: Vec::new(),
            },
            open: BinaryHeap::new(),
        }
    }
}

impl<Node: MulticriteriaDijkstraNode> Extend<NodeWithOptionalPredecessorAndIsTarget<Node, Index>>
    for PruneList<Node>
where
    Node::Cost: TwoCriteriaCost,
    Node::Identifier: Hash + Clone,
{
    fn extend<T: IntoIterator<Item = NodeWithOptionalPredecessorAndIsTarget<Node, Index>>>(
        &mut self,
        iter: T,
    ) {
        for node in iter {
            self.insert(node);
        }
    }
}

impl<Node: MulticriteriaDijkstraNode> Extend<NodeWithIsTarget<Node>> for PruneList<Node>
where
    Node::Cost: TwoCriteriaCost,
    Node::Identifier: Hash + Clone,
{
    fn extend<T: IntoIterator<Item = NodeWithIsTarget<Node>>>(&mut self, iter: T) {
        for NodeWithIsTarget { node, is_target } in iter {
            self.insert(NodeWithOptionalPredecessorAndIsTarget {
                node,
                predecessor_index: None,
                is_target,
            });
        }
    }
}
