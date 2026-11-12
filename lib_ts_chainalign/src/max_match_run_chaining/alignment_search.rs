//! Shared parts of the aligners that are based on the multicriteria Dijkstra.
//!
//! The aligners minimise the alignment cost and the length of the current match run.
//! Target nodes are not reported to the prune list, since target pruning would also prune nodes with
//! the same alignment cost as a target, but a larger match run.
//! Instead, nodes whose alignment cost exceeds the best target cost (or the cost limit) are not inserted.

use std::{fmt::Display, hash::Hash};

use generic_a_star::cost::AStarCost;
use multicriteria_dijkstra::{
    MulticriteriaDijkstra, MulticriteriaDijkstraContext, MulticriteriaDijkstraDatastructures,
    node::{MulticriteriaDijkstraCost, MulticriteriaDijkstraNode},
    prune_list::{
        MulticriteriaDijkstraPruneList, NodeWithOptionalPredecessor,
        NodeWithOptionalPredecessorAndIsTarget,
    },
};
use num_traits::bounds::UpperBounded;

use crate::{
    alignment::{Alignment, AlignmentType},
    max_match_run_chaining::multicriteria_dijkstra::{PruneList, TwoCriteriaCost},
};

/// The cost of a node, consisting of the alignment cost and the length of the current run of matches.
///
/// A node dominates another node if it is not worse in both criteria.
/// This is sound, because a smaller match run allows to extend the alignment with at least the same operations.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NodeCost<Cost> {
    pub cost: Cost,
    pub match_run: u32,
}

/// A node of an aligner.
pub trait AlignmentNode:
    MulticriteriaDijkstraNode<Cost = NodeCost<Self::AlignmentCost>, Identifier: Hash + Clone>
{
    type AlignmentCost: AStarCost;

    /// The alignment operation from the predecessor to this node.
    ///
    /// This is `None` for root nodes and for nodes that are connected to their predecessor without an alignment operation.
    fn predecessor_alignment_type(&self) -> Option<AlignmentType>;

    fn alignment_cost(&self) -> Self::AlignmentCost {
        self.cost().cost
    }
}

/// The context of an aligner.
pub trait AlignmentContext: MulticriteriaDijkstraContext<Node: AlignmentNode> {
    /// Returns all identifiers that are targets.
    fn target_identifiers(
        &self,
    ) -> impl Iterator<Item = <Self::Node as MulticriteriaDijkstraNode>::Identifier>;
}

pub type Datastructures<Node> = MulticriteriaDijkstraDatastructures<PruneList<Node>>;
pub type Dijkstra<Context> =
    MulticriteriaDijkstra<PruneList<<Context as MulticriteriaDijkstraContext>::Node>, Context>;
type AlignmentCostOf<Context> =
    <<Context as MulticriteriaDijkstraContext>::Node as AlignmentNode>::AlignmentCost;

/// Creates a Dijkstra that reuses the given datastructures.
///
/// The datastructures must be returned with [`MulticriteriaDijkstra::into_datastructures`].
pub fn new_dijkstra<Context: AlignmentContext>(
    datastructures: &mut Option<Datastructures<Context::Node>>,
    context: Context,
) -> Dijkstra<Context> {
    let mut dijkstra = MulticriteriaDijkstra::from_datastructures(
        datastructures
            .take()
            .expect("Datastructures were not returned"),
        context,
    );
    dijkstra.reset(|_| {});
    dijkstra
}

/// Search all nodes with an alignment cost of at most the given cost limit.
pub fn search_until_cost_limit<Context: AlignmentContext>(
    dijkstra: &mut Dijkstra<Context>,
    cost_limit: AlignmentCostOf<Context>,
) {
    // The match run is not limited, such that only the alignment cost is limited.
    dijkstra.search_until_cost_limit(&NodeCost {
        cost: cost_limit,
        match_run: u32::MAX,
    });
}

/// Search all nodes with an alignment cost of at most the cost of the best target.
///
/// Returns the cost of the best target and its alignment,
/// or `Cost::max_value()` and an empty alignment if no target was found.
///
/// Some nodes exceeding the cost of the best target may be inserted before the best target is found,
/// so the cost of nodes needs to be checked when iterating them after the search.
pub fn search_best_target<Context: AlignmentContext>(
    dijkstra: &mut Dijkstra<Context>,
) -> (AlignmentCostOf<Context>, Alignment) {
    // Nodes with a cost exceeding the best target cost found so far are not inserted.
    let mut best_target_cost = AlignmentCostOf::<Context>::max_value();
    dijkstra.search_with_prune_condition(|context, node| {
        if context.is_target(node) {
            best_target_cost = best_target_cost.min(node.alignment_cost());
        }
        node.alignment_cost() > best_target_cost
    });

    let prune_list = dijkstra.prune_list();
    let Some(best_target) = dijkstra
        .context()
        .target_identifiers()
        .filter_map(|identifier| prune_list.min_first_node_at(&identifier))
        .min_by_key(|node| node.alignment_cost())
    else {
        return (AlignmentCostOf::<Context>::max_value(), Vec::new().into());
    };

    // Backtrack the alignment.
    let mut alignment = Vec::new();
    let mut current = prune_list
        .iter_nodes_with_predecessor_at_identifier(best_target.identifier())
        .find(|node| std::ptr::eq(node.node, best_target))
        .unwrap();
    loop {
        let NodeWithOptionalPredecessor {
            node,
            predecessor_index,
        } = current;
        alignment.extend(node.predecessor_alignment_type());
        let Some(predecessor_index) = predecessor_index else {
            break;
        };
        current = prune_list.get(predecessor_index).unwrap();
    }
    alignment.reverse();

    (best_target.alignment_cost(), alignment.into())
}

/// Iterate over the node with minimum alignment cost of each identifier,
/// if that cost is at most the cost limit.
pub fn iter_min_cost_nodes<Node: AlignmentNode>(
    prune_list: &PruneList<Node>,
    cost_limit: Node::AlignmentCost,
) -> impl Iterator<Item = &Node> {
    prune_list
        .iter_fronts()
        // The front is ordered by alignment cost, so the first node has minimum cost.
        .filter_map(|(_, mut front)| front.next())
        .filter(move |node| node.alignment_cost() <= cost_limit)
}

impl<Cost: Ord> MulticriteriaDijkstraCost for NodeCost<Cost> {
    fn dominates(&self, other: &Self) -> bool {
        self.cost <= other.cost && self.match_run <= other.match_run
    }
}

impl<Cost: Ord + Copy> TwoCriteriaCost for NodeCost<Cost> {
    type First = Cost;
    type Second = u32;

    fn first(&self) -> Self::First {
        self.cost
    }

    fn second(&self) -> Self::Second {
        self.match_run
    }
}

impl<Cost: Display> Display for NodeCost<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.cost, self.match_run)
    }
}

/// Wraps the successor output of [`MulticriteriaDijkstraContext::generate_successors`],
/// attaching the given predecessor index to each node.
///
/// Targets are not reported to the prune list, see the module documentation.
pub struct SuccessorOutput<'output, Output, Index> {
    pub output: &'output mut Output,
    pub predecessor_index: Index,
}

impl<Node, Index: Copy, Output: Extend<NodeWithOptionalPredecessorAndIsTarget<Node, Index>>>
    Extend<Node> for SuccessorOutput<'_, Output, Index>
{
    fn extend<T: IntoIterator<Item = Node>>(&mut self, iter: T) {
        let predecessor_index = self.predecessor_index;
        self.output.extend(
            iter.into_iter()
                .map(|node| NodeWithOptionalPredecessorAndIsTarget {
                    node,
                    predecessor_index: Some(predecessor_index),
                    is_target: false,
                }),
        );
    }
}
