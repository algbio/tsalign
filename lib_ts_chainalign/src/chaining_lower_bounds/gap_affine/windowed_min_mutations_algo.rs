//! The algorithm for computing lower bounds for the case where the anchors are inexact matches.
//!
//! This is a multi-criteria search: nodes at the same coordinates are compared by cost and by the freedom their alignment history gives to future alignments.
//! Dominated nodes are not expanded, see [`pareto_front`](crate::alignment_history::pareto_front).

use std::fmt::Display;

use generic_a_star::{
    AStarContext, AStarIdentifier, AStarNode, closed_lists::AStarClosedList, cost::AStarCost,
    reset::Reset,
};

use crate::{
    alignment::{GapType, coordinates::AlignmentCoordinates},
    alignment_history::{
        extension_graph::{HistoryExtensionGraph, HistoryExtensionGraphNodeIndex},
        history_alignment_operations::AlignmentHistoryOperation,
        history_vec::AlignmentHistory,
        pareto_front::{ParetoFront, ParetoLabel},
    },
    costs::GapAffineCosts,
};



pub struct Context<'a, Cost, AlignmentHistoryVec> {
    costs: &'a GapAffineCosts<Cost>,
    anchor_k: u8,
    max_anchor_mutations: u8,
    max_n: usize,
    history_graph: HistoryExtensionGraph<AlignmentHistoryVec>,

    /// The Pareto fronts of expanded nodes, indexed by coordinates.
    ///
    /// Since a node identical to an expanded node is dominated, these replace the closed list,
    /// so the search can be run with a [`LastNodeClosedList`].
    /// If `pareto_pruning` is false, then these are empty and a real closed list is needed.
    pareto_fronts: Vec<ParetoFront<Cost>>,
    pareto_pruning: bool,
}

/// A closed list that only remembers the most recently closed node.
///
/// This is used when the context itself detects nodes that do not need to be expanded.
/// The A* algorithm only needs the closed list to retrieve the target node, so the search must abort on the first target.
pub struct LastNodeClosedList<Node> {
    last_node: Option<Node>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    pub history: HistoryExtensionGraphNodeIndex,
    gap_type: GapType,
}

impl<'a, Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory>
    Context<'a, Cost, AlignmentHistoryVec>
{
    /// Creates a new context.
    ///
    /// If `pareto_pruning` is false, then dominated nodes are expanded as well, which requires a closed list.
    /// This is only useful for testing.
    pub fn new(
        costs: &'a GapAffineCosts<Cost>,
        anchor_k: u8,
        max_anchor_mutations: u8,
        max_n: usize,
        pareto_pruning: bool,
    ) -> Self {
        Self {
            costs,
            anchor_k,
            max_anchor_mutations,
            max_n,
            history_graph: HistoryExtensionGraph::new(),
            pareto_fronts: if pareto_pruning {
                (0..(max_n + 1) * (max_n + 1))
                    .map(|_| ParetoFront::new())
                    .collect()
            } else {
                Vec::new()
            },
            pareto_pruning,
        }
    }

    fn pareto_front_index(&self, coordinates: AlignmentCoordinates) -> usize {
        coordinates.primary_ordinate_a().unwrap() * (self.max_n + 1)
            + coordinates.primary_ordinate_b().unwrap()
    }

    fn pareto_label(&self, node: &Node<Cost>) -> ParetoLabel<Cost> {
        ParetoLabel::new(
            node.cost,
            node.identifier.gap_type,
            node.identifier.history,
            &self.history_graph,
        )
    }

    /// Outputs the given successor, unless it is dominated by an expanded node.
    fn output_successor(&self, successor: Node<Cost>, output: &mut impl Extend<Node<Cost>>) {
        if !self.pareto_pruning
            || !self.pareto_fronts[self.pareto_front_index(successor.identifier.coordinates)]
                .is_dominated(&self.pareto_label(&successor), self.costs)
        {
            output.extend(std::iter::once(successor));
        }
    }
}

impl<Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory> AStarContext
    for Context<'_, Cost, AlignmentHistoryVec>
{
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: AlignmentCoordinates::new_primary(0, 0),
                history: self.history_graph.empty_node_id(),
                gap_type: GapType::None,
            },
            cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        if self.pareto_pruning {
            let label = self.pareto_label(node);
            let index = self.pareto_front_index(node.identifier.coordinates);
            if !self.pareto_fronts[index].insert(label, self.costs) {
                // The node is dominated by an expanded node, so its successors are dominated as well.
                return;
            }
        }

        let Node {
            identifier:
                Identifier {
                    coordinates,
                    history,
                    gap_type,
                },
            cost,
        } = node;
        let end = AlignmentCoordinates::new_primary(self.max_n, self.max_n);

        if coordinates.can_increment_both(end, None) {
            if let Some(history) = self.history_graph.try_extend(
                *history,
                AlignmentHistoryOperation::Match,
                self.anchor_k,
                self.max_anchor_mutations,
            ) {
                // Match
                let new_cost = *cost;
                self.output_successor(
                    Node {
                        identifier: Identifier {
                            coordinates: coordinates.increment_both(),
                            history,
                            gap_type: GapType::None,
                        },
                        cost: new_cost,
                    },
                    output,
                );
            }

            if let Some(history) = self.history_graph.try_extend(
                *history,
                AlignmentHistoryOperation::Substitution,
                self.anchor_k,
                self.max_anchor_mutations,
            ) {
                // Substitution
                let new_cost = *cost + self.costs.substitution;
                self.output_successor(
                    Node {
                        identifier: Identifier {
                            coordinates: coordinates.increment_both(),
                            history,
                            gap_type: GapType::None,
                        },
                        cost: new_cost,
                    },
                    output,
                );
            }
        }

        if coordinates.can_increment_a_or_ancestor(end, None) {
            if let Some(history) = self.history_graph.try_extend(
                *history,
                AlignmentHistoryOperation::GapInB,
                self.anchor_k,
                self.max_anchor_mutations,
            ) {
                // Gap in b
                let new_cost = *cost
                    + match gap_type {
                        GapType::InB => self.costs.gap_extend,
                        _ => self.costs.gap_open,
                    };
                self.output_successor(
                    Node {
                        identifier: Identifier {
                            coordinates: coordinates.increment_a(),
                            history,
                            gap_type: GapType::InB,
                        },
                        cost: new_cost,
                    },
                    output,
                );
            }
        }

        if coordinates.can_increment_b_or_descendant(end, None) {
            if let Some(history) = self.history_graph.try_extend(
                *history,
                AlignmentHistoryOperation::GapInA,
                self.anchor_k,
                self.max_anchor_mutations,
            ) {
                // Gap in a
                let new_cost = *cost
                    + match gap_type {
                        GapType::InA => self.costs.gap_extend,
                        _ => self.costs.gap_open,
                    };
                self.output_successor(
                    Node {
                        identifier: Identifier {
                            coordinates: coordinates.increment_b(),
                            history,
                            gap_type: GapType::InA,
                        },
                        cost: new_cost,
                    },
                    output,
                );
            }
        }
    }

    fn is_target(&self, _node: &Self::Node) -> bool {
        // Run until whole matrix is filled
        false
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<Cost, AlignmentHistoryVec> Reset for Context<'_, Cost, AlignmentHistoryVec> {
    fn reset(&mut self) {
        unimplemented!()
    }
}

impl<Node: AStarNode> AStarClosedList<Node> for LastNodeClosedList<Node> {
    fn new() -> Self {
        Self { last_node: None }
    }

    fn len(&self) -> usize {
        usize::from(self.last_node.is_some())
    }

    fn insert(&mut self, identifier: <Node as AStarNode>::Identifier, node: Node) -> Option<Node> {
        self.last_node
            .replace(node)
            .filter(|previous_node| *previous_node.identifier() == identifier)
    }

    fn get(&self, identifier: &<Node as AStarNode>::Identifier) -> Option<&Node> {
        self.last_node
            .as_ref()
            .filter(|node| node.identifier() == identifier)
    }

    fn iter<'this: 'node, 'node>(
        &'this self,
    ) -> impl use<'this, 'node, Node> + Iterator<Item = &'node Node>
    where
        Node: 'node,
    {
        self.last_node.iter()
    }
}

impl<Node> Reset for LastNodeClosedList<Node> {
    fn reset(&mut self) {
        self.last_node = None;
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = ();

    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        Cost::zero()
    }

    fn secondary_maximisable_score(&self) -> usize {
        0
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        // Backtracking not supported
        None
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        // Backtracking not supported
        None
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {}, {})",
            self.coordinates.primary_ordinate_a().unwrap(),
            self.coordinates.primary_ordinate_b().unwrap(),
            self.history,
            self.gap_type,
        )
    }
}

impl<Cost: Ord> PartialOrd for Node<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Cost: Ord> Ord for Node<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.cmp(&other.cost)
    }
}

impl AStarIdentifier for Identifier {}
