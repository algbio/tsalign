//! The algorithm for computing lower bounds for the case where the anchors are inexact matches.

use std::fmt::Display;

use generic_a_star::{AStarContext, AStarIdentifier, AStarNode, cost::AStarCost, reset::Reset};

use crate::{
    alignment::{GapType, coordinates::AlignmentCoordinates},
    alignment_history::{
        history_alignment_operations::AlignmentHistoryOperation,
        profile::{HistoryProfile, HistoryProfileParameters},
    },
    costs::GapAffineCosts,
};

pub struct Context<'a, Cost> {
    costs: &'a GapAffineCosts<Cost>,
    history_parameters: HistoryProfileParameters,
    max_n: usize,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    pub history: HistoryProfile,
    gap_type: GapType,
}

impl<'a, Cost> Context<'a, Cost> {
    pub fn new(
        costs: &'a GapAffineCosts<Cost>,
        anchor_k: u8,
        max_anchor_mutations: u8,
        max_n: usize,
    ) -> Self {
        Self {
            costs,
            history_parameters: HistoryProfileParameters::new(anchor_k, max_anchor_mutations),
            max_n,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: AlignmentCoordinates::new_primary(0, 0),
                history: HistoryProfile::empty(),
                gap_type: GapType::None,
            },
            cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
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
            if let Some(history) =
                history.try_extend(AlignmentHistoryOperation::Match, self.history_parameters)
            {
                // Match
                let new_cost = *cost;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        history,
                        gap_type: GapType::None,
                    },
                    cost: new_cost,
                }));
            }

            if let Some(history) = history.try_extend(
                AlignmentHistoryOperation::Substitution,
                self.history_parameters,
            ) {
                // Substitution
                let new_cost = *cost + self.costs.substitution;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        history,
                        gap_type: GapType::None,
                    },
                    cost: new_cost,
                }));
            }
        }

        if coordinates.can_increment_a_or_ancestor(end, None) {
            if let Some(history) =
                history.try_extend(AlignmentHistoryOperation::GapInB, self.history_parameters)
            {
                // Gap in b
                let new_cost = *cost
                    + match gap_type {
                        GapType::InB => self.costs.gap_extend,
                        _ => self.costs.gap_open,
                    };
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_a(),
                        history,
                        gap_type: GapType::InB,
                    },
                    cost: new_cost,
                }));
            }
        }

        if coordinates.can_increment_b_or_descendant(end, None) {
            if let Some(history) =
                history.try_extend(AlignmentHistoryOperation::GapInA, self.history_parameters)
            {
                // Gap in a
                let new_cost = *cost
                    + match gap_type {
                        GapType::InA => self.costs.gap_extend,
                        _ => self.costs.gap_open,
                    };
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_b(),
                        history,
                        gap_type: GapType::InA,
                    },
                    cost: new_cost,
                }));
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

impl<Cost> Reset for Context<'_, Cost> {
    fn reset(&mut self) {
        unimplemented!()
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
