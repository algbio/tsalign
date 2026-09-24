//! The algorithm for computing lower bounds for the case where the anchors are exact matches.

use std::fmt::Display;

use generic_a_star::{AStarContext, AStarIdentifier, AStarNode, cost::AStarCost, reset::Reset};

use crate::{
    alignment::{GapType, coordinates::AlignmentCoordinates},
    costs::GapAffineCosts,
};

pub struct Context<'a, Cost> {
    costs: &'a GapAffineCosts<Cost>,
    max_match_run: u32,
    max_n: usize,
    allow_start_match: bool,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    pub match_run: u32,
    gap_type: GapType,
}

impl<'a, Cost> Context<'a, Cost> {
    pub fn new(
        costs: &'a GapAffineCosts<Cost>,
        max_match_run: u32,
        max_n: usize,
        allow_start_match: bool,
    ) -> Self {
        Self {
            costs,
            max_match_run,
            max_n,
            allow_start_match,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: AlignmentCoordinates::new_primary(0, 0),
                // If we want to enfore the first alignment to be a non-match, we set the match run to u32::MAX.
                match_run: if self.allow_start_match { 0 } else { u32::MAX },
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
                    match_run,
                    gap_type,
                },
            cost,
        } = node;
        let end = AlignmentCoordinates::new_primary(self.max_n, self.max_n);

        if coordinates.can_increment_both(end, None) {
            if *match_run < self.max_match_run {
                // Match
                let new_cost = *cost;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        match_run: match_run + 1,
                        gap_type: GapType::None,
                    },
                    cost: new_cost,
                }));
            }

            // Substitution
            let new_cost = *cost + self.costs.substitution;
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_both(),
                    match_run: 0,
                    gap_type: GapType::None,
                },
                cost: new_cost,
            }));
        }

        if coordinates.can_increment_a_or_ancestor(end, None) {
            // Gap in b
            let new_cost = *cost
                + match gap_type {
                    GapType::InB => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_a(),
                    match_run: 0,
                    gap_type: GapType::InB,
                },
                cost: new_cost,
            }));
        }

        if coordinates.can_increment_b_or_descendant(end, None) {
            // Gap in a
            let new_cost = *cost
                + match gap_type {
                    GapType::InA => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_b(),
                    match_run: 0,
                    gap_type: GapType::InA,
                },
                cost: new_cost,
            }));
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
            self.match_run,
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
