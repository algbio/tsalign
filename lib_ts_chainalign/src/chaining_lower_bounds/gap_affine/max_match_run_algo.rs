//! The algorithm for computing lower bounds for the case where the match runs of the alignment are limited.
//!
//! This uses a bicriteria Dijkstra, minimising the alignment cost and the length of the current match run,
//! like the aligners in [`max_match_run_chaining`](crate::max_match_run_chaining).

use std::fmt::Display;

use generic_a_star::cost::AStarCost;
use multicriteria_dijkstra::{
    MulticriteriaDijkstraContext,
    node::{MulticriteriaDijkstraIdentifier, MulticriteriaDijkstraNode},
    prune_list::{IndexedNode, NodeWithIsTarget, NodeWithOptionalPredecessorAndIsTarget},
};

use crate::{
    alignment::{GapType, coordinates::PrimaryAlignmentCoordinates},
    costs::GapAffineCosts,
    max_match_run_chaining::alignment_search::{NodeCost, SuccessorOutput},
};

pub struct Context<'a, Cost> {
    costs: &'a GapAffineCosts<Cost>,
    max_match_run: u32,
    /// Start enforcing the max match run only when at least one sequence has reached this length.
    enforcement_offset: u8,
    max_n: usize,
    allow_start_match: bool,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub cost: NodeCost<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Identifier {
    pub coordinates: PrimaryAlignmentCoordinates,
    gap_type: GapType,
}

impl<'a, Cost> Context<'a, Cost> {
    pub fn new(
        costs: &'a GapAffineCosts<Cost>,
        max_match_run: u32,
        enforcement_offset: u8,
        max_n: usize,
        allow_start_match: bool,
    ) -> Self {
        // Forbidding a start match sets the match_run of the root to u32::MAX,
        // and a non-zero enforcement offset would allow to extend it with a match, overflowing the match_run.
        assert!(
            allow_start_match || enforcement_offset == 0,
            "Forbidding a start match is not supported with a non-zero enforcement offset."
        );

        Self {
            costs,
            max_match_run,
            enforcement_offset,
            max_n,
            allow_start_match,
        }
    }
}

impl<Cost: AStarCost> MulticriteriaDijkstraContext for Context<'_, Cost> {
    type Node = Node<Cost>;

    fn generate_root_nodes(&self, output: &mut impl Extend<NodeWithIsTarget<Self::Node>>) {
        output.extend(std::iter::once(NodeWithIsTarget {
            node: Node {
                identifier: Identifier {
                    coordinates: PrimaryAlignmentCoordinates::new(0, 0),
                    gap_type: GapType::None,
                },
                cost: NodeCost {
                    cost: Cost::zero(),
                    // If we want to enforce the first alignment operation to be a non-match, we set the match run to u32::MAX.
                    match_run: if self.allow_start_match { 0 } else { u32::MAX },
                },
            },
            // There are no targets, the whole matrix is filled.
            is_target: false,
        }));
    }

    fn generate_successors<Index: Copy>(
        &self,
        IndexedNode { node, index }: IndexedNode<Self::Node, Index>,
        output: &mut impl Extend<NodeWithOptionalPredecessorAndIsTarget<Self::Node, Index>>,
    ) {
        let output = &mut SuccessorOutput {
            output,
            predecessor_index: index,
        };
        let Node {
            identifier:
                Identifier {
                    coordinates,
                    gap_type,
                },
            cost: NodeCost { cost, match_run },
        } = node;
        let end = PrimaryAlignmentCoordinates::new(self.max_n, self.max_n);

        if coordinates.can_increment_both_primary(end) {
            if match_run < self.max_match_run
                || (coordinates.a() < self.enforcement_offset.into()
                    && coordinates.b() < self.enforcement_offset.into())
            {
                // Match
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(1),
                        gap_type: GapType::None,
                    },
                    cost: NodeCost {
                        cost,
                        match_run: match_run + 1,
                    },
                }));
            }

            // Substitution
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_both(1),
                    gap_type: GapType::None,
                },
                cost: NodeCost {
                    cost: cost + self.costs.substitution,
                    match_run: 0,
                },
            }));
        }

        if coordinates.can_increment_a_primary(end) {
            // Gap in b
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_a(),
                    gap_type: GapType::InB,
                },
                cost: NodeCost {
                    cost: cost
                        + match gap_type {
                            GapType::InB => self.costs.gap_extend,
                            _ => self.costs.gap_open,
                        },
                    match_run: 0,
                },
            }));
        }

        if coordinates.can_increment_b_primary(end) {
            // Gap in a
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_b(),
                    gap_type: GapType::InA,
                },
                cost: NodeCost {
                    cost: cost
                        + match gap_type {
                            GapType::InA => self.costs.gap_extend,
                            _ => self.costs.gap_open,
                        },
                    match_run: 0,
                },
            }));
        }
    }

    fn is_target(&self, _node: &Self::Node) -> bool {
        // Run until whole matrix is filled.
        false
    }
}

impl<Cost: AStarCost> MulticriteriaDijkstraNode for Node<Cost> {
    type Identifier = Identifier;

    type Cost = NodeCost<Cost>;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> &Self::Cost {
        &self.cost
    }
}

impl MulticriteriaDijkstraIdentifier for Identifier {}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {})",
            self.coordinates.a(),
            self.coordinates.b(),
            self.gap_type,
        )
    }
}
