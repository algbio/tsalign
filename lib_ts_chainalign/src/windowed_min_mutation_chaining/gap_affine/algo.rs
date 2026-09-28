use std::fmt::Display;

use generic_a_star::{AStarContext, AStarIdentifier, AStarNode, cost::AStarCost, reset::Reset};

use crate::{
    alignment::{
        AlignmentType, GapType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
    },
    alignment_history::{
        extension_graph::{HistoryExtensionGraph, HistoryExtensionGraphNodeIndex},
        history_alignment_operations::AlignmentHistoryOperation,
        history_vec::AlignmentHistory,
    },
    costs::GapAffineCosts,
};

pub struct Context<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec> {
    costs: &'costs GapAffineCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    history_graph: &'history mut HistoryExtensionGraph<AlignmentHistoryVec>,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    anchor_k: u8,
    max_anchor_mutations: u8,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub predecessor: Option<Identifier>,
    pub predecessor_alignment_type: Option<AlignmentType>,
    pub cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    gap_type: GapType,
    history: HistoryExtensionGraphNodeIndex,
}

impl<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec>
    Context<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec>
{
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        costs: &'costs GapAffineCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        history_graph: &'history mut HistoryExtensionGraph<AlignmentHistoryVec>,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self {
        Self {
            costs,
            sequences,
            rc_fn,
            history_graph,
            start,
            end,
            anchor_k,
            max_anchor_mutations,
        }
    }
}

impl<Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory> AStarContext
    for Context<'_, '_, '_, '_, Cost, AlignmentHistoryVec>
{
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: self.start,
                gap_type: GapType::None,
                history: self.history_graph.empty_node_id(),
            },
            predecessor: None,
            predecessor_alignment_type: None,
            cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let Node {
            identifier, cost, ..
        } = node;
        let predecessor = Some(*identifier);
        let Identifier {
            coordinates,
            gap_type,
            history,
        } = *identifier;

        if coordinates.can_increment_both(self.end, Some(self.sequences)) {
            let (ca, cb) = self.sequences.characters(coordinates, self.rc_fn);
            let is_match = ca == cb;

            if is_match {
                if let Some(history) = self.history_graph.try_extend(
                    history,
                    AlignmentHistoryOperation::Match,
                    self.anchor_k,
                    self.max_anchor_mutations,
                ) {
                    // Match
                    let new_cost = *cost;
                    output.extend(std::iter::once(Node {
                        identifier: Identifier {
                            coordinates: coordinates.increment_both(),
                            gap_type: GapType::None,
                            history,
                        },
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::Match),
                        cost: new_cost,
                    }));
                }
            } else if let Some(history) = self.history_graph.try_extend(
                history,
                AlignmentHistoryOperation::Substitution,
                self.anchor_k,
                self.max_anchor_mutations,
            ) {
                // Substitution
                let new_cost = *cost + self.costs.substitution;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        gap_type: GapType::None,
                        history,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::Substitution),
                    cost: new_cost,
                }));
            }
        }

        if coordinates.can_increment_a_or_ancestor(self.end, Some(self.sequences)) {
            if let Some(history) = self.history_graph.try_extend(
                history,
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
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_a(),
                        gap_type: GapType::InB,
                        history,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::GapB),
                    cost: new_cost,
                }));
            }
        }

        if coordinates.can_increment_b_or_descendant(self.end, Some(self.sequences)) {
            if let Some(history) = self.history_graph.try_extend(
                history,
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
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_b(),
                        gap_type: GapType::InA,
                        history,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::GapA),
                    cost: new_cost,
                }));
            }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.coordinates == self.end
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }

    fn is_label_setting(&self) -> bool {
        !self.costs.has_zero_cost()
    }
}

impl<Cost, AlignmentHistoryVec> Reset for Context<'_, '_, '_, '_, Cost, AlignmentHistoryVec> {
    fn reset(&mut self) {
        unimplemented!()
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = AlignmentType;

    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        Self::Cost::zero()
    }

    fn secondary_maximisable_score(&self) -> usize {
        0
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        self.predecessor_alignment_type
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.coordinates, self.gap_type)
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
