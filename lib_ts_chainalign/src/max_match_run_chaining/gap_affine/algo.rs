use std::fmt::Display;

use generic_a_star::cost::AStarCost;
use multicriteria_dijkstra::{
    MulticriteriaDijkstraContext,
    node::{MulticriteriaDijkstraIdentifier, MulticriteriaDijkstraNode},
    prune_list::{IndexedNode, NodeWithIsTarget, NodeWithOptionalPredecessorAndIsTarget},
};

use crate::{
    alignment::{
        AlignmentType, GapType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
    },
    costs::GapAffineCosts,
    max_match_run_chaining::alignment_search::{AlignmentContext, AlignmentNode, NodeCost},
};

pub struct Context<'costs, 'sequences, 'rc_fn, Cost> {
    costs: &'costs GapAffineCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    allow_direct_chaining: bool,
    allow_all_matches: bool,
    max_match_run: u32,
    enforcement_offset: u8,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub predecessor_alignment_type: Option<AlignmentType>,
    pub cost: NodeCost<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    pub gap_type: GapType,
    pub has_non_match: bool,
}

impl<'costs, 'sequences, 'rc_fn, Cost> Context<'costs, 'sequences, 'rc_fn, Cost> {
    /// # Parameters:
    ///
    /// * `allow_direct_chaining`: allow reaching the target with zero alignment operations.
    /// * `allow_all_matches`: allow reaching the target before any non-match was used. This still obeys the `max_match_run` parameter.
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        costs: &'costs GapAffineCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        allow_direct_chaining: bool,
        allow_all_matches: bool,
        max_match_run: u32,
        enforcement_offset: u8,
    ) -> Self {
        Self {
            costs,
            sequences,
            rc_fn,
            start,
            end,
            allow_direct_chaining,
            allow_all_matches,
            max_match_run,
            enforcement_offset,
        }
    }

    /// Returns true if nodes with the given identifier are targets.
    pub fn is_target_identifier(&self, identifier: &Identifier) -> bool {
        if identifier.coordinates != self.end {
            false
        } else if self.start == self.end {
            self.allow_direct_chaining
        } else {
            identifier.has_non_match || self.allow_all_matches
        }
    }
}

impl<Cost: AStarCost> MulticriteriaDijkstraContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn generate_root_nodes(&self, output: &mut impl Extend<NodeWithIsTarget<Self::Node>>) {
        output.extend(std::iter::once(NodeWithIsTarget {
            node: Node {
                identifier: Identifier {
                    coordinates: self.start,
                    gap_type: GapType::None,
                    has_non_match: false,
                },
                predecessor_alignment_type: None,
                cost: NodeCost {
                    cost: Cost::zero(),
                    match_run: 0,
                },
            },
            // Targets are not used for pruning, see `GapAffineAligner`.
            is_target: false,
        }));
    }

    fn generate_successors<Index: Copy>(
        &self,
        IndexedNode { node, index }: IndexedNode<Self::Node, Index>,
        output: &mut impl Extend<NodeWithOptionalPredecessorAndIsTarget<Self::Node, Index>>,
    ) {
        let Node {
            identifier:
                Identifier {
                    coordinates,
                    gap_type,
                    has_non_match,
                },
            cost: NodeCost { cost, match_run },
            ..
        } = node;
        let mut emit = |identifier, alignment_type, cost, match_run| {
            output.extend(std::iter::once(NodeWithOptionalPredecessorAndIsTarget {
                node: Node {
                    identifier,
                    predecessor_alignment_type: Some(alignment_type),
                    cost: NodeCost { cost, match_run },
                },
                predecessor_index: Some(index),
                // Targets are not used for pruning, see `GapAffineAligner`.
                is_target: false,
            }));
        };

        let enforce_max_match_run = if let Some(primary) = coordinates.into_primary() {
            let start = self.start.into_primary().unwrap();
            primary.a() - start.a() >= self.enforcement_offset.into()
                || primary.b() - start.b() >= self.enforcement_offset.into()
        } else {
            let secondary = coordinates.into_secondary().unwrap();
            let start = self.start.into_secondary().unwrap();
            start.ancestor() - secondary.ancestor() >= self.enforcement_offset.into()
                || secondary.descendant() - start.descendant() >= self.enforcement_offset.into()
        };

        if coordinates.can_increment_both(self.end, Some(self.sequences)) {
            let (ca, cb) = self.sequences.characters(coordinates, self.rc_fn);
            let is_match = ca == cb;

            if is_match {
                // Disallow runs of matches longer than the maximum.
                // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                if match_run < self.max_match_run || !enforce_max_match_run {
                    // Match
                    emit(
                        Identifier {
                            coordinates: coordinates.increment_both(),
                            gap_type: GapType::None,
                            has_non_match,
                        },
                        AlignmentType::Match,
                        cost,
                        match_run + 1,
                    );
                }
            } else {
                // Substitution
                emit(
                    Identifier {
                        coordinates: coordinates.increment_both(),
                        gap_type: GapType::None,
                        has_non_match: true,
                    },
                    AlignmentType::Substitution,
                    cost + self.costs.substitution,
                    0,
                );
            }
        }

        if coordinates.can_increment_a_or_ancestor(self.end, Some(self.sequences)) {
            // Gap in b
            emit(
                Identifier {
                    coordinates: coordinates.increment_a(),
                    gap_type: GapType::InB,
                    has_non_match: true,
                },
                AlignmentType::GapB,
                cost + match gap_type {
                    GapType::InB => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                },
                0,
            );
        }

        if coordinates.can_increment_b_or_descendant(self.end, Some(self.sequences)) {
            // Gap in a
            emit(
                Identifier {
                    coordinates: coordinates.increment_b(),
                    gap_type: GapType::InA,
                    has_non_match: true,
                },
                AlignmentType::GapA,
                cost + match gap_type {
                    GapType::InA => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                },
                0,
            );
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        self.is_target_identifier(&node.identifier)
    }
}

impl<Cost: AStarCost> AlignmentContext for Context<'_, '_, '_, Cost> {
    fn target_identifiers(&self) -> impl Iterator<Item = Identifier> {
        [GapType::None, GapType::InA, GapType::InB]
            .into_iter()
            .flat_map(|gap_type| {
                [false, true].map(|has_non_match| Identifier {
                    coordinates: self.end,
                    gap_type,
                    has_non_match,
                })
            })
            .filter(|identifier| self.is_target_identifier(identifier))
    }
}

impl<Cost: AStarCost> AlignmentNode for Node<Cost> {
    type AlignmentCost = Cost;

    fn predecessor_alignment_type(&self) -> Option<AlignmentType> {
        self.predecessor_alignment_type
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
            self.coordinates, self.gap_type, self.has_non_match,
        )
    }
}
