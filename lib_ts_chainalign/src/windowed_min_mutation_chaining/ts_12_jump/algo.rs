use std::fmt::Display;

use generic_a_star::{AStarContext, AStarIdentifier, AStarNode, cost::AStarCost, reset::Reset};

use crate::{
    alignment::{
        AlignmentType, GapType,
        coordinates::{
            AlignmentCoordinates, AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            SpecificSecondaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    alignment_history::{
        extension_graph::{HistoryExtensionGraph, HistoryExtensionGraphNodeIndex},
        history_alignment_operations::AlignmentHistoryOperation,
        history_vec::AlignmentHistory,
    },
    costs::AlignmentCosts,
};

pub struct Context<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec> {
    costs: &'costs AlignmentCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    history_graph: &'history mut HistoryExtensionGraph<AlignmentHistoryVec>,
    start: PrimaryAlignmentCoordinates,
    end: SpecificSecondaryAlignmentCoordinates,
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
pub enum Identifier {
    Primary {
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        history: HistoryExtensionGraphNodeIndex,
    },
    Jump12 {
        coordinates: AnySecondaryAlignmentCoordinates,
    },
    Secondary {
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        history: HistoryExtensionGraphNodeIndex,
    },
}

impl<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec>
    Context<'costs, 'sequences, 'rc_fn, 'history, Cost, AlignmentHistoryVec>
{
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        costs: &'costs AlignmentCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        history_graph: &'history mut HistoryExtensionGraph<AlignmentHistoryVec>,
        start: PrimaryAlignmentCoordinates,
        end: SpecificSecondaryAlignmentCoordinates,
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

    pub fn ts_kind(&self) -> TsKind {
        self.end.ts_kind()
    }
}

impl<Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory> AStarContext
    for Context<'_, '_, '_, '_, Cost, AlignmentHistoryVec>
{
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier::Primary {
                coordinates: self.start,
                gap_type: GapType::None,
                history: self.history_graph.empty_node_id(),
            },
            predecessor: None,
            predecessor_alignment_type: None,
            cost: self.costs.ts_base_cost.get(self.end.ts_kind()),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let Node {
            identifier, cost, ..
        } = node;
        let predecessor = Some(*identifier);

        let is_primary = matches!(identifier, Identifier::Primary { .. });
        let gap_affine_costs = if is_primary {
            &self.costs.primary_costs
        } else {
            &self.costs.secondary_costs
        };

        match *identifier {
            Identifier::Primary {
                coordinates,
                gap_type,
                history,
            } => {
                // Generate gap-affine successors.
                if coordinates.can_increment_both_secondary(self.end, Some(self.sequences)) {
                    let (ca, cb) = self.sequences.primary_characters(coordinates);
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
                                identifier: Identifier::new_primary(
                                    coordinates.increment_both(1),
                                    GapType::None,
                                    history,
                                ),
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
                        let new_cost = *cost + gap_affine_costs.substitution;
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_primary(
                                coordinates.increment_both(1),
                                GapType::None,
                                history,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::Substitution),
                            cost: new_cost,
                        }));
                    }
                }

                if coordinates.can_increment_a_secondary(self.end, Some(self.sequences)) {
                    if let Some(history) = self.history_graph.try_extend(
                        history,
                        AlignmentHistoryOperation::GapInB,
                        self.anchor_k,
                        self.max_anchor_mutations,
                    ) {
                        // Gap in b
                        let new_cost = *cost
                            + match gap_type {
                                GapType::InB => gap_affine_costs.gap_extend,
                                _ => gap_affine_costs.gap_open,
                            };
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_primary(
                                coordinates.increment_a(),
                                GapType::InB,
                                history,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::GapB),
                            cost: new_cost,
                        }));
                    }
                }

                if coordinates.can_increment_b_secondary(self.end, Some(self.sequences)) {
                    if let Some(history) = self.history_graph.try_extend(
                        history,
                        AlignmentHistoryOperation::GapInA,
                        self.anchor_k,
                        self.max_anchor_mutations,
                    ) {
                        // Gap in a
                        let new_cost = *cost
                            + match gap_type {
                                GapType::InA => gap_affine_costs.gap_extend,
                                _ => gap_affine_costs.gap_open,
                            };
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_primary(
                                coordinates.increment_b(),
                                GapType::InA,
                                history,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::GapA),
                            cost: new_cost,
                        }));
                    }
                }

                // Generate jump successors.
                // We do not count the jump costs here, because all paths anyways need to jump at some point.
                // Instead, we count them at the start.
                let new_cost = *cost;

                // This generates too many jumps, most of these are gonna be much too far.
                output.extend(
                    coordinates
                        .generate_12_jumps(self.end, self.sequences.primary_end())
                        .map(|(jump, coordinates)| Node {
                            identifier: Identifier::Jump12 { coordinates },
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::TsStart {
                                jump,
                                ts_kind: self.ts_kind(),
                            }),
                            cost: new_cost,
                        }),
                );
            }
            Identifier::Jump12 { coordinates } | Identifier::Secondary { coordinates, .. } => {
                let coordinates = coordinates.into_specific(self.ts_kind());
                let gap_type = identifier.gap_type();
                let history = identifier
                    .history()
                    .unwrap_or_else(|| self.history_graph.empty_node_id());

                // Generate gap-affine successors.
                if coordinates.can_increment_both_secondary(self.end) {
                    let (ca, cb) = self.sequences.secondary_characters(coordinates, self.rc_fn);
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
                                identifier: Identifier::new_secondary(
                                    coordinates.increment_both(1).into(),
                                    GapType::None,
                                    history,
                                ),
                                predecessor,
                                predecessor_alignment_type: Some(AlignmentType::Match),
                                cost: new_cost,
                            }));
                        }
                    } else {
                        if let Some(history) = self.history_graph.try_extend(
                            history,
                            AlignmentHistoryOperation::Substitution,
                            self.anchor_k,
                            self.max_anchor_mutations,
                        ) {
                            // Substitution
                            let new_cost = *cost + gap_affine_costs.substitution;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_secondary(
                                    coordinates.increment_both(1).into(),
                                    GapType::None,
                                    history,
                                ),
                                predecessor,
                                predecessor_alignment_type: Some(AlignmentType::Substitution),
                                cost: new_cost,
                            }));
                        }
                    }
                }

                if coordinates.can_increment_ancestor_secondary(self.end) {
                    if let Some(history) = self.history_graph.try_extend(
                        history,
                        AlignmentHistoryOperation::GapInB,
                        self.anchor_k,
                        self.max_anchor_mutations,
                    ) {
                        // Gap in b
                        let new_cost = *cost
                            + match gap_type {
                                GapType::InB => gap_affine_costs.gap_extend,
                                _ => gap_affine_costs.gap_open,
                            };
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_secondary(
                                coordinates.increment_ancestor().into(),
                                GapType::InB,
                                history,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::GapB),
                            cost: new_cost,
                        }));
                    }
                }

                if coordinates.can_increment_descendant_secondary(self.end) {
                    if let Some(history) = self.history_graph.try_extend(
                        history,
                        AlignmentHistoryOperation::GapInA,
                        self.anchor_k,
                        self.max_anchor_mutations,
                    ) {
                        // Gap in a
                        let new_cost = *cost
                            + match gap_type {
                                GapType::InA => gap_affine_costs.gap_extend,
                                _ => gap_affine_costs.gap_open,
                            };
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_secondary(
                                coordinates.increment_descendant().into(),
                                GapType::InA,
                                history,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::GapA),
                            cost: new_cost,
                        }));
                    }
                }
            }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.coordinates(self.ts_kind()) == self.end.into()
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }

    fn is_label_setting(&self) -> bool {
        false
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

impl Identifier {
    pub fn new_primary(
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        history: HistoryExtensionGraphNodeIndex,
    ) -> Self {
        Identifier::Primary {
            coordinates,
            gap_type,
            history,
        }
    }

    pub fn new_secondary(
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        history: HistoryExtensionGraphNodeIndex,
    ) -> Self {
        Identifier::Secondary {
            coordinates,
            gap_type,
            history,
        }
    }

    pub fn coordinates(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        match self {
            Identifier::Primary { coordinates, .. } => coordinates.into(),
            Identifier::Jump12 { coordinates, .. } => coordinates.into_specific(ts_kind).into(),
            Identifier::Secondary { coordinates, .. } => coordinates.into_specific(ts_kind).into(),
        }
    }

    pub fn gap_type(&self) -> GapType {
        match self {
            Identifier::Primary { gap_type, .. } => *gap_type,
            Identifier::Jump12 { .. } => GapType::None,
            Identifier::Secondary { gap_type, .. } => *gap_type,
        }
    }

    pub fn history(&self) -> Option<HistoryExtensionGraphNodeIndex> {
        match self {
            Identifier::Primary { history, .. } => Some(*history),
            Identifier::Jump12 { .. } => None,
            Identifier::Secondary { history, .. } => Some(*history),
        }
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}: {}",
            self.identifier,
            if let Some(predecessor) = &self.predecessor {
                format!("<-{predecessor}")
            } else {
                "".to_string()
            },
            self.cost,
        )
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Primary {
                coordinates,
                gap_type,
                ..
            } => write!(f, "P({}, {})", coordinates, gap_type),
            Self::Jump12 { coordinates } => write!(f, "J({})", coordinates),
            Self::Secondary {
                coordinates,
                gap_type,
                ..
            } => write!(f, "S({}, {})", coordinates, gap_type),
        }
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
