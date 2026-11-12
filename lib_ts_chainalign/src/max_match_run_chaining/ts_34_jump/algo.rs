use std::fmt::Display;

use generic_a_star::cost::AStarCost;
use multicriteria_dijkstra::{
    MulticriteriaDijkstraContext,
    node::{MulticriteriaDijkstraIdentifier, MulticriteriaDijkstraNode},
    prune_list::{IndexedNode, NodeWithIsTarget, NodeWithOptionalPredecessorAndIsTarget},
};

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
    costs::AlignmentCosts,
    max_match_run_chaining::alignment_search::{
        AlignmentContext, AlignmentNode, NodeCost, SuccessorOutput,
    },
};

const DEBUG_EXACT_34_JUMP: bool = false;

pub struct Context<'costs, 'sequences, 'rc_fn, Cost> {
    costs: &'costs AlignmentCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    start: SpecificSecondaryAlignmentCoordinates,
    end: PrimaryAlignmentCoordinates,
    enforce_non_match: bool,
    max_match_run: u32,
    enforcement_offset: u8,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    /// The alignment operation from the predecessor to this node.
    ///
    /// This is `None` for the root node and for nodes entering the final unenforced section.
    pub predecessor_alignment_type: Option<AlignmentType>,
    pub cost: NodeCost<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Identifier {
    Secondary {
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    },
    Jump34 {
        coordinates: PrimaryAlignmentCoordinates,
        has_non_match: bool,
    },
    Primary {
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        final_unenforced_section: bool,
    },
}

impl<'costs, 'sequences, 'rc_fn, Cost> Context<'costs, 'sequences, 'rc_fn, Cost> {
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        costs: &'costs AlignmentCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        start: SpecificSecondaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
        enforce_non_match: bool,
        max_match_run: u32,
        enforcement_offset: u8,
    ) -> Self {
        Self {
            costs,
            sequences,
            rc_fn,
            start,
            end,
            enforce_non_match,
            max_match_run,
            enforcement_offset,
        }
    }

    /// Returns true if nodes with the given identifier are targets.
    pub fn is_target_identifier(&self, identifier: &Identifier) -> bool {
        identifier.coordinates(self.ts_kind()) == self.end.into()
            && (identifier.has_non_match() || !self.enforce_non_match)
    }

    pub fn ts_kind(&self) -> TsKind {
        self.start.ts_kind()
    }
}

impl<Cost: AStarCost> MulticriteriaDijkstraContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn generate_root_nodes(&self, output: &mut impl Extend<NodeWithIsTarget<Self::Node>>) {
        output.extend(std::iter::once(NodeWithIsTarget {
            node: Node {
                identifier: Identifier::Secondary {
                    coordinates: self.start.into(),
                    gap_type: GapType::None,
                    has_non_match: false,
                },
                predecessor_alignment_type: None,
                cost: NodeCost {
                    cost: Cost::zero(),
                    match_run: 0,
                },
            },
            // Targets are not used for pruning, see `alignment_search`.
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
        let node = &node;
        let Node {
            identifier,
            cost:
                NodeCost {
                    cost,
                    match_run: current_match_run,
                },
            ..
        } = node;
        let match_run = *current_match_run;

        let is_secondary = matches!(identifier, Identifier::Secondary { .. });
        let gap_affine_costs = if is_secondary {
            &self.costs.secondary_costs
        } else {
            &self.costs.primary_costs
        };

        match *identifier {
            Identifier::Secondary {
                coordinates,
                gap_type,
                has_non_match,
            } => {
                let coordinates = coordinates.into_specific(self.ts_kind());

                // Generate gap-affine successors.
                if coordinates.can_increment_both_primary(self.end, Some(self.sequences)) {
                    let (ca, cb) = self.sequences.secondary_characters(coordinates, self.rc_fn);
                    let is_match = ca == cb;

                    if is_match {
                        // Disallow runs of matches longer than the maximum.
                        // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                        if match_run < self.max_match_run
                            || (self.start.ancestor() - coordinates.ancestor()
                                < self.enforcement_offset.into()
                                && coordinates.descendant() - self.start.descendant()
                                    < self.enforcement_offset.into())
                        {
                            // Match
                            let new_cost = *cost;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_secondary(
                                    coordinates.increment_both(1).into(),
                                    GapType::None,
                                    has_non_match,
                                ),
                                predecessor_alignment_type: Some(AlignmentType::Match),
                                cost: NodeCost {
                                    cost: new_cost,
                                    match_run: match_run + 1,
                                },
                            }));
                        }
                    } else {
                        // Substitution
                        let new_cost = *cost + gap_affine_costs.substitution;
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_secondary(
                                coordinates.increment_both(1).into(),
                                GapType::None,
                                true,
                            ),
                            predecessor_alignment_type: Some(AlignmentType::Substitution),
                            cost: NodeCost {
                                cost: new_cost,
                                match_run: 0,
                            },
                        }));
                    }
                }

                if coordinates.can_increment_ancestor_primary(Some(self.sequences)) {
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
                            true,
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }

                if coordinates.can_increment_descendant_primary(self.end, Some(self.sequences)) {
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
                            true,
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapA),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }

                // Generate jump successors.
                let new_cost = *cost;

                // This generates too many jumps, most of these are gonna be much too far.
                output.extend(coordinates.generate_34_jumps(self.end).map(
                    |(jump, coordinates)| {
                        if DEBUG_EXACT_34_JUMP {
                            println!(
                                "Jump from {} to {coordinates} by {jump}",
                                node.identifier.coordinates(self.ts_kind())
                            );
                        }
                        Node {
                            identifier: Identifier::Jump34 {
                                coordinates,
                                has_non_match,
                            },
                            predecessor_alignment_type: Some(AlignmentType::TsEnd { jump }),
                            cost: NodeCost {
                                cost: new_cost,
                                match_run: 0,
                            },
                        }
                    },
                ));
            }

            Identifier::Jump34 {
                coordinates,
                has_non_match,
            }
            | Identifier::Primary {
                coordinates,
                has_non_match,
                ..
            } => {
                let gap_type = identifier.gap_type();

                // Generate a successor into a possible final unenforced section.
                if let Identifier::Primary {
                    coordinates,
                    gap_type,
                    has_non_match,
                    final_unenforced_section: false,
                } = *identifier
                {
                    let identifier = Identifier::Primary {
                        coordinates,
                        gap_type,
                        has_non_match,
                        final_unenforced_section: true,
                    };
                    output.extend(std::iter::once(Node {
                        identifier,
                        predecessor_alignment_type: None,
                        cost: node.cost,
                    }));
                }

                let final_unenforced_section = identifier.final_unenforced_section();

                // Generate gap-affine successors.
                if coordinates.can_increment_both_primary(self.end) {
                    let (ca, cb) = self.sequences.primary_characters(coordinates);
                    let is_match = ca == cb;

                    if is_match {
                        // Disallow runs of matches longer than the maximum.
                        // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                        if match_run < self.max_match_run
                            || (final_unenforced_section
                                && match_run < self.enforcement_offset.into())
                        {
                            // Match
                            let new_cost = *cost;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_primary(
                                    coordinates.increment_both(1),
                                    GapType::None,
                                    has_non_match,
                                    final_unenforced_section,
                                ),
                                predecessor_alignment_type: Some(AlignmentType::Match),
                                cost: NodeCost {
                                    cost: new_cost,
                                    match_run: match_run + 1,
                                },
                            }));
                        }
                    } else if !final_unenforced_section {
                        // Substitution
                        let new_cost = *cost + gap_affine_costs.substitution;
                        output.extend(std::iter::once(Node {
                            identifier: Identifier::new_primary(
                                coordinates.increment_both(1),
                                GapType::None,
                                true,
                                final_unenforced_section,
                            ),
                            predecessor_alignment_type: Some(AlignmentType::Substitution),
                            cost: NodeCost {
                                cost: new_cost,
                                match_run: 0,
                            },
                        }));
                    }
                }

                if coordinates.can_increment_a_primary(self.end) && !final_unenforced_section {
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
                            true,
                            final_unenforced_section,
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }

                if coordinates.can_increment_b_primary(self.end) && !final_unenforced_section {
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
                            true,
                            final_unenforced_section,
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapA),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }
            }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        self.is_target_identifier(&node.identifier)
    }
}

impl<Cost: AStarCost> AlignmentContext for Context<'_, '_, '_, Cost> {
    fn target_identifiers(&self) -> impl Iterator<Item = Identifier> {
        let coordinates = self.end;
        [false, true]
            .into_iter()
            .flat_map(move |has_non_match| {
                std::iter::once(Identifier::Jump34 {
                    coordinates,
                    has_non_match,
                })
                .chain(
                    [GapType::None, GapType::InA, GapType::InB]
                        .into_iter()
                        .flat_map(move |gap_type| {
                            [false, true].map(|final_unenforced_section| Identifier::Primary {
                                coordinates,
                                gap_type,
                                has_non_match,
                                final_unenforced_section,
                            })
                        }),
                )
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

impl Identifier {
    pub fn new_primary(
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        final_unenforced_section: bool,
    ) -> Self {
        Identifier::Primary {
            coordinates,
            gap_type,
            has_non_match,
            final_unenforced_section,
        }
    }

    pub fn new_secondary(
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    ) -> Self {
        Identifier::Secondary {
            coordinates,
            gap_type,
            has_non_match,
        }
    }

    pub fn coordinates(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        match self {
            Identifier::Primary { coordinates, .. } => coordinates.into(),
            Identifier::Jump34 { coordinates, .. } => coordinates.into(),
            Identifier::Secondary { coordinates, .. } => coordinates.into_specific(ts_kind).into(),
        }
    }

    pub fn gap_type(&self) -> GapType {
        match self {
            Identifier::Primary { gap_type, .. } => *gap_type,
            Identifier::Jump34 { .. } => GapType::None,
            Identifier::Secondary { gap_type, .. } => *gap_type,
        }
    }

    pub fn has_non_match(&self) -> bool {
        match self {
            Identifier::Primary { has_non_match, .. } => *has_non_match,
            Identifier::Jump34 { has_non_match, .. } => *has_non_match,
            Identifier::Secondary { has_non_match, .. } => *has_non_match,
        }
    }

    pub fn final_unenforced_section(&self) -> bool {
        match self {
            Identifier::Primary {
                final_unenforced_section,
                ..
            } => *final_unenforced_section,
            _ => false,
        }
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Identifier::Primary {
                coordinates,
                gap_type,
                ..
            } => write!(f, "P({coordinates}, {gap_type})"),
            Identifier::Jump34 { coordinates, .. } => write!(f, "J({coordinates})"),
            Identifier::Secondary {
                coordinates,
                gap_type,
                ..
            } => write!(f, "S({coordinates}, {gap_type})"),
        }
    }
}

impl MulticriteriaDijkstraIdentifier for Identifier {}
