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

pub struct Context<'costs, 'sequences, 'rc_fn, Cost> {
    costs: &'costs AlignmentCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    start: PrimaryAlignmentCoordinates,
    end: SpecificSecondaryAlignmentCoordinates,
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
    Primary {
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    },
    Jump12 {
        coordinates: AnySecondaryAlignmentCoordinates,
        has_non_match: bool,
    },
    Secondary {
        coordinates: AnySecondaryAlignmentCoordinates,
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
        start: PrimaryAlignmentCoordinates,
        end: SpecificSecondaryAlignmentCoordinates,
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
        self.end.ts_kind()
    }
}

impl<Cost: AStarCost> MulticriteriaDijkstraContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn generate_root_nodes(&self, output: &mut impl Extend<NodeWithIsTarget<Self::Node>>) {
        output.extend(std::iter::once(NodeWithIsTarget {
            node: Node {
                identifier: Identifier::Primary {
                    coordinates: self.start,
                    gap_type: GapType::None,
                    has_non_match: false,
                },
                predecessor_alignment_type: None,
                cost: NodeCost {
                    cost: self.costs.ts_base_cost.get(self.end.ts_kind()),
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
                has_non_match,
            } => {
                // Generate gap-affine successors.
                if coordinates.can_increment_both_secondary(self.end, Some(self.sequences)) {
                    let (ca, cb) = self.sequences.primary_characters(coordinates);
                    let is_match = ca == cb;

                    if is_match {
                        // Disallow runs of matches longer than the maximum.
                        // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                        if match_run < self.max_match_run
                            || (coordinates.a() - self.start.a() < self.enforcement_offset.into()
                                && coordinates.b() - self.start.b()
                                    < self.enforcement_offset.into())
                        {
                            // Match
                            let new_cost = *cost;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_primary(
                                    coordinates.increment_both(1),
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
                            identifier: Identifier::new_primary(
                                coordinates.increment_both(1),
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

                if coordinates.can_increment_a_secondary(self.end, Some(self.sequences)) {
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
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }

                if coordinates.can_increment_b_secondary(self.end, Some(self.sequences)) {
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
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapA),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
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
                            identifier: Identifier::Jump12 {
                                coordinates,
                                has_non_match,
                            },
                            predecessor_alignment_type: Some(AlignmentType::TsStart {
                                jump,
                                ts_kind: self.ts_kind(),
                            }),
                            cost: NodeCost {
                                cost: new_cost,
                                match_run: 0,
                            },
                        }),
                );
            }
            Identifier::Jump12 {
                coordinates,
                has_non_match,
            }
            | Identifier::Secondary {
                coordinates,
                has_non_match,
                ..
            } => {
                let coordinates = coordinates.into_specific(self.ts_kind());
                let gap_type = identifier.gap_type();

                // Generate a successor into a possible final unenforced section.
                if let Identifier::Secondary {
                    coordinates,
                    gap_type,
                    has_non_match,
                    final_unenforced_section: false,
                } = *identifier
                {
                    let identifier = Identifier::Secondary {
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
                if coordinates.can_increment_both_secondary(self.end) {
                    let (ca, cb) = self.sequences.secondary_characters(coordinates, self.rc_fn);
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
                                identifier: Identifier::new_secondary(
                                    coordinates.increment_both(1).into(),
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
                            identifier: Identifier::new_secondary(
                                coordinates.increment_both(1).into(),
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

                if coordinates.can_increment_ancestor_secondary(self.end)
                    && !final_unenforced_section
                {
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
                            final_unenforced_section,
                        ),
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: NodeCost {
                            cost: new_cost,
                            match_run: 0,
                        },
                    }));
                }

                if coordinates.can_increment_descendant_secondary(self.end)
                    && !final_unenforced_section
                {
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
        let coordinates = self.end.into();
        [false, true]
            .into_iter()
            .flat_map(move |has_non_match| {
                std::iter::once(Identifier::Jump12 {
                    coordinates,
                    has_non_match,
                })
                .chain(
                    [GapType::None, GapType::InA, GapType::InB]
                        .into_iter()
                        .flat_map(move |gap_type| {
                            [false, true].map(|final_unenforced_section| Identifier::Secondary {
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
    ) -> Self {
        Identifier::Primary {
            coordinates,
            gap_type,
            has_non_match,
        }
    }

    pub fn new_secondary(
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        final_unenforced_section: bool,
    ) -> Self {
        Identifier::Secondary {
            coordinates,
            gap_type,
            has_non_match,
            final_unenforced_section,
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

    pub fn has_non_match(&self) -> bool {
        match self {
            Identifier::Primary { has_non_match, .. } => *has_non_match,
            Identifier::Jump12 { has_non_match, .. } => *has_non_match,
            Identifier::Secondary { has_non_match, .. } => *has_non_match,
        }
    }

    pub fn final_unenforced_section(&self) -> bool {
        match self {
            Identifier::Secondary {
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
            Self::Primary {
                coordinates,
                gap_type,
                has_non_match,
            } => write!(f, "P({}, {}, {})", coordinates, gap_type, has_non_match,),
            Self::Jump12 {
                coordinates,
                has_non_match,
            } => write!(f, "J({}, {})", coordinates, has_non_match),
            Self::Secondary {
                coordinates,
                gap_type,
                has_non_match,
                final_unenforced_section,
            } => write!(
                f,
                "S({}, {}, {}, {})",
                coordinates, gap_type, has_non_match, final_unenforced_section,
            ),
        }
    }
}

impl MulticriteriaDijkstraIdentifier for Identifier {}
