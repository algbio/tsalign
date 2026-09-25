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
    costs::AlignmentCosts,
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
    Secondary {
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        match_run: u32,
    },
    Jump34 {
        coordinates: PrimaryAlignmentCoordinates,
        has_non_match: bool,
    },
    Primary {
        coordinates: PrimaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        match_run: u32,
    },
}

impl<'costs, 'sequences, 'rc_fn, Cost> Context<'costs, 'sequences, 'rc_fn, Cost> {
    pub fn new(
        costs: &'costs AlignmentCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        start: SpecificSecondaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
        enforce_non_match: bool,
        max_match_run: u32,
    ) -> Self {
        Self {
            costs,
            sequences,
            rc_fn,
            start,
            end,
            enforce_non_match,
            max_match_run,
        }
    }

    pub fn ts_kind(&self) -> TsKind {
        self.start.ts_kind()
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier::Secondary {
                coordinates: self.start.into(),
                gap_type: GapType::None,
                has_non_match: false,
                match_run: 0,
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
                match_run,
            } => {
                let coordinates = coordinates.into_specific(self.ts_kind());

                // Generate gap-affine successors.
                if coordinates.can_increment_both_primary(self.end, Some(self.sequences)) {
                    let (ca, cb) = self.sequences.secondary_characters(coordinates, self.rc_fn);
                    let is_match = ca == cb;

                    if is_match {
                        // Disallow runs of matches longer than the maximum.
                        // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                        if match_run < self.max_match_run {
                            // Match
                            let new_cost = *cost;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_secondary(
                                    coordinates.increment_both(1).into(),
                                    GapType::None,
                                    has_non_match,
                                    match_run + 1,
                                ),
                                predecessor,
                                predecessor_alignment_type: Some(AlignmentType::Match),
                                cost: new_cost,
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
                                0,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::Substitution),
                            cost: new_cost,
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
                            0,
                        ),
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: new_cost,
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
                            0,
                        ),
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::GapA),
                        cost: new_cost,
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
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::TsEnd { jump }),
                            cost: new_cost,
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
                let match_run = identifier.match_run();

                // Generate gap-affine successors.
                if coordinates.can_increment_both_primary(self.end) {
                    let (ca, cb) = self.sequences.primary_characters(coordinates);
                    let is_match = ca == cb;

                    if is_match {
                        // Disallow runs of matches longer than the maximum.
                        // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                        if match_run < self.max_match_run {
                            // Match
                            let new_cost = *cost;
                            output.extend(std::iter::once(Node {
                                identifier: Identifier::new_primary(
                                    coordinates.increment_both(1),
                                    GapType::None,
                                    has_non_match,
                                    match_run + 1,
                                ),
                                predecessor,
                                predecessor_alignment_type: Some(AlignmentType::Match),
                                cost: new_cost,
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
                                0,
                            ),
                            predecessor,
                            predecessor_alignment_type: Some(AlignmentType::Substitution),
                            cost: new_cost,
                        }));
                    }
                }

                if coordinates.can_increment_a_primary(self.end) {
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
                            0,
                        ),
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::GapB),
                        cost: new_cost,
                    }));
                }

                if coordinates.can_increment_b_primary(self.end) {
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
                            0,
                        ),
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::GapA),
                        cost: new_cost,
                    }));
                }
            }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.coordinates(self.ts_kind()) == self.end.into()
            && (node.identifier.has_non_match() || !self.enforce_non_match)
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

impl<Cost> Reset for Context<'_, '_, '_, Cost> {
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
        has_non_match: bool,
        match_run: u32,
    ) -> Self {
        Identifier::Primary {
            coordinates,
            gap_type,
            has_non_match,
            match_run,
        }
    }

    pub fn new_secondary(
        coordinates: AnySecondaryAlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
        match_run: u32,
    ) -> Self {
        Identifier::Secondary {
            coordinates,
            gap_type,
            has_non_match,
            match_run,
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

    pub fn match_run(&self) -> u32 {
        match self {
            Identifier::Primary { match_run, .. } => *match_run,
            Identifier::Jump34 { .. } => 0,
            Identifier::Secondary { match_run, .. } => *match_run,
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
            Identifier::Primary {
                coordinates,
                gap_type,
                match_run,
                ..
            } => write!(f, "P({coordinates}, {gap_type}, {match_run})"),
            Identifier::Jump34 { coordinates, .. } => write!(f, "J({coordinates})"),
            Identifier::Secondary {
                coordinates,
                gap_type,
                match_run,
                ..
            } => write!(f, "S({coordinates}, {gap_type}, {match_run})"),
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
