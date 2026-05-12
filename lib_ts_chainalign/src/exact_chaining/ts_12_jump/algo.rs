use std::fmt::Display;

use generic_a_star::{
    AStarContext, AStarIdentifier, AStarNode,
    cost::{AStarCost, OrderedPairCost, U32Cost},
    reset::Reset,
};
use num_traits::Zero;

use crate::{
    alignment::{
        AlignmentType, GapType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
    },
    costs::AlignmentCosts,
};

pub struct Context<'costs, 'sequences, 'rc_fn, Cost> {
    costs: &'costs AlignmentCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    enforce_non_match: bool,
    max_match_run: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub predecessor: Option<Identifier>,
    pub predecessor_alignment_type: Option<AlignmentType>,
    pub cost: Cost,
    pub match_run: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Identifier {
    Primary {
        coordinates: AlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    },
    Jump12 {
        coordinates: AlignmentCoordinates,
        has_non_match: bool,
    },
    Secondary {
        coordinates: AlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    },
}

impl<'costs, 'sequences, 'rc_fn, Cost> Context<'costs, 'sequences, 'rc_fn, Cost> {
    pub fn new(
        costs: &'costs AlignmentCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        enforce_non_match: bool,
        max_match_run: u32,
    ) -> Self {
        assert!(start.is_primary());
        assert!(end.is_secondary());

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
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier::Primary {
                coordinates: self.start,
                gap_type: GapType::None,
                has_non_match: false,
            },
            predecessor: None,
            predecessor_alignment_type: None,
            cost: self.costs.ts_base_cost.get(self.end.ts_kind().unwrap()),
            match_run: 0,
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let Node {
            identifier,
            cost,
            match_run,
            ..
        } = node;
        let predecessor = Some(*identifier);

        let coordinates = identifier.coordinates();
        let has_non_match = identifier.has_non_match();
        let gap_type = identifier.gap_type();
        let is_primary = matches!(identifier, Identifier::Primary { .. });
        let gap_affine_costs = if is_primary {
            &self.costs.primary_costs
        } else {
            &self.costs.secondary_costs
        };

        // Generate gap-affine successors.
        if coordinates.can_increment_both(self.end, Some(self.sequences)) {
            let (ca, cb) = self.sequences.characters(coordinates, self.rc_fn);
            let is_match = ca == cb;

            if is_match {
                // Disallow runs of matches longer than the maximum.
                // This is because we do not want the exact chaining to find new anchors (which actually already exist).
                if *match_run < self.max_match_run {
                    // Match
                    let new_cost = *cost;
                    output.extend(std::iter::once(Node {
                        identifier: Identifier::new_primary_secondary(
                            is_primary,
                            coordinates.increment_both(),
                            GapType::None,
                            has_non_match,
                        ),
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::Match),
                        cost: new_cost,
                        match_run: match_run + 1,
                    }));
                }
            } else {
                // Substitution
                let new_cost = *cost + gap_affine_costs.substitution;
                output.extend(std::iter::once(Node {
                    identifier: Identifier::new_primary_secondary(
                        is_primary,
                        coordinates.increment_both(),
                        GapType::None,
                        true,
                    ),
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::Substitution),
                    cost: new_cost,
                    match_run: 0,
                }));
            }
        }

        if coordinates.can_increment_a(self.end, Some(self.sequences)) {
            // Gap in b
            let new_cost = *cost
                + match gap_type {
                    GapType::InB => gap_affine_costs.gap_extend,
                    _ => gap_affine_costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier::new_primary_secondary(
                    is_primary,
                    coordinates.increment_a(),
                    GapType::InB,
                    true,
                ),
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::GapB),
                cost: new_cost,
                match_run: 0,
            }));
        }

        if coordinates.can_increment_b(self.end, Some(self.sequences)) {
            // Gap in a
            let new_cost = *cost
                + match gap_type {
                    GapType::InA => gap_affine_costs.gap_extend,
                    _ => gap_affine_costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier::new_primary_secondary(
                    is_primary,
                    coordinates.increment_b(),
                    GapType::InA,
                    true,
                ),
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::GapA),
                cost: new_cost,
                match_run: 0,
            }));
        }

        // Generate jump successors.
        if is_primary {
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
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::TsStart {
                            jump,
                            ts_kind: coordinates.ts_kind().unwrap(),
                        }),
                        cost: new_cost,
                        match_run: 0,
                    }),
            );
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.coordinates() == self.end
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

    // Use match run as secondary cost
    type Cost = OrderedPairCost<Cost, U32Cost>;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        OrderedPairCost(self.cost, U32Cost::from_primitive(self.match_run))
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        OrderedPairCost(Cost::zero(), U32Cost::zero())
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
    pub fn new_primary_secondary(
        is_primary: bool,
        coordinates: AlignmentCoordinates,
        gap_type: GapType,
        has_non_match: bool,
    ) -> Self {
        if is_primary {
            Identifier::Primary {
                coordinates,
                gap_type,
                has_non_match,
            }
        } else {
            Identifier::Secondary {
                coordinates,
                gap_type,
                has_non_match,
            }
        }
    }

    pub fn coordinates(&self) -> AlignmentCoordinates {
        match self {
            Identifier::Primary { coordinates, .. } => *coordinates,
            Identifier::Jump12 { coordinates, .. } => *coordinates,
            Identifier::Secondary { coordinates, .. } => *coordinates,
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
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}: {}, {}",
            self.identifier,
            if let Some(predecessor) = &self.predecessor {
                format!("<-{predecessor}")
            } else {
                "".to_string()
            },
            self.cost,
            self.match_run
        )
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}({}, {})",
            match self {
                Identifier::Primary { .. } => "P".to_string(),
                Identifier::Jump12 { .. } => "J".to_string(),
                Identifier::Secondary { .. } => "S".to_string(),
            },
            self.coordinates(),
            self.gap_type(),
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
        self.cost
            .cmp(&other.cost)
            .then_with(|| self.match_run.cmp(&other.match_run))
    }
}

impl AStarIdentifier for Identifier {}
