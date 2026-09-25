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
    costs::GapAffineCosts,
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
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    gap_type: GapType,
    pub has_non_match: bool,
}

impl<'costs, 'sequences, 'rc_fn, Cost> Context<'costs, 'sequences, 'rc_fn, Cost> {
    ///
    /// **Parameters:**
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
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: self.start,
                gap_type: GapType::None,
                has_non_match: false,
            },
            predecessor: None,
            predecessor_alignment_type: None,
            cost: Cost::zero(),
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
        let Identifier {
            coordinates,
            gap_type,
            has_non_match,
        } = *identifier;

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
                        identifier: Identifier {
                            coordinates: coordinates.increment_both(),
                            gap_type: GapType::None,
                            has_non_match,
                        },
                        predecessor,
                        predecessor_alignment_type: Some(AlignmentType::Match),
                        cost: new_cost,
                        match_run: match_run + 1,
                    }));
                }
            } else {
                // Substitution
                let new_cost = *cost + self.costs.substitution;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        gap_type: GapType::None,
                        has_non_match: true,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::Substitution),
                    cost: new_cost,
                    match_run: 0,
                }));
            }
        }

        if coordinates.can_increment_a_or_ancestor(self.end, Some(self.sequences)) {
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
                    has_non_match: true,
                },
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::GapB),
                cost: new_cost,
                match_run: 0,
            }));
        }

        if coordinates.can_increment_b_or_descendant(self.end, Some(self.sequences)) {
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
                    has_non_match: true,
                },
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::GapA),
                cost: new_cost,
                match_run: 0,
            }));
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        if node.identifier.coordinates != self.end {
            false
        } else if self.start == self.end {
            self.allow_direct_chaining
        } else {
            node.identifier.has_non_match || self.allow_all_matches
        }
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

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}, {}", self.identifier, self.cost, self.match_run)
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
        self.cost
            .cmp(&other.cost)
            .then_with(|| self.match_run.cmp(&other.match_run))
    }
}

impl AStarIdentifier for Identifier {}
