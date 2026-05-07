use std::fmt::Display;

use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::{AStarNode, cost::AStarCost};
use num_traits::{Bounded, Zero};
use strategies::{
    AlignmentStrategiesNodeMemory, AlignmentStrategySelector, node_ord::NodeOrdStrategy,
    primary_match::PrimaryMatchStrategy,
    template_switch_min_length::TemplateSwitchMinLengthStrategy,
    template_switch_total_length::TemplateSwitchTotalLengthStrategy,
};

mod alignment_type;
pub mod context;
pub mod display;
mod identifier;
pub mod lower_bounds;
pub mod strategies;

pub use alignment_type::{AlignmentType, equal_cost_range::EqualCostRange};
pub use context::Context;
pub use identifier::{
    GapType, Identifier, TemplateSwitchDirection, TemplateSwitchPrimary, TemplateSwitchSecondary,
};

use crate::config::BaseCost;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node<Strategies: AlignmentStrategySelector> {
    node_data: NodeData<
        <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
            <Strategies as AlignmentStrategySelector>::Cost,
        >>::IdentifierPrimaryExtraData,
        <Strategies as AlignmentStrategySelector>::Cost,
    >,
    strategies: AlignmentStrategiesNodeMemory<Strategies>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct NodeData<PrimaryExtraData: Copy, Cost> {
    identifier: Identifier<PrimaryExtraData>,
    predecessor: Option<Identifier<PrimaryExtraData>>,
    predecessor_edge_type: AlignmentType,
    cost: Cost,
    a_star_lower_bound: Cost,
}

impl<Strategies: AlignmentStrategySelector> AStarNode for Node<Strategies> {
    type Identifier = Identifier<
        <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
            <Strategies as AlignmentStrategySelector>::Cost,
        >>::IdentifierPrimaryExtraData,
    >;

    type EdgeType = AlignmentType;

    type Cost = Strategies::Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.node_data.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.node_data.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        self.node_data.a_star_lower_bound
    }

    fn secondary_maximisable_score(&self) -> usize {
        self.strategies
            .template_switch_total_length
            .template_switch_total_length()
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.node_data.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        Some(self.node_data.predecessor_edge_type)
    }
}

impl<Strategies: AlignmentStrategySelector> Node<Strategies> {
    pub fn new_root_at<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        reference_index: usize,
        query_index: usize,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            node_data: NodeData::create_root(Identifier::Primary {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: 0,
                data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::create_root_identifier_primary_extra_data( context),
            }),
            strategies: AlignmentStrategiesNodeMemory::create_root(context),
        }
    }

    fn generate_primary_diagonal_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Strategies::Cost,
        is_match: bool,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { flank_index, .. }
        | Identifier::PrimaryReentry { flank_index, .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        let alignment_type = match (
            is_match,
            flank_index == successor_flank_index && successor_flank_index == 0,
        ) {
            (true, true) => AlignmentType::PrimaryMatch,
            (true, false) => AlignmentType::PrimaryFlankMatch,
            (false, true) => AlignmentType::PrimarySubstitution,
            (false, false) => AlignmentType::PrimaryFlankSubstitution,
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_diagonal_successor(
                successor_flank_index,
                alignment_type,
                context,
            ),
            cost_increment,
            alignment_type,
            context,
        ))
    }

    fn generate_primary_deletion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { flank_index, .. }
        | Identifier::PrimaryReentry { flank_index, .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        let alignment_type = if flank_index == successor_flank_index && successor_flank_index == 0 {
            AlignmentType::PrimaryDeletion
        } else {
            AlignmentType::PrimaryFlankDeletion
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_deletion_successor(
                successor_flank_index,
                alignment_type,
                context,
            ),
            cost_increment,
            alignment_type,
            context,
        ))
    }

    fn generate_primary_insertion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { flank_index, .. }
        | Identifier::PrimaryReentry { flank_index, .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        let alignment_type = if flank_index == successor_flank_index && successor_flank_index == 0 {
            AlignmentType::PrimaryInsertion
        } else {
            AlignmentType::PrimaryFlankInsertion
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_insertion_successor(
                successor_flank_index,
                alignment_type,
                context,
            ),
            cost_increment,
            alignment_type,
            context,
        ))
    }

    fn generate_initial_template_switch_entrance_successors<
        'result,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &'result self,
        rq_qr_cost_increment: Strategies::Cost,
        rr_qq_cost_increment: Strategies::Cost,
        base_cost: &'result BaseCost<Strategies::Cost>,
        context: &'result Context<SubsequenceType, Strategies>,
    ) -> impl 'result + Iterator<Item = Self> {
        if !matches!(
            self.node_data.identifier,
            Identifier::Primary { .. } | Identifier::PrimaryReentry { .. }
        ) {
            unreachable!("This method is only called on primary nodes.")
        }
        debug_assert!(
            rq_qr_cost_increment < Strategies::Cost::max_value()
                || rr_qq_cost_increment < Strategies::Cost::max_value(),
        );

        self.node_data
            .identifier
            .generate_initial_template_switch_entrance_successors()
            .filter_map(move |identifier| {
                let Identifier::TemplateSwitchEntrance {
                    template_switch_primary,
                    template_switch_secondary,
                    template_switch_direction,
                    template_switch_first_offset,
                    ..
                } = &identifier
                else {
                    unreachable!("This closure is only called on template switch entrances.")
                };

                let base_cost = base_cost.get(
                    *template_switch_primary,
                    *template_switch_secondary,
                    *template_switch_direction,
                );
                let cost_increment = match (*template_switch_primary, *template_switch_secondary) {
                    (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Query)
                    | (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Reference) => {
                        rq_qr_cost_increment
                    }
                    (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Reference)
                    | (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Query) => {
                        rr_qq_cost_increment
                    }
                };

                (base_cost != Strategies::Cost::max_value()
                    && cost_increment != Strategies::Cost::max_value())
                .then(|| {
                    self.generate_successor(
                        identifier,
                        cost_increment + base_cost,
                        AlignmentType::TemplateSwitchEntrance {
                            primary: *template_switch_primary,
                            secondary: *template_switch_secondary,
                            direction: *template_switch_direction,
                            equal_cost_range: EqualCostRange::new_invalid(),
                            first_offset: *template_switch_first_offset,
                        },
                        context,
                    )
                })
            })
    }

    fn generate_template_switch_entrance_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        successor_template_switch_first_offset: isize,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_direction,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_direction,
                template_switch_first_offset: successor_template_switch_first_offset,
            },
            cost_increment,
            AlignmentType::TemplateSwitchEntrance {
                primary: template_switch_primary,
                secondary: template_switch_secondary,
                direction: template_switch_direction,
                equal_cost_range: EqualCostRange::new_invalid(),
                first_offset: successor_template_switch_first_offset,
            },
            context,
        ))
    }

    fn generate_secondary_root_node<
        'this,
        'reference,
        'query,
        'context,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &'this self,
        context: &'context mut Context<'reference, 'query, SubsequenceType, Strategies>,
    ) -> impl use<'this, 'reference, 'query, 'context, SubsequenceType, Strategies>
    + IntoIterator<Item = Self> {
        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_direction,
            template_switch_first_offset,
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        let primary_index = match template_switch_primary {
            TemplateSwitchPrimary::Reference => entrance_reference_index,
            TemplateSwitchPrimary::Query => entrance_query_index,
        };

        let secondary_index = (match template_switch_secondary {
            TemplateSwitchSecondary::Reference => entrance_reference_index,
            TemplateSwitchSecondary::Query => entrance_query_index,
        } as isize
            + template_switch_first_offset) as usize;

        match template_switch_secondary {
            TemplateSwitchSecondary::Reference => {
                debug_assert!(secondary_index <= context.reference.len(), "{self}")
            }
            TemplateSwitchSecondary::Query => {
                debug_assert!(secondary_index <= context.query.len(), "{self}")
            }
        }

        let secondary_root_node = self.generate_successor(
            Identifier::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_direction,
                length: 0,
                primary_index,
                secondary_index,
                gap_type: GapType::None,
            },
            Strategies::Cost::zero(),
            AlignmentType::SecondaryRoot,
            context,
        );

        self.strategies
            .template_switch_min_length_strategy
            .template_switch_min_length_lookahead(secondary_root_node, context)
    }

    fn generate_secondary_diagonal_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        is_match: bool,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_diagonal_successor(),
            cost_increment,
            if is_match {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            },
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_deletion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_deletion_successor(),
            cost_increment,
            AlignmentType::SecondaryDeletion,
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_insertion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_insertion_successor(),
            cost_increment,
            AlignmentType::SecondaryInsertion,
            context,
        ))
    }

    fn generate_initial_template_switch_exit_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let Identifier::Secondary {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_direction,
            primary_index,
            length,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_direction,
                primary_index,
                anti_primary_gap: length.try_into().unwrap(),
            },
            cost_increment,
            AlignmentType::TemplateSwitchExit {
                anti_primary_gap: length.try_into().unwrap(),
            },
            context,
        ))
    }

    fn generate_template_switch_exit_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Strategies::Cost,
        successor_anti_primary_gap: isize,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_direction,
            primary_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_direction,
                primary_index,
                anti_primary_gap: successor_anti_primary_gap,
            },
            cost_increment,
            AlignmentType::TemplateSwitchExit {
                anti_primary_gap: successor_anti_primary_gap,
            },
            context,
        ))
    }

    pub fn generate_primary_reentry_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        context: &Context<SubsequenceType, Strategies>,
        cost_increment: Strategies::Cost,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let identifier @ Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            primary_index,
            anti_primary_gap,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        let (reference_index, query_index) = match template_switch_primary {
            TemplateSwitchPrimary::Reference => {
                let query_index = entrance_query_index as isize + anti_primary_gap;

                // TODO the latter condition should never be true if we generate TS exit nodes correctly.
                if query_index < 0 || query_index as usize >= context.query.len() {
                    return None;
                }

                (primary_index, query_index as usize)
            }
            TemplateSwitchPrimary::Query => {
                let reference_index = entrance_reference_index as isize + anti_primary_gap;

                // TODO the latter condition should never be true if we generate TS exit nodes correctly.
                if reference_index < 0 || reference_index as usize >= context.reference.len() {
                    return None;
                }

                (reference_index as usize, primary_index)
            }
        };

        debug_assert!(reference_index != usize::MAX, "{self:?}");
        debug_assert!(query_index != usize::MAX, "{self:?}");
        debug_assert!(reference_index < isize::MAX as usize, "{self:?}");
        debug_assert!(query_index < isize::MAX as usize, "{self:?}");

        Some(self.generate_successor(
            Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: -context.config.right_flank_length,
                data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::generate_successor_identifier_primary_extra_data(identifier, AlignmentType::PrimaryReentry, context),
            },
            cost_increment,
            AlignmentType::PrimaryReentry,
            context,
        ))
    }

    fn generate_template_switch_shortcut_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        delta_reference: isize,
        delta_query: isize,
        cost_increment: Strategies::Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Strategies::Cost::max_value() {
            return None;
        }

        let identifier @ (Identifier::Primary {
            reference_index,
            query_index,
            flank_index,
            ..
        }
        | Identifier::PrimaryReentry {
            reference_index,
            query_index,
            flank_index,
            ..
        }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        assert_eq!(flank_index, context.config.left_flank_length);

        let reference_index =
            usize::try_from(isize::try_from(reference_index).unwrap() + delta_reference).ok()?;
        let query_index =
            usize::try_from(isize::try_from(query_index).unwrap() + delta_query).ok()?;

        if reference_index >= context.reference.len() || query_index >= context.query.len() {
            return None;
        }

        Some(self.generate_successor(
            Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: -context.config.right_flank_length,
                data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::generate_successor_identifier_primary_extra_data(identifier, AlignmentType::PrimaryReentry, context),
            },
            cost_increment,
            AlignmentType::PrimaryShortcut {
                delta_reference,
                delta_query,
            },
            context,
        ))
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        cost_increment: Strategies::Cost,
        alignment_type: AlignmentType,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            node_data: self.node_data.generate_successor(
                identifier,
                cost_increment,
                alignment_type,
            ),
            strategies: self
                .strategies
                .generate_successor(identifier, alignment_type, context),
        }
    }
}

impl<PrimaryExtraData: Copy, Cost: AStarCost> NodeData<PrimaryExtraData, Cost> {
    fn create_root(identifier: Identifier<PrimaryExtraData>) -> Self {
        Self {
            identifier,
            predecessor: None,
            predecessor_edge_type: AlignmentType::Root,
            cost: Cost::zero(),
            a_star_lower_bound: Cost::zero(),
        }
    }

    fn generate_successor(
        &self,
        identifier: Identifier<PrimaryExtraData>,
        cost_increment: Cost,
        alignment_type: AlignmentType,
    ) -> Self {
        let cost = self.cost.checked_add(&cost_increment).unwrap();
        let a_star_lower_bound = self.a_star_lower_bound.saturating_sub(&cost_increment);
        Self {
            identifier,
            predecessor: Some(self.identifier),
            predecessor_edge_type: alignment_type,
            cost,
            a_star_lower_bound,
        }
    }

    fn lower_bound_cost(&self) -> Cost {
        self.cost + self.a_star_lower_bound
    }
}

impl<Strategies: AlignmentStrategySelector> Ord for Node<Strategies> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.strategies.node_ord_strategy.cmp(self, other)
    }
}

impl<Strategies: AlignmentStrategySelector> PartialOrd for Node<Strategies> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Strategies: AlignmentStrategySelector> Display for Node<Strategies>
where
    AlignmentStrategiesNodeMemory<Strategies>: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            node_data:
                NodeData {
                    identifier,
                    predecessor,
                    predecessor_edge_type,
                    cost,
                    a_star_lower_bound,
                },
            strategies,
        } = self;
        write!(f, "{identifier}; ")?;
        if let Some(predecessor) = predecessor {
            write!(f, "predecessor: {predecessor}; ")?;
        }
        write!(f, "alignment_type: {predecessor_edge_type}; ")?;
        write!(f, "cost: {cost} + {a_star_lower_bound}; ")?;
        write!(f, "strategies: {strategies}")
    }
}
