use std::{fmt::Debug, hash::Hash};

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::AStarIdentifier;

use super::{
    AlignmentType, Context,
    strategies::{AlignmentStrategySelector, primary_match::PrimaryMatchStrategy},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub enum Identifier<PrimaryExtraData> {
    Primary {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
        /// Positive for left flank, negative for right flank.
        flank_index: isize,
        data: PrimaryExtraData,
    },
    PrimaryReentry {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
        /// Positive for left flank, negative for right flank.
        flank_index: isize,
        data: PrimaryExtraData,
    },
    TemplateSwitchEntrance {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_descendant: TemplateSwitchDescendant,
        template_switch_ancestor: TemplateSwitchAncestor,
        template_switch_direction: TemplateSwitchDirection,
        template_switch_first_offset: isize,
    },
    Secondary {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_descendant: TemplateSwitchDescendant,
        template_switch_ancestor: TemplateSwitchAncestor,
        template_switch_direction: TemplateSwitchDirection,
        length: usize,
        /// The index that does not jump.
        descendant_index: usize,
        /// The index that jumps.
        ancestor_index: usize,
        gap_type: GapType,
    },
    TemplateSwitchExit {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_descendant: TemplateSwitchDescendant,
        template_switch_ancestor: TemplateSwitchAncestor,
        template_switch_direction: TemplateSwitchDirection,
        /// The index that does not jump.
        descendant_index: usize,
        anti_descendant_gap: isize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub enum GapType {
    Insertion,
    Deletion,
    None,
}

/// The descendant sequence is the sequence for which the template switch does not jump.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TemplateSwitchDescendant {
    Reference,
    Query,
}

/// The ancestor sequence is the sequence for which the template switch jumps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TemplateSwitchAncestor {
    Reference,
    Query,
}

/// A reverse TSM is a TSM, while a forward TSM is a repeat.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TemplateSwitchDirection {
    Forward,
    Reverse,
}

impl<PrimaryExtraData> Identifier<PrimaryExtraData> {
    pub const fn new_primary(
        reference_index: usize,
        query_index: usize,
        flank_index: isize,
        gap_type: GapType,
        data: PrimaryExtraData,
    ) -> Self {
        debug_assert!(reference_index != usize::MAX);
        debug_assert!(query_index != usize::MAX);
        debug_assert!(reference_index < isize::MAX as usize);
        debug_assert!(query_index < isize::MAX as usize);

        Self::Primary {
            reference_index,
            query_index,
            flank_index,
            gap_type,
            data,
        }
    }

    pub fn generate_primary_diagonal_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>,
        PrimaryMatch: PrimaryMatchStrategy<Strategies::Cost, IdentifierPrimaryExtraData = PrimaryExtraData>,
    >(
        self,
        flank_index: isize,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index: reference_index + 1,
                    query_index: query_index + 1,
                    flank_index,
                    gap_type: GapType::None,
                    data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                    <Strategies as AlignmentStrategySelector>::Cost,
                >>::generate_successor_identifier_primary_extra_data(self, alignment_type, context),
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_primary_deletion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>,
        PrimaryMatch: PrimaryMatchStrategy<Strategies::Cost, IdentifierPrimaryExtraData = PrimaryExtraData>,
    >(
        self,
        flank_index: isize,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index: reference_index + 1,
                    query_index,
                    flank_index,
                    gap_type: GapType::Deletion,
                    data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                    <Strategies as AlignmentStrategySelector>::Cost,
                >>::generate_successor_identifier_primary_extra_data(self, alignment_type, context),
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_primary_insertion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>,
        PrimaryMatch: PrimaryMatchStrategy<Strategies::Cost, IdentifierPrimaryExtraData = PrimaryExtraData>,
    >(
        self,
        flank_index: isize,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index,
                    query_index: query_index + 1,
                    flank_index,
                    gap_type: GapType::Insertion,
                    data: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                    <Strategies as AlignmentStrategySelector>::Cost,
                >>::generate_successor_identifier_primary_extra_data(self, alignment_type, context),
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_initial_template_switch_entrance_successors(
        self,
    ) -> impl Iterator<Item = Self> {
        match self {
            Identifier::Primary {
                reference_index: entrance_reference_index,
                query_index: entrance_query_index,
                ..
            }
            | Identifier::PrimaryReentry {
                reference_index: entrance_reference_index,
                query_index: entrance_query_index,
                ..
            } => [
                TemplateSwitchDescendant::Reference,
                TemplateSwitchDescendant::Reference,
                TemplateSwitchDescendant::Query,
                TemplateSwitchDescendant::Query,
                TemplateSwitchDescendant::Reference,
                TemplateSwitchDescendant::Reference,
                TemplateSwitchDescendant::Query,
                TemplateSwitchDescendant::Query,
            ]
            .into_iter()
            .zip([
                TemplateSwitchAncestor::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchAncestor::Query,
            ])
            .zip([
                TemplateSwitchDirection::Forward,
                TemplateSwitchDirection::Forward,
                TemplateSwitchDirection::Forward,
                TemplateSwitchDirection::Forward,
                TemplateSwitchDirection::Reverse,
                TemplateSwitchDirection::Reverse,
                TemplateSwitchDirection::Reverse,
                TemplateSwitchDirection::Reverse,
            ])
            .flat_map(
                move |(
                    (template_switch_descendant, template_switch_ancestor),
                    template_switch_direction,
                )| {
                    match template_switch_direction {
                        TemplateSwitchDirection::Forward => vec![
                            Identifier::TemplateSwitchEntrance {
                                entrance_reference_index,
                                entrance_query_index,
                                template_switch_descendant,
                                template_switch_ancestor,
                                template_switch_direction,
                                template_switch_first_offset: -1,
                            },
                            Identifier::TemplateSwitchEntrance {
                                entrance_reference_index,
                                entrance_query_index,
                                template_switch_descendant,
                                template_switch_ancestor,
                                template_switch_direction,
                                template_switch_first_offset: 1,
                            },
                        ],
                        TemplateSwitchDirection::Reverse => {
                            vec![Identifier::TemplateSwitchEntrance {
                                entrance_reference_index,
                                entrance_query_index,
                                template_switch_descendant,
                                template_switch_ancestor,
                                template_switch_direction,
                                template_switch_first_offset: 0,
                            }]
                        }
                    }
                },
            ),

            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_secondary_diagonal_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length,
                descendant_index,
                ancestor_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length: length + 1,
                descendant_index: descendant_index + 1,
                ancestor_index: match template_switch_direction {
                    TemplateSwitchDirection::Forward => ancestor_index + 1,
                    TemplateSwitchDirection::Reverse => ancestor_index - 1,
                },
                gap_type: GapType::None,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The ancestor contains a base missing in the descendant.
    pub fn generate_secondary_deletion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length,
                descendant_index,
                ancestor_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length,
                descendant_index,
                ancestor_index: match template_switch_direction {
                    TemplateSwitchDirection::Forward => ancestor_index + 1,
                    TemplateSwitchDirection::Reverse => ancestor_index - 1,
                },
                gap_type: GapType::Deletion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The ancestor misses a base present in the descendant.
    pub fn generate_secondary_insertion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length,
                descendant_index,
                ancestor_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length: length + 1,
                descendant_index: descendant_index + 1,
                ancestor_index,
                gap_type: GapType::Insertion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// Returns the anti-diagonal for variants where it exists, or [`usize::MAX`](core::primitive::usize::MAX) otherwise.
    pub fn anti_diagonal(self) -> usize {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                reference_index + query_index
            }
            _ => usize::MAX,
        }
    }
}

impl TemplateSwitchDescendant {
    pub fn inverted(&self) -> Self {
        match self {
            Self::Reference => Self::Query,
            Self::Query => Self::Reference,
        }
    }

    pub fn get<
        'reference: 'result,
        'query: 'result,
        'result,
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &'reference SubsequenceType,
        query: &'query SubsequenceType,
    ) -> &'result SubsequenceType {
        match self {
            TemplateSwitchDescendant::Reference => reference,
            TemplateSwitchDescendant::Query => query,
        }
    }
}

impl TemplateSwitchAncestor {
    pub fn inverted(&self) -> Self {
        match self {
            Self::Reference => Self::Query,
            Self::Query => Self::Reference,
        }
    }

    pub fn get<
        'reference: 'result,
        'query: 'result,
        'result,
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &'reference SubsequenceType,
        query: &'query SubsequenceType,
    ) -> &'result SubsequenceType {
        match self {
            TemplateSwitchAncestor::Reference => reference,
            TemplateSwitchAncestor::Query => query,
        }
    }
}

impl TemplateSwitchDirection {
    pub fn inverted(&self) -> Self {
        *self
    }
}

impl<PrimaryExtraData: Debug + Clone + Eq + Hash> AStarIdentifier for Identifier<PrimaryExtraData> {}
