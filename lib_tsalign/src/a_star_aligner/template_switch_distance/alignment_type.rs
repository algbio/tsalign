use equal_cost_range::EqualCostRange;

use crate::a_star_aligner::alignment_result::IAlignmentType;

use super::identifier::{
    TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection,
};

pub mod equal_cost_range;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AlignmentType {
    /// The query contains a base that is missing from the reference.
    PrimaryInsertion,
    /// The query is missing a base present in the reference.
    PrimaryDeletion,
    /// The query contains a different base than the reference.
    PrimarySubstitution,
    /// The query contains the same base as the reference.
    PrimaryMatch,
    /// The query contains a base that is missing from the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankInsertion,
    /// The query is missing a base present in the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankDeletion,
    /// The query contains a different base than the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankSubstitution,
    /// The query contains the same base as the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankMatch,
    /// The TS ancestor is missing a base present in the TS descendant.
    SecondaryInsertion,
    /// The TS ancestor contains a base that is missing from the TS descendant.
    SecondaryDeletion,
    /// The TS ancestor contains a different base than the TS descendant.
    SecondarySubstitution,
    /// The TS ancestor contains the same base as the TS descendant.
    SecondaryMatch,
    /// A template switch entrance.
    TemplateSwitchEntrance {
        first_offset: isize,
        equal_cost_range: EqualCostRange,
        descendant: TemplateSwitchDescendant,
        ancestor: TemplateSwitchAncestor,
        direction: TemplateSwitchDirection,
    },
    /// A template switch exit.
    TemplateSwitchExit {
        /// The number of characters that are skipped on the anti-descendant sequence.
        /// If negative, it is the number of characters that are repeated on the anti-descendant sequence.
        ///
        /// In terms of switchpoints, this is the difference `SP4 - SP1`.
        ///
        /// Note that the anti-descendant sequence is not necessarily equal to the ancestor sequence.
        anti_descendant_gap: isize,
    },
    /// This node is the root node, hence it was not generated via alignment.
    Root,
    /// The root node of a secondary graph.
    SecondaryRoot,
    /// A reentry node into the primary graph, treated like a root.
    PrimaryReentry,
    /// A shortcut in the primary matrix.
    ///
    /// Used only for computing lower bounds.
    PrimaryShortcut {
        delta_reference: isize,
        delta_query: isize,
    },
}

impl IAlignmentType for AlignmentType {
    fn is_repeatable(&self) -> bool {
        match self {
            Self::PrimaryInsertion
            | Self::PrimaryFlankInsertion
            | Self::SecondaryInsertion
            | Self::PrimaryDeletion
            | Self::PrimaryFlankDeletion
            | Self::SecondaryDeletion
            | Self::PrimarySubstitution
            | Self::PrimaryFlankSubstitution
            | Self::SecondarySubstitution
            | Self::PrimaryMatch
            | Self::PrimaryFlankMatch
            | Self::SecondaryMatch
            | Self::Root
            | Self::SecondaryRoot
            | Self::PrimaryReentry => true,
            Self::TemplateSwitchEntrance { .. }
            | Self::TemplateSwitchExit { .. }
            | Self::PrimaryShortcut { .. } => false,
        }
    }

    fn is_repeated(&self, previous: &Self) -> bool {
        if (matches!(self, Self::PrimaryInsertion | Self::PrimaryFlankInsertion)
            && matches!(
                previous,
                Self::PrimaryInsertion | Self::PrimaryFlankInsertion
            ))
            || (matches!(self, Self::PrimaryDeletion | Self::PrimaryFlankDeletion)
                && matches!(previous, Self::PrimaryDeletion | Self::PrimaryFlankDeletion))
            || (matches!(
                self,
                Self::PrimarySubstitution | Self::PrimaryFlankSubstitution
            ) && matches!(
                previous,
                Self::PrimarySubstitution | Self::PrimaryFlankSubstitution
            ))
            || (matches!(self, Self::PrimaryMatch | Self::PrimaryFlankMatch)
                && matches!(previous, Self::PrimaryMatch | Self::PrimaryFlankMatch))
        {
            return true;
        }

        match (self, previous) {
            (
                Self::TemplateSwitchEntrance {
                    descendant: descendant_a,
                    ancestor: ancestor_a,
                    ..
                },
                Self::TemplateSwitchEntrance {
                    descendant: descendant_b,
                    ancestor: ancestor_b,
                    ..
                },
            ) => descendant_a == descendant_b && ancestor_a == ancestor_b,
            (Self::TemplateSwitchExit { .. }, Self::TemplateSwitchExit { .. }) => true,
            (Self::PrimaryShortcut { .. }, Self::PrimaryShortcut { .. }) => false,
            (a, b) => a == b,
        }
    }

    fn is_internal(&self) -> bool {
        matches!(
            self,
            Self::Root | Self::SecondaryRoot | Self::PrimaryReentry
        )
    }

    fn is_template_switch_entrance(&self) -> bool {
        matches!(self, Self::TemplateSwitchEntrance { .. })
    }

    fn is_template_switch_exit(&self) -> bool {
        matches!(self, Self::TemplateSwitchExit { .. })
    }
}

impl AlignmentType {
    pub fn inverted(&self) -> Self {
        match self {
            Self::PrimaryInsertion => Self::PrimaryDeletion,
            Self::PrimaryDeletion => Self::PrimaryInsertion,
            Self::PrimaryFlankInsertion => Self::PrimaryFlankDeletion,
            Self::PrimaryFlankDeletion => Self::PrimaryFlankInsertion,
            Self::SecondaryInsertion => Self::SecondaryDeletion,
            Self::SecondaryDeletion => Self::SecondaryInsertion,
            Self::TemplateSwitchEntrance {
                descendant,
                ancestor,
                direction,
                equal_cost_range,
                first_offset,
            } => Self::TemplateSwitchEntrance {
                descendant: descendant.inverted(),
                ancestor: ancestor.inverted(),
                direction: direction.inverted(),
                equal_cost_range: *equal_cost_range,
                first_offset: *first_offset,
            },
            Self::PrimaryShortcut {
                delta_reference,
                delta_query,
            } => Self::PrimaryShortcut {
                delta_reference: *delta_query,
                delta_query: *delta_reference,
            },

            symmetric @ (Self::PrimarySubstitution
            | Self::PrimaryMatch
            | Self::PrimaryFlankSubstitution
            | Self::PrimaryFlankMatch
            | Self::SecondarySubstitution
            | Self::SecondaryMatch
            | Self::TemplateSwitchExit { .. }
            | Self::Root
            | Self::SecondaryRoot
            | Self::PrimaryReentry) => *symmetric,
        }
    }
}
