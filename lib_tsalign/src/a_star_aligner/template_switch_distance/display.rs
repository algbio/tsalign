use std::fmt::{Display, Formatter, Result};

use super::{
    AlignmentType, GapType, Identifier, TemplateSwitchAncestor, TemplateSwitchDescendant,
    alignment_type::equal_cost_range::EqualCostRange, identifier::TemplateSwitchDirection,
};

impl Display for AlignmentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::PrimaryInsertion | Self::PrimaryFlankInsertion | Self::SecondaryInsertion => {
                write!(f, "I")
            }
            Self::PrimaryDeletion | Self::PrimaryFlankDeletion | Self::SecondaryDeletion => {
                write!(f, "D")
            }
            Self::PrimarySubstitution
            | Self::PrimaryFlankSubstitution
            | Self::SecondarySubstitution => write!(f, "X"),
            Self::PrimaryMatch | Self::PrimaryFlankMatch | Self::SecondaryMatch => write!(f, "="),
            Self::TemplateSwitchEntrance {
                descendant,
                ancestor,
                direction,
                equal_cost_range,
                first_offset,
            } => write!(
                f,
                "[TS{descendant}{ancestor}{direction}:{equal_cost_range}:{first_offset}:"
            ),
            Self::TemplateSwitchExit {
                anti_descendant_gap,
            } => write!(f, ":{anti_descendant_gap}]"),
            Self::Root => Ok(()),
            Self::SecondaryRoot => Ok(()),
            Self::PrimaryReentry => Ok(()),
            Self::PrimaryShortcut {
                delta_reference,
                delta_query,
            } => write!(f, "[PS:R{delta_reference}Q{delta_query}]"),
        }
    }
}

impl Display for GapType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Insertion => write!(f, "I"),
            Self::Deletion => write!(f, "D"),
            Self::None => write!(f, "=/X"),
        }
    }
}

impl Display for TemplateSwitchDescendant {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Reference => write!(f, "R"),
            Self::Query => write!(f, "Q"),
        }
    }
}

impl Display for TemplateSwitchAncestor {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Reference => write!(f, "R"),
            Self::Query => write!(f, "Q"),
        }
    }
}

impl Display for TemplateSwitchDirection {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Forward => write!(f, "F"),
            Self::Reverse => write!(f, "R"),
        }
    }
}

impl Display for EqualCostRange {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.is_valid() {
            let Self {
                min_start,
                max_start,
                min_end,
                max_end,
            } = self;
            write!(f, "[{min_start},{max_start}]:[{min_end},{max_end}]")
        } else {
            write!(f, "[-]:[-]")
        }
    }
}

impl<PrimaryExtraData> Display for Identifier<PrimaryExtraData> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
                ..
            } => write!(
                f,
                "Primary({reference_index}R, {query_index}Q, {flank_index}F, {gap_type})",
            ),

            Self::PrimaryReentry {
                reference_index,
                query_index,
                flank_index,
                gap_type,
                ..
            } => write!(
                f,
                "PrimaryReentry({reference_index}R, {query_index}Q, {flank_index}F, {gap_type})",
            ),

            Self::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                template_switch_first_offset,
            } => {
                write!(
                    f,
                    "TemplateSwitchEntrance({entrance_reference_index}R, {entrance_query_index}Q, {template_switch_descendant}D, {template_switch_ancestor}A, {template_switch_direction}D, {template_switch_first_offset}O)",
                )
            }

            #[allow(clippy::uninlined_format_args)]
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                length,
                descendant_index,
                ancestor_index,
                gap_type,
            } => write!(
                f,
                "Secondary({}R, {}Q, {}L, {}P, {}S, {}, {}, {}, {})",
                entrance_reference_index,
                entrance_query_index,
                length,
                descendant_index,
                ancestor_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                gap_type,
            ),

            #[allow(clippy::uninlined_format_args)]
            Self::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
                descendant_index: primary_index,
                anti_descendant_gap: anti_primary_gap,
            } => write!(
                f,
                "TemplateSwitchExit({}R, {}Q, {}P, {}G, {}, {}, {})",
                entrance_reference_index,
                entrance_query_index,
                primary_index,
                anti_primary_gap,
                template_switch_descendant,
                template_switch_ancestor,
                template_switch_direction,
            ),
        }
    }
}
