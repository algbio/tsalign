use std::iter;

use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchSecondary,
};
use log::{trace, warn};
use tagged_vec::TaggedVec;

use crate::ts_arrangement::character::Char;

use super::{
    complement::TsComplementArrangement,
    index_types::{ArrangementColumn, SourceColumn, TsInnerIdentifier},
    source::{RemovedHiddenChars, SourceChar, TsSourceArrangement},
    template_switch::TemplateSwitch,
};

pub struct TsInnerArrangement {
    inners: TaggedVec<TsInnerIdentifier, TsInner>,
}

pub struct TsInner {
    sequence: TaggedVec<ArrangementColumn, InnerChar>,
    template_switch: TemplateSwitch,
    reference: bool,
    complement: bool,
}

#[derive(Debug, Clone, Copy)]
pub enum InnerChar {
    Inner {
        column: SourceColumn,
        lower_case: bool,
        copy_depth: Option<usize>,
    },
    OptionalInner {
        column: SourceColumn,
        lower_case: bool,
        copy_depth: Option<usize>,
    },
    Gap {
        copy_depth: Option<usize>,
    },
    Blank,
}

impl TsInnerArrangement {
    pub fn new(
        source_arrangement: &mut TsSourceArrangement,
        complement_arrangement: &mut TsComplementArrangement,
        template_switches: Vec<TemplateSwitch>,
        visualise_equal_cost_ranges: bool,
    ) -> Self {
        let mut result = Self {
            inners: Default::default(),
        };

        trace!(
            "source width: {}; complement width: {}",
            source_arrangement.width(),
            complement_arrangement.width(),
        );

        for ts in template_switches {
            trace!(
                "source_inner: {:?}",
                ts.inner
                    .iter()
                    .map(|c| format!("{}", c.source_column()))
                    .collect::<Vec<_>>()
            );
            trace!("inner_alignment: {:?}", ts.inner_alignment.cigar());

            let (mut sp2_secondary, mut sp3_secondary) = match ts.secondary {
                TemplateSwitchSecondary::Reference => (
                    source_arrangement
                        .try_reference_source_to_arrangement_column(ts.sp2_secondary)
                        .unwrap_or_else(|| {
                            trace!("SP2 is at reference end");
                            source_arrangement.reference().len().into()
                        }),
                    source_arrangement
                        .try_reference_source_to_arrangement_column(ts.sp3_secondary)
                        .unwrap_or_else(|| {
                            trace!("SP3 is at reference end");
                            source_arrangement.reference().len().into()
                        }),
                ),
                TemplateSwitchSecondary::Query => (
                    source_arrangement
                        .try_query_source_to_arrangement_column(ts.sp2_secondary)
                        .unwrap_or_else(|| {
                            trace!("SP2 is at query end");
                            source_arrangement.query().len().into()
                        }),
                    source_arrangement
                        .try_query_source_to_arrangement_column(ts.sp3_secondary)
                        .unwrap_or_else(|| {
                            trace!("SP3 is at query end");
                            source_arrangement.query().len().into()
                        }),
                ),
            };
            let forward = sp2_secondary < sp3_secondary;
            trace!(
                "sp2_secondary: {} -> {sp2_secondary}; sp3_secondary: {} -> {sp3_secondary}; {}",
                ts.sp2_secondary,
                ts.sp3_secondary,
                if forward { "forward" } else { "reverse" }
            );

            let mut source_inner = ts.inner.iter().copied();
            let mut inner = TaggedVec::<ArrangementColumn, _>::default();
            inner.extend(iter::repeat_n(
                InnerChar::Blank,
                sp3_secondary.min(sp2_secondary).into(),
            ));
            let mut current_arrangement_column = sp3_secondary.min(sp2_secondary);
            debug_assert!(current_arrangement_column.primitive() < source_arrangement.width());
            debug_assert!(current_arrangement_column.primitive() < complement_arrangement.width());

            if forward {
                // Align inner against source.
                for alignment_type in ts.inner_alignment.iter_flat_cloned() {
                    match alignment_type {
                        AlignmentType::SecondaryInsertion => {
                            let is_gap = loop {
                                if source_arrangement.secondary(ts.secondary).len()
                                    <= current_arrangement_column.primitive()
                                {
                                    break false;
                                }

                                let c = source_arrangement.secondary(ts.secondary)
                                    [current_arrangement_column];

                                if c.is_gap() || c.is_source_char() {
                                    break c.is_gap();
                                }

                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            };

                            if !is_gap {
                                source_arrangement.insert_secondary_gap_with_minimum_copy_depth(
                                    ts.secondary,
                                    current_arrangement_column,
                                );

                                complement_arrangement.insert_blank(current_arrangement_column);
                                for existing_inner in result.inners.iter_values_mut() {
                                    existing_inner
                                        .sequence
                                        .insert(current_arrangement_column, InnerChar::Blank);
                                }

                                sp3_secondary += 1;
                            }

                            inner.push(source_inner.next().unwrap().into());
                            current_arrangement_column += 1;
                        }
                        AlignmentType::SecondaryDeletion => {
                            while !source_arrangement.secondary(ts.secondary)
                                [current_arrangement_column]
                                .is_source_char()
                            {
                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            }

                            inner.push(InnerChar::Gap {
                                copy_depth: source_arrangement.secondary(ts.secondary)
                                    [current_arrangement_column]
                                    .copy_depth(),
                            });
                            current_arrangement_column += 1;
                        }
                        AlignmentType::SecondarySubstitution | AlignmentType::SecondaryMatch => {
                            while !source_arrangement.secondary(ts.secondary)
                                [current_arrangement_column]
                                .is_source_char()
                            {
                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            }

                            let mut inner_char: InnerChar = source_inner.next().unwrap().into();
                            if alignment_type == AlignmentType::SecondarySubstitution {
                                source_arrangement.secondary_to_lower_case(
                                    ts.secondary,
                                    current_arrangement_column,
                                );
                                inner_char.to_lower_case();
                            }

                            inner.push(inner_char);
                            current_arrangement_column += 1;
                        }
                        _ => unreachable!(),
                    }
                }

                // We skip further secondary non-source chars for the assertion below.
                while current_arrangement_column.primitive() < source_arrangement.width()
                    && !source_arrangement.secondary(ts.secondary)[current_arrangement_column]
                        .is_source_char()
                {
                    current_arrangement_column += 1;
                }
                assert_eq!(current_arrangement_column, sp3_secondary);
            } else {
                // Align inner against source complement in reverse.
                let mut source_inner = source_inner.rev();
                for alignment_type in ts.inner_alignment.iter_flat_cloned().rev() {
                    trace!("Processing inner alignment {alignment_type}");

                    match alignment_type {
                        AlignmentType::SecondaryInsertion => {
                            let is_gap = loop {
                                if complement_arrangement
                                    .secondary_complement(ts.secondary)
                                    .len()
                                    <= current_arrangement_column.primitive()
                                {
                                    break false;
                                }

                                let c = complement_arrangement.secondary_complement(ts.secondary)
                                    [current_arrangement_column];

                                if c.is_gap() || c.is_source_char() {
                                    break c.is_gap();
                                }

                                trace!("Skipping inner character");
                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            };

                            if !is_gap {
                                complement_arrangement.insert_secondary_complement_gap(
                                    ts.secondary,
                                    current_arrangement_column,
                                );

                                source_arrangement.insert_blank(current_arrangement_column);
                                for existing_inner in result.inners.iter_values_mut() {
                                    existing_inner
                                        .sequence
                                        .insert(current_arrangement_column, InnerChar::Blank);
                                }

                                sp2_secondary += 1;
                            }

                            inner.push(source_inner.next().unwrap().into());
                            current_arrangement_column += 1;
                        }
                        AlignmentType::SecondaryDeletion => {
                            while !complement_arrangement.secondary_complement(ts.secondary)
                                [current_arrangement_column]
                                .is_source_char()
                            {
                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            }

                            complement_arrangement
                                .show_secondary_character(ts.secondary, current_arrangement_column);
                            inner.push(InnerChar::Gap {
                                copy_depth: source_arrangement.secondary(ts.secondary)
                                    [current_arrangement_column]
                                    .copy_depth(),
                            });
                            current_arrangement_column += 1;
                        }
                        AlignmentType::SecondarySubstitution | AlignmentType::SecondaryMatch => {
                            while !source_arrangement.secondary(ts.secondary)
                                [current_arrangement_column]
                                .is_source_char()
                            {
                                inner.push(InnerChar::Blank);
                                current_arrangement_column += 1;
                            }

                            complement_arrangement
                                .show_secondary_character(ts.secondary, current_arrangement_column);

                            let mut inner_char: InnerChar = source_inner.next().unwrap().into();
                            if alignment_type == AlignmentType::SecondarySubstitution {
                                complement_arrangement.secondary_to_lower_case(
                                    ts.secondary,
                                    current_arrangement_column,
                                );
                                inner_char.to_lower_case();
                            }

                            inner.push(inner_char);
                            current_arrangement_column += 1;
                        }
                        _ => unreachable!(),
                    }
                }

                // We skip further secondary non-source chars for the assertion below.
                while current_arrangement_column.primitive() < source_arrangement.width()
                    && !source_arrangement.secondary(ts.secondary)[current_arrangement_column]
                        .is_source_char()
                {
                    current_arrangement_column += 1;
                }
                assert_eq!(current_arrangement_column, sp2_secondary);
            }

            let suffix_blanks =
                iter::repeat_n(InnerChar::Blank, source_arrangement.reference().len())
                    .skip(inner.len());
            inner.extend(suffix_blanks);

            if visualise_equal_cost_ranges {
                // Add characters to visualise TSM equal cost range.
                if forward {
                    warn!("TSM equal cost range visualisation is not implemented for forward TSMs.")
                } else {
                    // Insert range characters before point 2.

                    let last_initial_blank = inner
                        .iter()
                        .take_while(|(_, c)| c.is_blank())
                        .last()
                        .map(|(i, _)| i);
                    let first_final_blank = inner
                        .iter()
                        .rev()
                        .take_while(|(_, c)| c.is_blank())
                        .last()
                        .map(|(i, _)| i)
                        .unwrap_or(inner.len().into());
                    let first_non_blank = inner
                        .iter()
                        .find(|(_, c)| !c.is_blank())
                        .map(|(i, _)| i)
                        .unwrap();
                    let first_source_column = inner
                        .iter_values()
                        .filter(|c| c.is_source_char())
                        .map(|c| c.source_column())
                        .next()
                        .unwrap();
                    let last_source_column = inner
                        .iter_values()
                        .rev()
                        .filter(|c| c.is_source_char())
                        .map(|c| c.source_column())
                        .next()
                        .unwrap();
                    assert!(first_source_column >= last_source_column);

                    // Add prefix to extend to max_end.
                    let mut arrangement_column =
                        last_initial_blank.map(|i| i + 1usize).unwrap_or(0.into());
                    let mut source_column = first_source_column;
                    #[allow(clippy::explicit_counter_loop)]
                    for _ in 0..ts.equal_cost_range.max_end {
                        arrangement_column -= 1;
                        source_column += 1;

                        inner[arrangement_column] = InnerChar::OptionalInner {
                            column: source_column,
                            lower_case: false,
                            copy_depth: None,
                        };
                    }

                    // Add suffix to extend to min_start.
                    let mut arrangement_column = first_final_blank - 1usize;
                    let mut source_column = last_source_column;
                    for _ in 0..-ts.equal_cost_range.min_start {
                        arrangement_column += 1;
                        source_column -= 1;

                        inner[arrangement_column] = InnerChar::OptionalInner {
                            column: source_column,
                            lower_case: false,
                            copy_depth: None,
                        };
                    }

                    // Convert prefix to extend to min_end.
                    let mut arrangement_column = first_non_blank;
                    for _ in 0..-ts.equal_cost_range.min_end {
                        while !inner[arrangement_column].is_source_char() {
                            arrangement_column += 1;
                        }

                        inner[arrangement_column].to_optional();
                        arrangement_column += 1;
                    }

                    // Convert suffix to extend to max_start.
                    let mut arrangement_column = first_final_blank;
                    for _ in 0..ts.equal_cost_range.max_start {
                        arrangement_column -= 1;

                        while !inner[arrangement_column].is_source_char() {
                            arrangement_column -= 1;
                        }

                        inner[arrangement_column].to_optional();
                    }
                }
            }

            let is_reference = match ts.secondary {
                TemplateSwitchSecondary::Reference => true,
                TemplateSwitchSecondary::Query => false,
            };
            result
                .inners
                .push(TsInner::new(inner, ts, is_reference, !forward));
        }

        result
    }

    pub fn remove_columns(
        &mut self,
        columns: impl IntoIterator<Item = ArrangementColumn> + Clone,
        removed_hidden_chars: &RemovedHiddenChars,
    ) {
        for inner in self.inners.iter_values_mut() {
            inner.sequence.remove_multi(columns.clone());
            inner
                .template_switch
                .remove_hidden_chars(removed_hidden_chars);
        }
    }

    pub fn inners(&self) -> &TaggedVec<TsInnerIdentifier, TsInner> {
        &self.inners
    }

    pub fn reference_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inners
            .iter()
            .filter(|inner| inner.1.reference && !inner.1.complement)
    }

    pub fn query_inners(&self) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inners
            .iter()
            .filter(|inner| !inner.1.reference && !inner.1.complement)
    }

    pub fn reference_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inners
            .iter()
            .filter(|inner| inner.1.reference && inner.1.complement)
    }

    pub fn query_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inners
            .iter()
            .filter(|inner| !inner.1.reference && inner.1.complement)
    }

    pub fn inner_first_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        let sequence = &self.inners[inner_identifier].sequence;

        sequence
            .iter()
            .find(|(_, c)| !c.is_blank())
            .map(|(i, _)| i)
            .unwrap_or(sequence.len().into())
    }

    pub fn inner_last_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        let sequence = &self.inners[inner_identifier].sequence;

        sequence
            .iter()
            .rev()
            .find(|(_, c)| !c.is_blank())
            .map(|(i, _)| i)
            .unwrap() // If None, would need to return -1 here, but return value is unsigned.
    }
}

impl TsInner {
    fn new(
        sequence: TaggedVec<ArrangementColumn, InnerChar>,
        template_switch: TemplateSwitch,
        reference: bool,
        complement: bool,
    ) -> Self {
        Self {
            sequence,
            template_switch,
            reference,
            complement,
        }
    }

    pub fn sequence(&self) -> &TaggedVec<ArrangementColumn, InnerChar> {
        &self.sequence
    }

    pub fn template_switch(&self) -> &TemplateSwitch {
        &self.template_switch
    }
}

impl InnerChar {
    pub fn to_lower_case(&mut self) {
        match self {
            Self::Inner { lower_case, .. } | Self::OptionalInner { lower_case, .. } => {
                *lower_case = true
            }
            Self::Gap { .. } | Self::Blank => panic!("Not lowercasable"),
        }
    }

    pub fn to_optional(&mut self) {
        match *self {
            Self::Inner {
                column,
                lower_case,
                copy_depth,
            }
            | Self::OptionalInner {
                column,
                lower_case,
                copy_depth,
            } => {
                *self = Self::OptionalInner {
                    column,
                    lower_case,
                    copy_depth,
                }
            }
            Self::Gap { .. } | Self::Blank => panic!("Not optionalisable"),
        }
    }
}

impl From<SourceChar> for InnerChar {
    fn from(value: SourceChar) -> Self {
        match value {
            SourceChar::Source {
                column,
                lower_case,
                copy_depth,
            } => Self::Inner {
                column,
                lower_case,
                copy_depth,
            },
            SourceChar::Hidden { .. } => {
                panic!("Cannot be translated into InnerChar")
            }
            SourceChar::Gap { copy_depth } => Self::Gap { copy_depth },
            SourceChar::Separator | SourceChar::Spacer | SourceChar::Blank => Self::Blank,
        }
    }
}

impl Char for InnerChar {
    fn source_column(&self) -> SourceColumn {
        match self {
            Self::Inner { column, .. } | Self::OptionalInner { column, .. } => *column,
            Self::Gap { .. } | Self::Blank => panic!("Has no source column"),
        }
    }

    fn is_char(&self) -> bool {
        matches!(self, Self::Inner { .. } | Self::OptionalInner { .. })
    }

    fn is_gap(&self) -> bool {
        matches!(self, Self::Gap { .. })
    }

    fn is_spacer(&self) -> bool {
        false
    }

    fn is_blank(&self) -> bool {
        matches!(self, Self::Blank)
    }

    fn is_source_char(&self) -> bool {
        self.is_char()
    }

    fn is_hidden(&self) -> bool {
        false
    }
}
