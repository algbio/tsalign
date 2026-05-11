use std::{
    fmt::Debug,
    ops::{Bound, Range, RangeBounds},
};

use character::Char;
use complement::{ComplementChar, TsComplementArrangement};
use index_types::{ArrangementCharColumn, ArrangementColumn, SourceColumn, TsInnerIdentifier};
use inner::{TsInner, TsInnerArrangement};
use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchSecondary,
};
use log::debug;
use source::{SourceChar, TsSourceArrangement};
use tagged_vec::TaggedVec;
use template_switch::TemplateSwitch;

use crate::error::Result;

pub mod character;
pub mod complement;
pub mod index_types;
pub mod inner;
pub mod row;
pub mod source;
pub mod template_switch;

pub struct TsArrangement {
    source: TsSourceArrangement,
    complement: TsComplementArrangement,
    inner: TsInnerArrangement,
}

impl TsArrangement {
    pub fn new(
        reference_alignment_offset: usize,
        query_alignment_offset: usize,
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
    ) -> Result<Self> {
        let mut template_switches = Vec::new();
        let mut source = TsSourceArrangement::new(
            reference_alignment_offset,
            query_alignment_offset,
            reference_length,
            query_length,
            alignment,
            &mut template_switches,
        )?;
        let mut complement = TsComplementArrangement::new(&source);
        let inner = TsInnerArrangement::new(&mut source, &mut complement, template_switches);

        Ok(Self {
            source,
            complement,
            inner,
        })
    }

    /// Removes columns that contain only blanks and hidden characters.
    pub fn remove_empty_columns(&mut self) {
        let mut remove_columns = Vec::new();

        'column_iter: for column in self.reference().iter_indices() {
            if !self.reference()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.query()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.reference_complement()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.query_complement()[column].is_blank_or_hidden() {
                continue;
            }
            for inner in self.inner.inners().iter_values() {
                if !inner.sequence()[column].is_blank_or_hidden() {
                    continue 'column_iter;
                }
            }

            remove_columns.push(column);
        }

        let remove_columns = remove_columns;
        debug!("Removing columns: {remove_columns:?}");

        let removed_hidden_chars = self.source.remove_columns(remove_columns.iter().copied());
        self.complement
            .remove_columns(remove_columns.iter().copied());
        self.inner
            .remove_columns(remove_columns.iter().copied(), &removed_hidden_chars);
    }

    /// Unhides all characters of a complement, if any character of the complement is visible.
    pub fn show_complete_complements_if_used(&mut self) {
        let show = |sequence: &mut TaggedVec<ArrangementColumn, ComplementChar>| {
            if sequence.iter_values().any(Char::is_visible_char) {
                sequence.iter_values_mut().for_each(|c| {
                    if c.is_char() {
                        c.show()
                    }
                })
            }
        };

        show(self.complement.reference_complement_mut());
        show(self.complement.query_complement_mut());
    }

    /// Removes the columns from the arrangement that are further than `context_character_amount` characters away from the first and last TSM-related columns.
    /// Returns the source column ranges that are still included in the arrangement after the removal.
    /// The first range is the reference source column range, the second range is the query source column range.
    pub fn limit_context_to(
        &mut self,
        context_character_amount: usize,
    ) -> (Range<SourceColumn>, Range<SourceColumn>) {
        let first_inner_column = self.first_interesting_column();
        let last_inner_column = self.last_interesting_column();
        let first_printed_source_column_inclusive =
            first_inner_column.saturating_sub(context_character_amount);
        let last_printed_source_column_exclusive = last_inner_column
            .saturating_add(1usize)
            .saturating_add(context_character_amount);

        let result = (
            first_printed_source_column_inclusive
                ..last_printed_source_column_exclusive.min(self.source.reference_length().into()),
            first_printed_source_column_inclusive
                ..last_printed_source_column_exclusive.min(self.source.query_length().into()),
        );

        let first_printed_arrangement_column_inclusive = self
            .reference_source_to_arrangement_column(first_printed_source_column_inclusive)
            .min(self.query_source_to_arrangement_column(first_printed_source_column_inclusive));
        let last_printed_arrangement_column_exclusive = self
            .reference_source_to_arrangement_column(
                last_printed_source_column_exclusive.min(self.source.reference_length().into()),
            )
            .max(self.query_source_to_arrangement_column(
                last_printed_source_column_exclusive.min(self.source.query_length().into()),
            ));
        self.remove_column_range(last_printed_arrangement_column_exclusive..);
        self.remove_column_range(..first_printed_arrangement_column_inclusive);

        result
    }

    pub fn remove_column_range(&mut self, range: impl RangeBounds<ArrangementColumn> + Debug) {
        debug!("Removing column range: {range:?}");
        let range_start = match range.start_bound() {
            Bound::Included(inclusive) => *inclusive,
            Bound::Excluded(exclusive) => *exclusive + 1usize,
            Bound::Unbounded => 0usize.into(),
        };
        let range_end = match range.end_bound() {
            Bound::Included(inclusive) => *inclusive + 1usize,
            Bound::Excluded(exclusive) => *exclusive,
            Bound::Unbounded => self.width().into(),
        };
        let remove_columns = || (range_start.primitive()..range_end.primitive()).map(Into::into);

        let removed_hidden_chars = self.source.remove_columns(remove_columns());
        self.complement.remove_columns(remove_columns());
        self.inner
            .remove_columns(remove_columns(), &removed_hidden_chars);
    }

    pub fn reference(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        self.source.reference()
    }

    pub fn query(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        self.source.query()
    }

    pub fn reference_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        self.complement.reference_complement()
    }

    pub fn query_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        self.complement.query_complement()
    }

    pub fn inners(&self) -> &TaggedVec<TsInnerIdentifier, TsInner> {
        self.inner.inners()
    }

    pub fn reference_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.reference_inners()
    }

    pub fn query_inners(&self) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.query_inners()
    }

    pub fn reference_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.reference_complement_inners()
    }

    pub fn query_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.query_complement_inners()
    }

    pub fn template_switches(&self) -> impl Iterator<Item = (TsInnerIdentifier, &TemplateSwitch)> {
        self.inners()
            .iter()
            .map(|(identifier, inner)| (identifier, inner.template_switch()))
    }

    pub fn width(&self) -> usize {
        self.source.reference().len()
    }

    pub fn reference_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        self.source
            .reference_arrangement_char_to_arrangement_column(column)
    }

    pub fn query_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        self.source
            .query_arrangement_char_to_arrangement_column(column)
    }

    pub fn reference_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> ArrangementColumn {
        self.source.reference_source_to_arrangement_column(column)
    }

    pub fn query_source_to_arrangement_column(&self, column: SourceColumn) -> ArrangementColumn {
        self.source.query_source_to_arrangement_column(column)
    }

    pub fn reference_arrangement_to_source_column(
        &self,
        column: ArrangementColumn,
    ) -> SourceColumn {
        self.source.reference_arrangement_to_source_column(column)
    }

    pub fn query_arrangement_to_source_column(&self, column: ArrangementColumn) -> SourceColumn {
        self.source.query_arrangement_to_source_column(column)
    }

    pub fn try_reference_arrangement_to_source_column(
        &self,
        column: ArrangementColumn,
    ) -> Option<SourceColumn> {
        self.source
            .try_reference_arrangement_to_source_column(column)
    }

    pub fn try_query_arrangement_to_source_column(
        &self,
        column: ArrangementColumn,
    ) -> Option<SourceColumn> {
        self.source.try_query_arrangement_to_source_column(column)
    }

    pub fn secondary_source_to_arrangement_column(
        &self,
        column: SourceColumn,
        secondary: TemplateSwitchSecondary,
    ) -> ArrangementColumn {
        match secondary {
            TemplateSwitchSecondary::Reference => {
                self.source.reference_source_to_arrangement_column(column)
            }
            TemplateSwitchSecondary::Query => {
                self.source.query_source_to_arrangement_column(column)
            }
        }
    }

    pub fn reference_arrangement_char_to_source_column(
        &self,
        column: ArrangementCharColumn,
    ) -> SourceColumn {
        self.source
            .reference_arrangement_char_to_source_column(column)
    }

    pub fn query_arrangement_char_to_source_column(
        &self,
        column: ArrangementCharColumn,
    ) -> SourceColumn {
        self.source.query_arrangement_char_to_source_column(column)
    }

    pub fn inner_first_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        self.inner.inner_first_non_blank_column(inner_identifier)
    }

    pub fn inner_last_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        self.inner.inner_last_non_blank_column(inner_identifier)
    }

    /// Returns the index of the first column that is related to a TSM.
    pub fn first_interesting_column(&self) -> SourceColumn {
        self.inners()
            .iter_values()
            .map(|inner| {
                [
                    self.reference_arrangement_char_to_source_column(
                        inner.template_switch().sp1_reference,
                    ),
                    self.query_arrangement_char_to_source_column(inner.template_switch().sp1_query),
                    inner.template_switch().sp2_secondary,
                    inner.template_switch().sp3_secondary,
                    self.reference_arrangement_char_to_source_column(
                        inner.template_switch().sp4_reference,
                    ),
                    self.query_arrangement_char_to_source_column(inner.template_switch().sp4_query),
                    inner
                        .sequence()
                        .iter_values()
                        .filter_map(|c| {
                            if c.is_gap_or_blank() {
                                None
                            } else {
                                Some(c.source_column())
                            }
                        })
                        .next()
                        .unwrap(),
                ]
                .into_iter()
                .min()
                .unwrap()
            })
            .min()
            .unwrap_or(0.into())
    }

    /// Returns the index of the last column that is related to a TSM.
    pub fn last_interesting_column(&self) -> SourceColumn {
        self.inners()
            .iter_values()
            .map(|inner| {
                [
                    self.reference_arrangement_char_to_source_column(
                        inner.template_switch().sp1_reference,
                    )
                    .saturating_sub(1),
                    self.query_arrangement_char_to_source_column(inner.template_switch().sp1_query)
                        .saturating_sub(1),
                    inner.template_switch().sp2_secondary.saturating_sub(1),
                    inner.template_switch().sp3_secondary.saturating_sub(1),
                    self.reference_arrangement_char_to_source_column(
                        inner.template_switch().sp4_reference,
                    )
                    .saturating_sub(1),
                    self.query_arrangement_char_to_source_column(inner.template_switch().sp4_query)
                        .saturating_sub(1),
                    inner
                        .sequence()
                        .iter_values()
                        .rev()
                        .filter_map(|c| {
                            if c.is_gap_or_blank() {
                                None
                            } else {
                                Some(c.source_column())
                            }
                        })
                        .next()
                        .unwrap(),
                ]
                .into_iter()
                .max()
                .unwrap()
            })
            .max()
            .unwrap_or(
                self.source
                    .reference_length()
                    .max(self.source.query_length())
                    .into(),
            )
    }
}
