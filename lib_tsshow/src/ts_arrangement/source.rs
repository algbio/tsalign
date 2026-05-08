use std::{cmp::Ordering, iter};

use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{
        AlignmentType, TemplateSwitchDirection, TemplateSwitchPrimary, TemplateSwitchSecondary,
    },
};
use log::{debug, trace};
use tagged_vec::TaggedVec;

use super::{
    character::Char,
    index_types::{ArrangementCharColumn, ArrangementColumn, SourceColumn},
    template_switch::TemplateSwitch,
};
use crate::error::Result;

pub struct TsSourceArrangement {
    reference: TaggedVec<ArrangementColumn, SourceChar>,
    query: TaggedVec<ArrangementColumn, SourceChar>,
    reference_length: usize,
    query_length: usize,
}

pub struct RemovedHiddenChars {
    reference: Vec<ArrangementCharColumn>,
    query: Vec<ArrangementCharColumn>,
}

#[derive(Debug, Clone, Copy)]
pub enum SourceChar {
    Source {
        column: SourceColumn,
        lower_case: bool,
        copy_depth: Option<usize>,
    },
    Hidden {
        column: SourceColumn,
        copy_depth: Option<usize>,
    },
    Gap {
        copy_depth: Option<usize>,
    },
    /// Separates the aligned section of the sequences from the unaligned section.
    ///
    /// This is used when the alignment is only computed on a substring of the reference and the query.
    Separator,
    /// Like a blank, but not treated as an empty column.
    Spacer,
    Blank,
}

impl TsSourceArrangement {
    pub fn new(
        reference_alignment_offset: usize,
        query_alignment_offset: usize,
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
        template_switches_out: &mut impl Extend<TemplateSwitch>,
    ) -> Result<Self> {
        let mut ts_index = 0;
        let mut alignment = alignment.into_iter();

        // Insert the source sequences with blanks to the left to make sure their alignments start at the same arrangement column.
        let reference_left_blank_count =
            query_alignment_offset.saturating_sub(reference_alignment_offset);
        let query_left_blank_count =
            reference_alignment_offset.saturating_sub(query_alignment_offset);
        let mut result = Self {
            reference: FromIterator::from_iter(
                iter::repeat_n(SourceChar::Blank, reference_left_blank_count)
                    .chain((0..reference_length).map(SourceChar::new_source)),
            ),
            query: FromIterator::from_iter(
                iter::repeat_n(SourceChar::Blank, query_left_blank_count)
                    .chain((0..query_length).map(SourceChar::new_source)),
            ),
            reference_length,
            query_length,
        };

        let mut current_reference_index =
            ArrangementColumn::from(reference_left_blank_count + reference_alignment_offset);
        let mut current_query_index =
            ArrangementColumn::from(query_left_blank_count + query_alignment_offset);
        debug_assert_eq!(current_reference_index, current_query_index);

        // If there are further characters to the left of the alignment, we insert a separator into both sequences.
        if reference_alignment_offset > 0 || query_alignment_offset > 0 {
            result
                .reference
                .insert(current_reference_index, SourceChar::Separator);
            current_reference_index += 1;

            result
                .query
                .insert(current_query_index, SourceChar::Separator);
            current_query_index += 1;
        }
        debug_assert_eq!(current_reference_index, current_query_index);

        while let Some(alignment_type) = alignment.next() {
            trace!(
                "Source alignment type: {alignment_type}; R/Q indices: {current_reference_index}/{current_query_index} {}/{}",
                result.reference_arrangement_to_source_column(current_reference_index),
                result.query_arrangement_to_source_column(current_query_index),
            );

            match alignment_type {
                AlignmentType::PrimaryInsertion | AlignmentType::PrimaryFlankInsertion => {
                    result.reference.insert(
                        current_reference_index,
                        SourceChar::Gap {
                            copy_depth: result.query[current_query_index].copy_depth(),
                        },
                    );
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => {
                    result.query.insert(
                        current_query_index,
                        SourceChar::Gap {
                            copy_depth: result.reference[current_reference_index].copy_depth(),
                        },
                    );
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimarySubstitution | AlignmentType::PrimaryFlankSubstitution => {
                    result.reference[current_reference_index].to_lower_case();
                    result.query[current_query_index].to_lower_case();
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimaryMatch | AlignmentType::PrimaryFlankMatch => {
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::TemplateSwitchEntrance {
                    primary,
                    secondary,
                    direction,
                    first_offset,
                    ..
                } => {
                    template_switches_out.extend([result.align_ts(
                        ts_index,
                        primary,
                        secondary,
                        direction,
                        first_offset,
                        &mut alignment,
                        &mut current_reference_index,
                        &mut current_query_index,
                    )]);
                    ts_index += 1;
                }

                AlignmentType::Root | AlignmentType::PrimaryReentry => { /* Do nothing */ }
                AlignmentType::TemplateSwitchExit { .. }
                | AlignmentType::SecondaryInsertion
                | AlignmentType::SecondaryDeletion
                | AlignmentType::SecondarySubstitution
                | AlignmentType::SecondaryMatch
                | AlignmentType::SecondaryRoot
                | AlignmentType::PrimaryShortcut { .. } => unreachable!(),
            }
        }

        debug_assert_eq!(current_reference_index, current_query_index);

        // If there are further source chars right of the alignment, we insert a separator into both sequences.
        if result
            .try_reference_arrangement_to_source_column(current_reference_index)
            .map(|source_current_reference_index| {
                source_current_reference_index < SourceColumn::new(reference_length - 1)
            })
            .unwrap_or(false)
            || result
                .try_query_arrangement_to_source_column(current_query_index)
                .map(|source_current_query_index| {
                    source_current_query_index < SourceColumn::new(query_length - 1)
                })
                .unwrap_or(false)
        {
            result
                .reference
                .insert(current_reference_index, SourceChar::Separator);
            current_reference_index += 1;

            result
                .query
                .insert(current_query_index, SourceChar::Separator);
            current_query_index += 1;
        }

        debug_assert_eq!(current_reference_index, current_query_index);

        // Fill the back with blanks to ensure equal length.
        while result.reference.len() < result.query.len() {
            result.reference.push(SourceChar::Blank);
        }
        while result.query.len() < result.reference.len() {
            result.query.push(SourceChar::Blank);
        }

        debug_assert_eq!(result.reference.len(), result.query.len());
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn align_ts(
        &mut self,
        ts_index: usize,
        ts_primary: TemplateSwitchPrimary,
        ts_secondary: TemplateSwitchSecondary,
        ts_direction: TemplateSwitchDirection,
        first_offset: isize,
        mut alignment: impl Iterator<Item = AlignmentType>,
        current_reference_index: &mut ArrangementColumn,
        current_query_index: &mut ArrangementColumn,
    ) -> TemplateSwitch {
        debug!("Aligning template switch #{ts_index}");
        let sp1_reference =
            self.reference_arrangement_to_arrangement_char_column(*current_reference_index);
        let sp1_query = self.query_arrangement_to_arrangement_char_column(*current_query_index);
        let sp2_secondary = match ts_secondary {
            TemplateSwitchSecondary::Reference => {
                let source_current_reference_index =
                    self.reference_arrangement_to_source_column(*current_reference_index);
                let adjusted_source_current_reference_index =
                    source_current_reference_index
                        .checked_sub(self.count_reference_copy_chars_before_next_real_char(
                            *current_reference_index,
                        ))
                        .unwrap();
                trace!(
                    "current_reference_index: {current_reference_index} -> {source_current_reference_index} -> {adjusted_source_current_reference_index}"
                );
                adjusted_source_current_reference_index
            }
            TemplateSwitchSecondary::Query => {
                let source_current_query_index =
                    self.query_arrangement_to_source_column(*current_query_index);
                let adjusted_source_current_query_index = source_current_query_index
                    .checked_sub(
                        self.count_query_copy_chars_before_next_real_char(*current_query_index),
                    )
                    .unwrap();
                trace!(
                    "current_query_index: {current_query_index} -> {source_current_query_index} -> {adjusted_source_current_query_index}"
                );
                adjusted_source_current_query_index
            }
        } + first_offset;
        trace!("first_offset: {first_offset}");
        trace!("sp2_secondary: {sp2_secondary}");

        let mut sp3_secondary = sp2_secondary;
        let mut primary_inner_length = 0;
        let mut inner_alignment = Alignment::new();

        let anti_primary_gap = loop {
            let alignment_type = alignment.next().unwrap_or_else(|| unreachable!());
            trace!("Secondary alignment type: {alignment_type}");

            match alignment_type {
                AlignmentType::TemplateSwitchExit { anti_primary_gap } => {
                    break anti_primary_gap;
                }
                AlignmentType::SecondaryDeletion => {
                    match ts_direction {
                        TemplateSwitchDirection::Forward => sp3_secondary += 1,
                        TemplateSwitchDirection::Reverse => sp3_secondary -= 1,
                    }
                    inner_alignment.push(alignment_type);
                }

                AlignmentType::SecondarySubstitution | AlignmentType::SecondaryMatch => {
                    match ts_direction {
                        TemplateSwitchDirection::Forward => sp3_secondary += 1,
                        TemplateSwitchDirection::Reverse => sp3_secondary -= 1,
                    }
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
                }
                AlignmentType::SecondaryInsertion => {
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
                }
                AlignmentType::SecondaryRoot => { /* Do nothing */ }
                _ => unreachable!(),
            }
        };
        let sp3_secondary = sp3_secondary;
        let primary_inner_length = primary_inner_length;

        let (primary, anti_primary, current_primary_index, current_anti_primary_index) =
            match ts_primary {
                TemplateSwitchPrimary::Reference => (
                    &mut self.reference,
                    &mut self.query,
                    current_reference_index,
                    current_query_index,
                ),
                TemplateSwitchPrimary::Query => (
                    &mut self.query,
                    &mut self.reference,
                    current_query_index,
                    current_reference_index,
                ),
            };

        // Mark the TS inner in the source arrangement as hidden.
        // Take a copy of the unhidden characters at the same time.
        let inner = primary
            .iter_values_mut()
            .skip(current_primary_index.primitive())
            .take(primary_inner_length)
            .map(|c| {
                let result = *c;
                c.hide();
                result
            })
            .collect();

        *current_primary_index += primary_inner_length;

        let anti_primary_inner_length = if anti_primary_gap < 0 {
            // Insert copies.
            let duplicate_rev: Vec<_> = anti_primary
                .iter_values()
                .take(current_anti_primary_index.primitive())
                .rev()
                .filter_map(|c| {
                    if c.is_char() {
                        Some(c.make_visible_copy())
                    } else {
                        None
                    }
                })
                .take(usize::try_from(-anti_primary_gap).unwrap())
                .collect();

            anti_primary.splice(
                *current_anti_primary_index..*current_anti_primary_index,
                duplicate_rev.into_iter().rev(),
            );
            0
        } else {
            *current_anti_primary_index += usize::try_from(anti_primary_gap).unwrap();
            usize::try_from(anti_primary_gap).unwrap()
        };

        // Insert spacers.
        let mut required_spacer_count = 4usize.saturating_sub(anti_primary_inner_length);

        match primary_inner_length.cmp(&anti_primary_inner_length) {
            Ordering::Less => {
                let delta = anti_primary_inner_length
                    .checked_sub(primary_inner_length)
                    .unwrap();
                primary.splice(
                    *current_primary_index..*current_primary_index,
                    iter::repeat_n(SourceChar::Blank, delta),
                );
                *current_primary_index += delta;
            }
            Ordering::Equal => { /* Do nothing */ }
            Ordering::Greater => {
                let delta = primary_inner_length
                    .checked_sub(anti_primary_inner_length)
                    .unwrap();
                anti_primary.splice(
                    *current_anti_primary_index..*current_anti_primary_index,
                    iter::repeat_n(SourceChar::Spacer, required_spacer_count)
                        .chain(iter::repeat(SourceChar::Blank))
                        .take(delta),
                );
                required_spacer_count = required_spacer_count.saturating_sub(delta);
                *current_anti_primary_index += delta;
            }
        }

        primary.splice(
            *current_primary_index..*current_primary_index,
            iter::repeat_n(SourceChar::Blank, required_spacer_count),
        );
        anti_primary.splice(
            *current_anti_primary_index..*current_anti_primary_index,
            iter::repeat_n(SourceChar::Spacer, required_spacer_count),
        );
        *current_primary_index += required_spacer_count;
        *current_anti_primary_index += required_spacer_count;

        let (current_reference_index, current_query_index) = match ts_primary {
            TemplateSwitchPrimary::Reference => (current_primary_index, current_anti_primary_index),
            TemplateSwitchPrimary::Query => (current_anti_primary_index, current_primary_index),
        };

        let sp4_reference =
            self.reference_arrangement_to_arrangement_char_column(*current_reference_index);
        let sp4_query = self.query_arrangement_to_arrangement_char_column(*current_query_index);

        TemplateSwitch {
            index: ts_index,
            primary: ts_primary,
            secondary: ts_secondary,
            sp1_reference,
            sp1_query,
            sp4_reference,
            sp4_query,
            sp2_secondary,
            sp3_secondary,
            inner,
            inner_alignment,
        }
    }

    pub fn secondary(
        &self,
        secondary: TemplateSwitchSecondary,
    ) -> &TaggedVec<ArrangementColumn, SourceChar> {
        match secondary {
            TemplateSwitchSecondary::Reference => self.reference(),
            TemplateSwitchSecondary::Query => self.query(),
        }
    }

    pub fn reference(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        &self.reference
    }

    pub fn query(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        &self.query
    }

    pub fn width(&self) -> usize {
        debug_assert_eq!(self.reference.len(), self.query.len());
        self.reference.len()
    }

    pub fn secondary_to_lower_case(
        &mut self,
        secondary: TemplateSwitchSecondary,
        column: ArrangementColumn,
    ) {
        match secondary {
            TemplateSwitchSecondary::Reference => self.reference[column].to_lower_case(),
            TemplateSwitchSecondary::Query => self.query[column].to_lower_case(),
        }
    }

    pub fn insert_secondary_gap_with_minimum_copy_depth(
        &mut self,
        secondary: TemplateSwitchSecondary,
        column: ArrangementColumn,
    ) {
        let secondary_sequence = self.secondary(secondary);
        let copy_depth = if column == ArrangementColumn::ZERO {
            secondary_sequence[column].copy_depth()
        } else if column == ArrangementColumn::from(secondary_sequence.len()) {
            secondary_sequence[column - 1].copy_depth()
        } else {
            let copy_depth_1 = secondary_sequence[column - 1].copy_depth();
            let copy_depth_2 = secondary_sequence[column].copy_depth();

            if let (Some(copy_depth_1), Some(copy_depth_2)) = (copy_depth_1, copy_depth_2) {
                Some(copy_depth_1.min(copy_depth_2))
            } else {
                None
            }
        };

        self.insert_secondary_gap(secondary, column, copy_depth);
    }

    pub fn insert_secondary_gap(
        &mut self,
        secondary: TemplateSwitchSecondary,
        column: ArrangementColumn,
        copy_depth: Option<usize>,
    ) {
        match secondary {
            TemplateSwitchSecondary::Reference => self.insert_reference_gap(column, copy_depth),
            TemplateSwitchSecondary::Query => self.insert_query_gap(column, copy_depth),
        }
    }

    pub fn insert_reference_gap(&mut self, column: ArrangementColumn, copy_depth: Option<usize>) {
        self.reference
            .insert(column, SourceChar::Gap { copy_depth });
        self.query.insert(column, SourceChar::Blank);
    }

    pub fn insert_query_gap(&mut self, column: ArrangementColumn, copy_depth: Option<usize>) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Gap { copy_depth });
    }

    pub fn insert_blank(&mut self, column: ArrangementColumn) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Blank);
    }

    pub fn remove_columns(
        &mut self,
        columns: impl IntoIterator<Item = ArrangementColumn> + Clone,
    ) -> RemovedHiddenChars {
        let result = RemovedHiddenChars {
            reference: columns
                .clone()
                .into_iter()
                .filter_map(|c| {
                    if self.reference[c].is_char() {
                        Some(self.reference_arrangement_to_arrangement_char_column(c))
                    } else {
                        None
                    }
                })
                .collect(),
            query: columns
                .clone()
                .into_iter()
                .filter_map(|c| {
                    if self.query[c].is_char() {
                        Some(self.query_arrangement_to_arrangement_char_column(c))
                    } else {
                        None
                    }
                })
                .collect(),
        };

        trace!(
            "Removing reference columns: {:?}",
            columns
                .clone()
                .into_iter()
                .map(|c| self.reference[c])
                .collect::<Vec<_>>()
        );
        trace!(
            "Removing query columns: {:?}",
            columns
                .clone()
                .into_iter()
                .map(|c| self.query[c])
                .collect::<Vec<_>>()
        );

        self.reference.remove_multi(columns.clone());
        self.query.remove_multi(columns);
        result
    }

    pub fn reference_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> ArrangementColumn {
        Self::source_to_arrangement_column(&self.reference, column)
    }

    pub fn try_reference_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> Option<ArrangementColumn> {
        Self::try_source_to_arrangement_column(&self.reference, column)
    }

    pub fn query_source_to_arrangement_column(&self, column: SourceColumn) -> ArrangementColumn {
        Self::source_to_arrangement_column(&self.query, column)
    }

    pub fn try_query_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> Option<ArrangementColumn> {
        Self::try_source_to_arrangement_column(&self.query, column)
    }

    fn source_to_arrangement_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        source_column: SourceColumn,
    ) -> ArrangementColumn {
        Self::try_source_to_arrangement_column(sequence, source_column)
            .unwrap_or_else(|| panic!("Source column {source_column} has no matching arrangement column. There are {} source columns in the arrangement.", sequence.iter_values().filter(|c| matches!(c, SourceChar::Source {  .. } | SourceChar::Hidden {  .. })).count()))
    }

    fn try_source_to_arrangement_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        source_column: SourceColumn,
    ) -> Option<ArrangementColumn> {
        sequence
            .iter()
            .filter_map(|(i, c)| match c {
                SourceChar::Source { column, .. } | SourceChar::Hidden { column, .. }
                    if *column + 1usize == source_column =>
                {
                    Some(sequence.len().into())
                }
                SourceChar::Source { column, .. } | SourceChar::Hidden { column, .. }
                    if *column == source_column =>
                {
                    Some(i)
                }
                _ => None,
            })
            .min()
    }

    pub fn reference_arrangement_to_arrangement_char_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        Self::arrangement_to_arrangement_char_column(&self.reference, arrangement_column)
    }

    pub fn query_arrangement_to_arrangement_char_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        Self::arrangement_to_arrangement_char_column(&self.query, arrangement_column)
    }

    fn arrangement_to_arrangement_char_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        assert!(sequence[arrangement_column].is_char());
        sequence
            .iter_values()
            .take(arrangement_column.primitive())
            .filter(|c| c.is_char())
            .count()
            .into()
    }

    pub fn reference_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        Self::arrangement_to_source_column(&self.reference, arrangement_column)
    }

    pub fn try_reference_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> Option<SourceColumn> {
        Self::try_arrangement_to_source_column(&self.reference, arrangement_column)
    }

    pub fn query_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        Self::arrangement_to_source_column(&self.query, arrangement_column)
    }

    pub fn try_query_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> Option<SourceColumn> {
        Self::try_arrangement_to_source_column(&self.query, arrangement_column)
    }

    fn arrangement_to_source_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        Self::try_arrangement_to_source_column(sequence, arrangement_column).unwrap()
    }

    fn try_arrangement_to_source_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        arrangement_column: ArrangementColumn,
    ) -> Option<SourceColumn> {
        if arrangement_column.primitive() >= sequence.len() {
            return None;
        }

        // This may also be called on a non-source char.
        assert!(sequence[arrangement_column].is_char());
        Some(
            sequence
                .iter_values()
                .take(arrangement_column.into())
                .filter(|c| c.is_source_char())
                .count()
                .into(),
        )
    }

    pub fn reference_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        Self::arrangement_char_to_arrangement_column(&self.reference, column)
    }

    pub fn query_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        Self::arrangement_char_to_arrangement_column(&self.query, column)
    }

    fn arrangement_char_to_arrangement_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        sequence
            .iter()
            .filter_map(|(i, c)| if c.is_char() { Some(i) } else { None })
            .nth(column.primitive())
            .unwrap()
    }

    pub fn reference_arrangement_char_to_source_column(
        &self,
        column: ArrangementCharColumn,
    ) -> SourceColumn {
        Self::arrangement_char_to_source_column(&self.reference, column)
    }

    pub fn query_arrangement_char_to_source_column(
        &self,
        column: ArrangementCharColumn,
    ) -> SourceColumn {
        Self::arrangement_char_to_source_column(&self.query, column)
    }

    fn arrangement_char_to_source_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        column: ArrangementCharColumn,
    ) -> SourceColumn {
        sequence
            .iter_values()
            .filter_map(|c| {
                if c.is_char() {
                    Some(c.source_column())
                } else {
                    None
                }
            })
            .nth(column.primitive())
            .unwrap()
    }

    fn count_reference_copy_chars_before_next_real_char(&self, offset: ArrangementColumn) -> usize {
        Self::count_copy_chars_before_next_real_char(&self.reference, offset)
    }

    fn count_query_copy_chars_before_next_real_char(&self, offset: ArrangementColumn) -> usize {
        Self::count_copy_chars_before_next_real_char(&self.query, offset)
    }

    fn count_copy_chars_before_next_real_char(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        offset: ArrangementColumn,
    ) -> usize {
        sequence
            .iter_values()
            .skip(offset.into())
            .take_while(|c| !c.is_source_char())
            .filter(|c| c.is_char() && c.is_copy())
            .count()
    }

    pub fn reference_length(&self) -> usize {
        self.reference_length
    }

    pub fn query_length(&self) -> usize {
        self.query_length
    }
}

impl SourceChar {
    pub fn new_source(column: impl Into<SourceColumn>) -> Self {
        Self::Source {
            column: column.into(),
            lower_case: false,
            copy_depth: None,
        }
    }

    pub fn is_copy(&self) -> bool {
        match self {
            Self::Source { copy_depth, .. }
            | Self::Hidden { copy_depth, .. }
            | Self::Gap { copy_depth } => copy_depth.is_some(),
            Self::Separator | Self::Spacer | Self::Blank => panic!("Blank has no copy property"),
        }
    }

    pub fn copy_depth(&self) -> Option<usize> {
        match self {
            Self::Source { copy_depth, .. } | Self::Hidden { copy_depth, .. } => *copy_depth,
            Self::Gap { copy_depth } => *copy_depth,
            Self::Separator | Self::Spacer | Self::Blank => panic!("Blank has no copy property"),
        }
    }

    pub fn to_lower_case(&mut self) {
        match self {
            Self::Source { lower_case, .. } => *lower_case = true,
            Self::Hidden { .. }
            | Self::Gap { .. }
            | Self::Separator
            | Self::Spacer
            | Self::Blank => {
                panic!("Not lowercasable: {self:?}")
            }
        }
    }

    pub fn to_upper_case(&mut self) {
        match self {
            Self::Source { lower_case, .. } => *lower_case = false,
            Self::Hidden { .. }
            | Self::Gap { .. }
            | Self::Separator
            | Self::Spacer
            | Self::Blank => {
                panic!("Not uppercasable: {self:?}")
            }
        }
    }

    pub fn hide(&mut self) {
        match self {
            Self::Source {
                column,
                lower_case,
                copy_depth,
            } => {
                assert!(!*lower_case);
                *self = Self::Hidden {
                    column: *column,
                    copy_depth: *copy_depth,
                };
            }
            Self::Hidden { .. } => unreachable!("Already hidden"),
            Self::Gap { .. } | Self::Separator | Self::Spacer | Self::Blank => {
                unreachable!("Cannot be hidden: {self:?}")
            }
        }
    }

    pub fn make_copy(&self) -> Self {
        match self {
            Self::Source {
                column, copy_depth, ..
            } => Self::Source {
                column: *column,
                lower_case: false,
                copy_depth: Some(copy_depth.map(|copy_depth| copy_depth + 1).unwrap_or(0)),
            },
            Self::Hidden { column, copy_depth } => Self::Hidden {
                column: *column,
                copy_depth: Some(copy_depth.map(|copy_depth| copy_depth + 1).unwrap_or(0)),
            },
            Self::Gap { .. } | Self::Separator | Self::Spacer | Self::Blank => {
                panic!("Should never be copied: {self:?}")
            }
        }
    }

    pub fn make_visible_copy(&self) -> Self {
        match self {
            Self::Source {
                column, copy_depth, ..
            }
            | Self::Hidden { column, copy_depth } => Self::Source {
                column: *column,
                lower_case: false,
                copy_depth: Some(copy_depth.map(|copy_depth| copy_depth + 1).unwrap_or(0)),
            },
            Self::Gap { .. } | Self::Separator | Self::Spacer | Self::Blank => {
                panic!("Should never be copied: {self:?}")
            }
        }
    }
}

impl Char for SourceChar {
    fn source_column(&self) -> SourceColumn {
        match self {
            Self::Source { column, .. } | Self::Hidden { column, .. } => *column,
            Self::Gap { .. } | Self::Separator | Self::Spacer | Self::Blank => panic!("Not a char"),
        }
    }

    fn is_char(&self) -> bool {
        matches!(self, Self::Source { .. } | Self::Hidden { .. })
    }

    fn is_gap(&self) -> bool {
        matches!(self, Self::Gap { .. })
    }

    fn is_spacer(&self) -> bool {
        matches!(self, Self::Spacer)
    }

    fn is_blank(&self) -> bool {
        matches!(self, Self::Blank)
    }

    /// Returns true if this is a source char and not a copy of a source char.
    fn is_source_char(&self) -> bool {
        matches!(
            self,
            Self::Source {
                copy_depth: None,
                ..
            } | Self::Hidden {
                copy_depth: None,
                ..
            }
        )
    }

    fn is_hidden(&self) -> bool {
        matches!(self, Self::Hidden { .. })
    }
}

impl RemovedHiddenChars {
    pub fn reference(&self) -> &[ArrangementCharColumn] {
        &self.reference
    }

    pub fn query(&self) -> &[ArrangementCharColumn] {
        &self.query
    }
}
