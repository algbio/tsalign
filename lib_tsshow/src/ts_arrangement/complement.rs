use lib_tsalign::a_star_aligner::template_switch_distance::TemplateSwitchAncestor;
use tagged_vec::TaggedVec;

use super::{
    character::Char,
    index_types::{ArrangementColumn, SourceColumn},
    source::{SourceChar, TsSourceArrangement},
};

pub struct TsComplementArrangement {
    reference_c: TaggedVec<ArrangementColumn, ComplementChar>,
    query_c: TaggedVec<ArrangementColumn, ComplementChar>,
}

#[derive(Debug, Clone, Copy)]
pub enum ComplementChar {
    Complement {
        column: SourceColumn,
        lower_case: bool,
        source_hidden: bool,
    },
    Hidden {
        column: SourceColumn,
        source_hidden: bool,
    },
    Gap {
        source_hidden: bool,
    },
    Blank,
}

impl TsComplementArrangement {
    pub fn new(source_arrangement: &TsSourceArrangement) -> Self {
        let mut result = Self {
            reference_c: Default::default(),
            query_c: Default::default(),
        };

        for (sequence, sequence_c) in [
            (source_arrangement.reference(), &mut result.reference_c),
            (source_arrangement.query(), &mut result.query_c),
        ] {
            for source_character in sequence.iter_values() {
                match source_character {
                    SourceChar::Source {
                        column,
                        copy_depth: None,
                        ..
                    } => {
                        sequence_c.push(ComplementChar::new_hidden(*column, false));
                    }
                    SourceChar::Hidden {
                        column,
                        copy_depth: None,
                    } => {
                        sequence_c.push(ComplementChar::new_hidden(*column, true));
                    }
                    SourceChar::Source {
                        copy_depth: Some(_),
                        ..
                    }
                    | SourceChar::Hidden {
                        copy_depth: Some(_),
                        ..
                    }
                    | SourceChar::Gap { .. }
                    | SourceChar::Separator
                    | SourceChar::Spacer
                    | SourceChar::Blank => {
                        sequence_c.push(ComplementChar::Blank);
                    }
                }
            }
        }

        result
    }

    pub fn width(&self) -> usize {
        debug_assert_eq!(self.reference_c.len(), self.query_c.len());
        self.reference_c.len()
    }

    pub fn ancestor_complement(
        &self,
        ancestor: TemplateSwitchAncestor,
    ) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        match ancestor {
            TemplateSwitchAncestor::Reference => self.reference_complement(),
            TemplateSwitchAncestor::Query => self.query_complement(),
        }
    }

    pub fn reference_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        &self.reference_c
    }

    pub fn query_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        &self.query_c
    }

    pub fn reference_complement_mut(
        &mut self,
    ) -> &mut TaggedVec<ArrangementColumn, ComplementChar> {
        &mut self.reference_c
    }

    pub fn query_complement_mut(&mut self) -> &mut TaggedVec<ArrangementColumn, ComplementChar> {
        &mut self.query_c
    }

    pub fn show_ancestor_character(
        &mut self,
        ancestor: TemplateSwitchAncestor,
        column: ArrangementColumn,
    ) {
        match ancestor {
            TemplateSwitchAncestor::Reference => self.show_reference_character(column),
            TemplateSwitchAncestor::Query => self.show_query_character(column),
        }
    }

    pub fn show_reference_character(&mut self, column: ArrangementColumn) {
        self.reference_c[column].show();
    }

    pub fn show_query_character(&mut self, column: ArrangementColumn) {
        self.query_c[column].show();
    }

    pub fn ancestor_to_lower_case(
        &mut self,
        ancestor: TemplateSwitchAncestor,
        column: ArrangementColumn,
    ) {
        match ancestor {
            TemplateSwitchAncestor::Reference => self.reference_complement_to_lower_case(column),
            TemplateSwitchAncestor::Query => self.query_complement_to_lower_case(column),
        }
    }

    pub fn reference_complement_to_lower_case(&mut self, column: ArrangementColumn) {
        self.reference_c[column].to_lower_case();
    }

    pub fn query_complement_to_lower_case(&mut self, column: ArrangementColumn) {
        self.query_c[column].to_lower_case();
    }

    pub fn insert_ancestor_complement_gap(
        &mut self,
        ancestor: TemplateSwitchAncestor,
        column: ArrangementColumn,
    ) {
        match ancestor {
            TemplateSwitchAncestor::Reference => self.insert_reference_complement_gap(column),
            TemplateSwitchAncestor::Query => self.insert_query_complement_gap(column),
        }
    }

    pub fn insert_reference_complement_gap(&mut self, column: ArrangementColumn) {
        let source_hidden = Self::is_insert_gap_source_hidden(&self.reference_c, column);
        self.reference_c
            .insert(column, ComplementChar::Gap { source_hidden });
        self.query_c.insert(column, ComplementChar::Blank);
    }

    pub fn insert_query_complement_gap(&mut self, column: ArrangementColumn) {
        let source_hidden = Self::is_insert_gap_source_hidden(&self.query_c, column);
        self.reference_c.insert(column, ComplementChar::Blank);
        self.query_c
            .insert(column, ComplementChar::Gap { source_hidden });
    }

    pub fn insert_blank(&mut self, column: ArrangementColumn) {
        self.reference_c.insert(column, ComplementChar::Blank);
        self.query_c.insert(column, ComplementChar::Blank);
    }

    fn is_insert_gap_source_hidden(
        sequence: &TaggedVec<ArrangementColumn, ComplementChar>,
        column: ArrangementColumn,
    ) -> bool {
        sequence
            .iter_values()
            .skip(column.primitive())
            .filter_map(ComplementChar::source_hidden_option)
            .next()
            .unwrap_or(true)
            && sequence
                .iter_values()
                .take(column.primitive())
                .rev()
                .filter_map(ComplementChar::source_hidden_option)
                .next()
                .unwrap_or(true)
    }

    pub fn remove_columns(&mut self, columns: impl IntoIterator<Item = ArrangementColumn> + Clone) {
        self.reference_c.remove_multi(columns.clone());
        self.query_c.remove_multi(columns);
    }
}

impl ComplementChar {
    pub fn new_hidden(column: SourceColumn, source_hidden: bool) -> Self {
        Self::Hidden {
            column,
            source_hidden,
        }
    }

    pub fn source_hidden_option(&self) -> Option<bool> {
        match self {
            ComplementChar::Complement { source_hidden, .. }
            | ComplementChar::Hidden { source_hidden, .. }
            | ComplementChar::Gap { source_hidden } => Some(*source_hidden),
            ComplementChar::Blank => None,
        }
    }

    pub fn show(&mut self) {
        match self {
            ComplementChar::Complement { .. } => { /* Do nothing */ }
            ComplementChar::Hidden {
                column,
                source_hidden,
            } => {
                *self = ComplementChar::Complement {
                    column: *column,
                    lower_case: false,
                    source_hidden: *source_hidden,
                };
            }
            ComplementChar::Gap { .. } | ComplementChar::Blank => panic!("Not showable"),
        }
    }

    pub fn to_lower_case(&mut self) {
        match self {
            ComplementChar::Complement { lower_case, .. } => *lower_case = true,
            ComplementChar::Hidden { .. } | ComplementChar::Gap { .. } | ComplementChar::Blank => {
                panic!("Not lowercasable")
            }
        }
    }
}

impl Char for ComplementChar {
    fn source_column(&self) -> SourceColumn {
        match self {
            Self::Complement { column, .. } | Self::Hidden { column, .. } => *column,
            Self::Gap { .. } | Self::Blank => panic!("Not a char"),
        }
    }

    fn is_char(&self) -> bool {
        matches!(self, Self::Complement { .. } | Self::Hidden { .. })
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
        matches!(self, Self::Hidden { .. })
    }
}
