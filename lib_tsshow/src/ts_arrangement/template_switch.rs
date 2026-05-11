use super::{
    index_types::{ArrangementCharColumn, SourceColumn},
    source::{RemovedHiddenChars, SourceChar},
};
use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{
        AlignmentType, EqualCostRange, TemplateSwitchPrimary, TemplateSwitchSecondary,
    },
};

#[derive(Debug, Clone)]
pub struct TemplateSwitch {
    pub index: usize,
    pub primary: TemplateSwitchPrimary,
    pub secondary: TemplateSwitchSecondary,
    pub sp1_reference: ArrangementCharColumn,
    pub sp1_query: ArrangementCharColumn,
    pub sp4_reference: ArrangementCharColumn,
    pub sp4_query: ArrangementCharColumn,
    pub sp2_secondary: SourceColumn,
    pub sp3_secondary: SourceColumn,
    pub inner: Vec<SourceChar>,
    pub inner_alignment: Alignment<AlignmentType>,
    pub equal_cost_range: EqualCostRange,
}

impl TemplateSwitch {
    pub fn remove_hidden_chars(&mut self, columns: &RemovedHiddenChars) {
        Self::remove_hidden_chars_by_sequence(
            &mut self.sp1_reference,
            &mut self.sp4_reference,
            columns.reference(),
        );
        Self::remove_hidden_chars_by_sequence(
            &mut self.sp1_query,
            &mut self.sp4_query,
            columns.query(),
        );
    }

    fn remove_hidden_chars_by_sequence(
        sp1: &mut ArrangementCharColumn,
        sp4: &mut ArrangementCharColumn,
        columns: &[ArrangementCharColumn],
    ) {
        *sp1 -= columns.iter().filter(|c| **c < *sp1).count();
        *sp4 -= columns.iter().filter(|c| **c < *sp4).count();
    }
}
