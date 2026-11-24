use std::fmt::Display;

use crate::alignment_result::AlignmentType;
use iter::{
    CompactAlignmentIter, CompactAlignmentIterCloned, FlatAlignmentIter, FlatAlignmentIterCloned,
};

pub mod iter;
pub mod stream;
pub mod template_switch_specifics;

#[derive(Debug, Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Alignment {
    alignment: Vec<(usize, AlignmentType)>,
}

impl Alignment {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn push(&mut self, alignment_type: AlignmentType) {
        self.push_n(1, alignment_type);
    }

    pub fn push_n(&mut self, multiplicity: usize, alignment_type: AlignmentType) {
        if let Some((existing_multiplicity, last_alignment_type)) = self.alignment.last_mut() {
            if *last_alignment_type == alignment_type {
                *existing_multiplicity += multiplicity;
            } else {
                self.alignment.push((multiplicity, alignment_type));
            }
        } else {
            self.alignment.push((multiplicity, alignment_type));
        }
    }

    pub fn inner_mut(&mut self) -> &mut Vec<(usize, AlignmentType)> {
        &mut self.alignment
    }

    pub fn into_inner(self) -> Vec<(usize, AlignmentType)> {
        self.alignment
    }
}

impl Alignment {
    pub fn iter_compact(&self) -> CompactAlignmentIter<'_> {
        CompactAlignmentIter::new(&self.alignment)
    }

    pub fn iter_compact_cloned(&self) -> CompactAlignmentIterCloned<'_> {
        CompactAlignmentIterCloned::new(&self.alignment)
    }

    pub fn iter_flat(&self) -> FlatAlignmentIter<'_> {
        FlatAlignmentIter::new(&self.alignment)
    }

    pub fn iter_flat_cloned(&self) -> FlatAlignmentIterCloned<'_> {
        FlatAlignmentIterCloned::new(&self.alignment)
    }

    pub fn cigar(&self) -> String {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl std::fmt::Write) -> std::fmt::Result {
        for (amount, alignment_type) in &self.alignment {
            if alignment_type.is_repeatable() {
                write!(writer, "{amount}{alignment_type}")?;
            } else {
                write!(writer, "{alignment_type}")?;
            }
        }

        Ok(())
    }

    pub fn reverse(&self) -> Self {
        Self {
            alignment: self.alignment.iter().cloned().rev().collect(),
        }
    }
}

impl From<Vec<(usize, AlignmentType)>> for Alignment {
    fn from(value: Vec<(usize, AlignmentType)>) -> Self {
        Self { alignment: value }
    }
}

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.write_cigar(f)
    }
}
