use std::fmt::Display;

use iter::{
    CompactAlignmentIter, CompactAlignmentIterCloned, FlatAlignmentIter, FlatAlignmentIterCloned,
};

use super::IAlignmentType;

pub mod iter;
pub mod stream;
pub mod template_switch_specifics;

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Alignment<AlignmentType> {
    alignment: Vec<(usize, AlignmentType)>,
}

impl<AlignmentType> Alignment<AlignmentType> {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn push(&mut self, alignment_type: AlignmentType)
    where
        AlignmentType: Eq,
    {
        self.push_n(1, alignment_type);
    }

    pub fn push_n(&mut self, multiplicity: usize, alignment_type: AlignmentType)
    where
        AlignmentType: Eq,
    {
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

impl<AlignmentType: IAlignmentType> Alignment<AlignmentType> {
    pub fn iter_compact(&self) -> CompactAlignmentIter<'_, AlignmentType> {
        CompactAlignmentIter::new(&self.alignment)
    }

    pub fn iter_compact_cloned(&self) -> CompactAlignmentIterCloned<'_, AlignmentType>
    where
        AlignmentType: Clone,
    {
        CompactAlignmentIterCloned::new(&self.alignment)
    }

    pub fn iter_flat(&self) -> FlatAlignmentIter<'_, AlignmentType> {
        FlatAlignmentIter::new(&self.alignment)
    }

    pub fn iter_flat_cloned(&self) -> FlatAlignmentIterCloned<'_, AlignmentType>
    where
        AlignmentType: Clone,
    {
        FlatAlignmentIterCloned::new(&self.alignment)
    }

    pub fn cigar(&self) -> String
    where
        AlignmentType: Display,
    {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl std::fmt::Write) -> std::fmt::Result
    where
        AlignmentType: Display,
    {
        for (amount, alignment_type) in &self.alignment {
            if alignment_type.is_repeatable() {
                write!(writer, "{amount}{alignment_type}")?;
            } else {
                write!(writer, "{alignment_type}")?;
            }
        }

        Ok(())
    }

    pub fn reverse(&self) -> Self
    where
        AlignmentType: Clone,
    {
        Self {
            alignment: self.alignment.iter().cloned().rev().collect(),
        }
    }
}

impl<AlignmentType> From<Vec<(usize, AlignmentType)>> for Alignment<AlignmentType> {
    fn from(value: Vec<(usize, AlignmentType)>) -> Self {
        Self { alignment: value }
    }
}

impl<AlignmentType: IAlignmentType + Display> Display for Alignment<AlignmentType> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.write_cigar(f)
    }
}

impl<AlignmentType> Default for Alignment<AlignmentType> {
    fn default() -> Self {
        Self {
            alignment: Default::default(),
        }
    }
}
