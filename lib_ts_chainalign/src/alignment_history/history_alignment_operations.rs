use std::fmt::Display;

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlignmentHistoryOperation {
    Match,
    Substitution,
    /// A gap in sequence 1, meaning that sequence 2 has characters that are missing from sequence 1.
    GapInA,
    /// A gap in sequence 2, meaning that sequence 1 has characters that are missing from sequence 2.
    GapInB,
}

impl AlignmentHistoryOperation {
    pub fn to_bits(self) -> usize {
        match self {
            Self::Match => 0b00,
            Self::Substitution => 0b01,
            Self::GapInA => 0b10,
            Self::GapInB => 0b11,
        }
    }

    pub fn from_bits(bits: usize) -> Option<Self> {
        match bits {
            0b00 => Some(Self::Match),
            0b01 => Some(Self::Substitution),
            0b10 => Some(Self::GapInA),
            0b11 => Some(Self::GapInB),
            _ => None,
        }
    }

    pub fn len_a(&self) -> usize {
        match self {
            Self::Match | Self::Substitution | Self::GapInB => 1,
            Self::GapInA => 0,
        }
    }

    pub fn len_b(&self) -> usize {
        match self {
            Self::Match | Self::Substitution | Self::GapInA => 1,
            Self::GapInB => 0,
        }
    }

    /// The number of mutations that this alignment history operation represents.
    ///
    /// A match is not a mutation, while a substitution or a gap is a mutation.
    pub fn mutations(&self) -> usize {
        match self {
            Self::Match => 0,
            Self::Substitution | Self::GapInA | Self::GapInB => 1,
        }
    }
}

impl Display for AlignmentHistoryOperation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Match => write!(f, "M"),
            Self::Substitution => write!(f, "S"),
            Self::GapInA => write!(f, "A"),
            Self::GapInB => write!(f, "B"),
        }
    }
}
