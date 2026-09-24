#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlignmentHistoryOperation {
    Match,
    Substitution,
    GapA,
    GapB,
}

impl AlignmentHistoryOperation {
    pub fn to_bits(&self) -> usize {
        match self {
            Self::Match => 0b00,
            Self::Substitution => 0b01,
            Self::GapA => 0b10,
            Self::GapB => 0b11,
        }
    }

    pub fn from_bits(bits: usize) -> Option<Self> {
        match bits {
            0b00 => Some(Self::Match),
            0b01 => Some(Self::Substitution),
            0b10 => Some(Self::GapA),
            0b11 => Some(Self::GapB),
            _ => None,
        }
    }
}
