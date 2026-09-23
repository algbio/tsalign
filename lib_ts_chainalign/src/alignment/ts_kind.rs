use std::fmt::Display;

use lib_tsalign::a_star_aligner::template_switch_distance::{
    TemplateSwitchAncestor, TemplateSwitchDescendant,
};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct TsKind {
    pub ancestor: TsAncestor,
    pub descendant: TsDescendant,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum TsAncestor {
    Seq1,
    Seq2,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum TsDescendant {
    Seq1,
    Seq2,
}

impl TsKind {
    /// Ancestor is Seq1, descendant is Seq1.
    pub const TS11: Self = TsKind {
        ancestor: TsAncestor::Seq1,
        descendant: TsDescendant::Seq1,
    };
    /// Ancestor is Seq1, descendant is Seq2.
    pub const TS12: Self = TsKind {
        ancestor: TsAncestor::Seq1,
        descendant: TsDescendant::Seq2,
    };
    /// Ancestor is Seq2, descendant is Seq1.
    pub const TS21: Self = TsKind {
        ancestor: TsAncestor::Seq2,
        descendant: TsDescendant::Seq1,
    };
    /// Ancestor is Seq2, descendant is Seq2.
    pub const TS22: Self = TsKind {
        ancestor: TsAncestor::Seq2,
        descendant: TsDescendant::Seq2,
    };

    pub fn iter() -> impl Iterator<Item = Self> {
        [Self::TS11, Self::TS12, Self::TS21, Self::TS22].into_iter()
    }

    /// Index of this [`TsKind`] in the iterator produced by [`Self::iter()`].
    ///
    /// ```rust
    /// # use lib_ts_chainalign::alignment::ts_kind::TsKind;
    /// assert_eq!(
    ///     TsKind::iter().enumerate().collect::<Vec<_>>(),
    ///     TsKind::iter().map(|ts_kind| (ts_kind.index(), ts_kind)).collect::<Vec<_>>()
    /// );
    /// ```
    pub fn index(&self) -> usize {
        match (self.ancestor, self.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => 0,
            (TsAncestor::Seq1, TsDescendant::Seq2) => 1,
            (TsAncestor::Seq2, TsDescendant::Seq1) => 2,
            (TsAncestor::Seq2, TsDescendant::Seq2) => 3,
        }
    }

    pub fn digits(&self) -> &'static str {
        match (self.ancestor, self.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => "11",
            (TsAncestor::Seq1, TsDescendant::Seq2) => "12",
            (TsAncestor::Seq2, TsDescendant::Seq1) => "21",
            (TsAncestor::Seq2, TsDescendant::Seq2) => "22",
        }
    }

    pub fn from_index(index: usize) -> Self {
        match index {
            0 => Self::TS11,
            1 => Self::TS12,
            2 => Self::TS21,
            3 => Self::TS22,
            _ => panic!("Invalid index for TsKind: {index}"),
        }
    }

    pub fn swap(&self) -> Self {
        match (self.ancestor, self.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => Self::TS22,
            (TsAncestor::Seq1, TsDescendant::Seq2) => Self::TS21,
            (TsAncestor::Seq2, TsDescendant::Seq1) => Self::TS12,
            (TsAncestor::Seq2, TsDescendant::Seq2) => Self::TS11,
        }
    }

    pub fn swap_index(index: usize) -> usize {
        match index {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            _ => panic!("Invalid index for TsKind: {index}"),
        }
    }
}

impl TsAncestor {
    pub fn into_tsalign_secondary(self) -> TemplateSwitchAncestor {
        match self {
            TsAncestor::Seq1 => TemplateSwitchAncestor::Reference,
            TsAncestor::Seq2 => TemplateSwitchAncestor::Query,
        }
    }
}

impl TsDescendant {
    pub fn into_tsalign_primary(self) -> TemplateSwitchDescendant {
        match self {
            TsDescendant::Seq1 => TemplateSwitchDescendant::Reference,
            TsDescendant::Seq2 => TemplateSwitchDescendant::Query,
        }
    }
}

impl Display for TsKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TS{}{}",
            match self.ancestor {
                TsAncestor::Seq1 => "1",
                TsAncestor::Seq2 => "2",
            },
            match self.descendant {
                TsDescendant::Seq1 => "1",
                TsDescendant::Seq2 => "2",
            }
        )
    }
}
