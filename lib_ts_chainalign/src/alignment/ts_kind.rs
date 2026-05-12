use std::fmt::Display;

use lib_tsalign::a_star_aligner::template_switch_distance::{
    TemplateSwitchPrimary, TemplateSwitchSecondary,
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
}

impl TsAncestor {
    pub fn into_tsalign_secondary(self) -> TemplateSwitchSecondary {
        match self {
            TsAncestor::Seq1 => TemplateSwitchSecondary::Reference,
            TsAncestor::Seq2 => TemplateSwitchSecondary::Query,
        }
    }
}

impl TsDescendant {
    pub fn into_tsalign_primary(self) -> TemplateSwitchPrimary {
        match self {
            TsDescendant::Seq1 => TemplateSwitchPrimary::Reference,
            TsDescendant::Seq2 => TemplateSwitchPrimary::Query,
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
