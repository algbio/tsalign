use std::fmt::Display;

use crate::{
    alignment::{
        sequences::AlignmentSequences,
        ts_kind::{TsAncestor, TsDescendant, TsKind},
    },
    anchors::secondary::SecondaryAnchor,
};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlignmentCoordinates {
    Primary(PrimaryAlignmentCoordinates),
    Secondary(SpecificSecondaryAlignmentCoordinates),
}

/// Alignment coordinates in the primary sequence space.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct PrimaryAlignmentCoordinates {
    a: usize,
    b: usize,
}

/// Alignment coordinates in the secondary sequence space, without specifying which secondary sequence space.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct AnySecondaryAlignmentCoordinates {
    /// Ancestor right index in the forward sequence.
    ancestor: usize,
    /// Descendant left index in the forward sequence.
    descendant: usize,
}

/// Alignment coordinates in a specified secondary sequence space.
#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct SpecificSecondaryAlignmentCoordinates {
    coordinates: AnySecondaryAlignmentCoordinates,
    ts_kind: TsKind,
}

impl AlignmentCoordinates {
    pub fn new_primary(a: usize, b: usize) -> Self {
        Self::Primary(PrimaryAlignmentCoordinates::new(a, b))
    }

    pub fn new_secondary(ancestor: usize, descendant: usize, ts_kind: TsKind) -> Self {
        Self::Secondary(SpecificSecondaryAlignmentCoordinates::new(
            ancestor, descendant, ts_kind,
        ))
    }

    pub fn primary_ordinate_a(&self) -> Option<usize> {
        match self {
            Self::Primary(primary) => Some(primary.a()),
            Self::Secondary { .. } => None,
        }
    }

    pub fn primary_ordinate_b(&self) -> Option<usize> {
        match self {
            Self::Primary(primary) => Some(primary.b()),
            Self::Secondary { .. } => None,
        }
    }

    /// Ancestor right index in the forward sequence.
    pub fn secondary_ordinate_ancestor(&self) -> Option<usize> {
        match self {
            Self::Secondary(secondary) => Some(secondary.ancestor()),
            Self::Primary { .. } => None,
        }
    }

    /// Descendant left index in the forward sequence.
    pub fn secondary_ordinate_descendant(&self) -> Option<usize> {
        match self {
            Self::Secondary(secondary) => Some(secondary.descendant()),
            Self::Primary { .. } => None,
        }
    }

    pub fn ts_kind(&self) -> Option<TsKind> {
        match self {
            Self::Secondary(secondary) => Some(secondary.ts_kind()),
            Self::Primary { .. } => None,
        }
    }

    pub fn is_primary(&self) -> bool {
        matches!(self, Self::Primary { .. })
    }

    pub fn is_secondary(&self) -> bool {
        matches!(self, Self::Secondary { .. })
    }

    pub fn into_primary(self) -> Option<PrimaryAlignmentCoordinates> {
        match self {
            Self::Primary(primary) => Some(primary),
            Self::Secondary { .. } => None,
        }
    }

    pub fn into_secondary(self) -> Option<SpecificSecondaryAlignmentCoordinates> {
        match self {
            Self::Primary { .. } => None,
            Self::Secondary(secondary) => Some(secondary),
        }
    }

    /// Checks if ordinate a can be incremented.
    /// In secondary alignments, ordinate a is the ancestor.
    ///
    /// If ordinate a and the `end` coordinates are both primary or both secondary, then the check is performed normally.
    /// If they differ, then there is a jump before the `end` boundary.
    pub fn can_increment_a_or_ancestor(
        &self,
        end: Self,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        match self {
            Self::Primary(primary) => primary.can_increment_a(end, sequences),
            Self::Secondary(secondary) => secondary.can_increment_ancestor(end, sequences),
        }
    }

    /// Checks if ordinate b can be incremented.
    /// In secondary alignments, ordinate b is the descendant.
    ///
    /// If ordinate b and the `end` coordinates are both primary or both secondary, then the check is performed normally.
    /// If the they differ, then there is a jump before the `end` boundary.
    pub fn can_increment_b_or_descendant(
        &self,
        end: Self,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        match self {
            Self::Primary(primary) => primary.can_increment_b(end, sequences),
            Self::Secondary(secondary) => secondary.can_increment_descendant(end, sequences),
        }
    }

    pub fn can_increment_both(&self, end: Self, sequences: Option<&AlignmentSequences>) -> bool {
        self.can_increment_a_or_ancestor(end, sequences)
            && self.can_increment_b_or_descendant(end, sequences)
    }

    pub fn increment_a(&self) -> Self {
        match self {
            Self::Primary(primary) => Self::Primary(primary.increment_a()),
            Self::Secondary(secondary) => Self::Secondary(secondary.increment_ancestor()),
        }
    }

    pub fn increment_b(&self) -> Self {
        match self {
            Self::Primary(primary) => Self::Primary(primary.increment_b()),
            Self::Secondary(secondary) => Self::Secondary(secondary.increment_descendant()),
        }
    }

    pub fn increment_both(&self) -> Self {
        self.increment_a().increment_b()
    }

    /// Generate all possible 12-jumps.
    ///
    /// The TS kind is given by the `start` coordinates.
    /// The left and right limits of the jump are given by the `start` and `end` coordinates.
    /// The `end` coordinates must be in primary form and simply be the end of the aligned sequences.
    /// The `start` coordinates are in secondary form.
    pub fn generate_12_jumps(
        &self,
        start: SpecificSecondaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
    ) -> impl Iterator<Item = (isize, AnySecondaryAlignmentCoordinates)> {
        let Self::Primary(primary) = self else {
            panic!("Can only generate 12-jumps from primary coordinates");
        };

        primary.generate_12_jumps(start, end)
    }

    /// Generate all possible 34-jumps.
    ///
    /// The `end` coordinates are in primary form and limit the jump to the left of (or into) them.
    pub fn generate_34_jumps(
        &self,
        end: PrimaryAlignmentCoordinates,
    ) -> impl Iterator<Item = (isize, PrimaryAlignmentCoordinates)> {
        let Self::Secondary(secondary) = self else {
            panic!("Can only generate 34-jumps from secondary coordinates");
        };

        secondary.generate_34_jumps(end)
    }
}

impl PrimaryAlignmentCoordinates {
    pub fn new(a: usize, b: usize) -> Self {
        Self { a, b }
    }

    pub fn a(&self) -> usize {
        self.a
    }

    pub fn b(&self) -> usize {
        self.b
    }

    pub fn increment_a(self) -> Self {
        Self {
            a: self.a + 1,
            b: self.b,
        }
    }

    pub fn increment_b(self) -> Self {
        Self {
            a: self.a,
            b: self.b + 1,
        }
    }

    pub fn increment_both(self, increment: usize) -> Self {
        Self {
            a: self.a + increment,
            b: self.b + increment,
        }
    }

    /// Checks if ordinate a can be incremented.
    ///
    /// If the `end` coordinates are primary, then the check is performed normally.
    /// If the `end` coordinates are secondary, then there is a jump before the `end` boundary.
    pub fn can_increment_a(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        match end {
            AlignmentCoordinates::Primary(end) => self.can_increment_a_primary(end),
            AlignmentCoordinates::Secondary(end) => self.can_increment_a_secondary(end, sequences),
        }
    }

    /// Checks if the ordinate a can be incremented against primary coordinates.
    pub fn can_increment_a_primary(&self, end: PrimaryAlignmentCoordinates) -> bool {
        // Incrementing primary ordinate a is always a plus operation, so we only need to check the upper bound.
        self.a() < end.a()
    }

    /// Checks if the ordinate a can be incremented against secondary coordinates.
    pub fn can_increment_a_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        if let Some(sequences) = sequences {
            match end.ts_kind().descendant {
                // Descendant is a, so it is limited by the descendant ordinate.
                TsDescendant::Seq1 => self.a() < end.descendant(),
                // Descendant is b, so a can go until the end of the sequence.
                TsDescendant::Seq2 => self.a() < sequences.primary_end().a(),
            }
        } else {
            true
        }
    }

    /// Checks if ordinate b can be incremented.
    ///
    /// If the `end` coordinates are primary, then the check is performed normally.
    /// If the `end` coordinates are secondary, then there is a jump before the `end` boundary.
    pub fn can_increment_b(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        match end {
            AlignmentCoordinates::Primary(end) => self.can_increment_b_primary(end),
            AlignmentCoordinates::Secondary(end) => self.can_increment_b_secondary(end, sequences),
        }
    }

    /// Checks if the ordinate b can be incremented against primary coordinates.
    pub fn can_increment_b_primary(&self, end: PrimaryAlignmentCoordinates) -> bool {
        // Incrementing ordinate b is always a plus operation, so we only need to check the upper bound.
        self.b() < end.b()
    }

    /// Checks if the ordinate b can be incremented against secondary coordinates.
    pub fn can_increment_b_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        if let Some(sequences) = sequences {
            match end.ts_kind().descendant {
                // Descendant is a, so b can go until the end of the sequence.
                TsDescendant::Seq1 => self.b() < sequences.primary_end().b(),
                // Descendant is b, so it is limited by the descendant ordinate.
                TsDescendant::Seq2 => self.b() < end.descendant(),
            }
        } else {
            true
        }
    }

    pub fn can_increment_both(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        self.can_increment_a(end, sequences) && self.can_increment_b(end, sequences)
    }

    pub fn can_increment_both_primary(&self, end: PrimaryAlignmentCoordinates) -> bool {
        self.can_increment_a_primary(end) && self.can_increment_b_primary(end)
    }

    pub fn can_increment_both_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        self.can_increment_a_secondary(end, sequences)
            && self.can_increment_b_secondary(end, sequences)
    }

    /// Generate all possible 12-jumps.
    ///
    /// The TS kind is given by the `start` coordinates.
    /// The left and right limits of the jump are given by the `start` and `end` coordinates.
    /// The `end` coordinates must be in primary form and simply be the end of the aligned sequences.
    /// The `start` coordinates are in secondary form.
    pub fn generate_12_jumps(
        &self,
        start: SpecificSecondaryAlignmentCoordinates,
        end: Self,
    ) -> impl Iterator<Item = (isize, AnySecondaryAlignmentCoordinates)> {
        let Self { a, b, .. } = *self;
        let ts_kind = start.ts_kind();
        let ancestor_zero = match ts_kind.ancestor {
            TsAncestor::Seq1 => a,
            TsAncestor::Seq2 => b,
        } as isize;
        let ancestor_limit = match ts_kind.ancestor {
            TsAncestor::Seq1 => end.a(),
            TsAncestor::Seq2 => end.b(),
        };
        let descendant = match ts_kind.descendant {
            TsDescendant::Seq1 => a,
            TsDescendant::Seq2 => b,
        };

        (start.ancestor()..=ancestor_limit).map(move |ancestor| {
            (
                ancestor as isize - ancestor_zero,
                AnySecondaryAlignmentCoordinates::new(ancestor, descendant),
            )
        })
    }
}

impl AnySecondaryAlignmentCoordinates {
    pub fn new(ancestor: usize, descendant: usize) -> Self {
        Self {
            ancestor,
            descendant,
        }
    }

    pub fn into_specific(self, ts_kind: TsKind) -> SpecificSecondaryAlignmentCoordinates {
        SpecificSecondaryAlignmentCoordinates::new(self.ancestor, self.descendant, ts_kind)
    }

    /// Ancestor right index in the forward sequence.
    pub fn ancestor(&self) -> usize {
        self.ancestor
    }

    /// Descendant left index in the forward sequence.
    pub fn descendant(&self) -> usize {
        self.descendant
    }

    pub fn increment_ancestor(self) -> Self {
        Self {
            ancestor: self.ancestor.wrapping_sub(1),
            descendant: self.descendant,
        }
    }

    pub fn increment_descendant(self) -> Self {
        Self {
            ancestor: self.ancestor,
            descendant: self.descendant + 1,
        }
    }

    pub fn increment_both(self, increment: usize) -> Self {
        Self {
            // This was checked sub before, but since increment_a is wrapping, this most likely should also be wrapping.
            ancestor: self.ancestor.wrapping_sub(increment),
            descendant: self.descendant + increment,
        }
    }

    /// Checks if the ancestor can be incremented.
    ///
    /// If the `end` coordinates are secondary, then the check is performed normally.
    /// If the `end` coordinates are primary, then there is a jump before the `end` boundary.
    pub fn can_increment_ancestor(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
        ts_kind: TsKind,
    ) -> bool {
        match end {
            AlignmentCoordinates::Primary(_) => self.can_increment_ancestor_primary(sequences),
            AlignmentCoordinates::Secondary(end) => {
                self.can_increment_ancestor_secondary(end, ts_kind)
            }
        }
    }

    /// Checks if the ancestor can be incremented against primary coordinates.
    pub fn can_increment_ancestor_primary(&self, sequences: Option<&AlignmentSequences>) -> bool {
        if sequences.is_some() {
            0 < self.ancestor()
        } else {
            true
        }
    }

    /// Checks if the ancestor can be incremented against secondary coordinates.
    pub fn can_increment_ancestor_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> bool {
        assert_eq!(ts_kind, end.ts_kind());
        // Incrementing the secondary ancestor is always a minus operation, so we only need to check the lower bound.
        end.ancestor() < self.ancestor()
    }

    /// Checks if the descendant can be incremented.
    ///
    /// If the `end` coordinates are secondary, then the check is performed normally.
    /// If the `end` coordinates are primary, then there is a jump before the `end` boundary.
    pub fn can_increment_descendant(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
        ts_kind: TsKind,
    ) -> bool {
        match end {
            AlignmentCoordinates::Primary(end) => {
                self.can_increment_descendant_primary(end, sequences, ts_kind)
            }
            AlignmentCoordinates::Secondary(end) => {
                self.can_increment_descendant_secondary(end, ts_kind)
            }
        }
    }

    /// Checks if the descendant can be incremented against primary coordinates.
    pub fn can_increment_descendant_primary(
        &self,
        end: PrimaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
        ts_kind: TsKind,
    ) -> bool {
        if sequences.is_some() {
            match ts_kind.descendant {
                TsDescendant::Seq1 => self.descendant() < end.a(),
                TsDescendant::Seq2 => self.descendant() < end.b(),
            }
        } else {
            true
        }
    }

    /// Checks if the descendant can be incremented against secondary coordinates.
    pub fn can_increment_descendant_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> bool {
        assert_eq!(ts_kind, end.ts_kind());
        // Incrementing ordinate b is always a plus operation, so we only need to check the upper bound.
        self.descendant() < end.descendant()
    }

    pub fn can_increment_both(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
        ts_kind: TsKind,
    ) -> bool {
        self.can_increment_ancestor(end, sequences, ts_kind)
            && self.can_increment_descendant(end, sequences, ts_kind)
    }

    pub fn can_increment_both_primary(
        &self,
        end: PrimaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
        ts_kind: TsKind,
    ) -> bool {
        self.can_increment_ancestor_primary(sequences)
            && self.can_increment_descendant_primary(end, sequences, ts_kind)
    }

    pub fn can_increment_both_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> bool {
        self.can_increment_ancestor_secondary(end, ts_kind)
            && self.can_increment_descendant_secondary(end, ts_kind)
    }

    pub fn cmp<Cost>(&self, anchor: &SecondaryAnchor<Cost>) -> std::cmp::Ordering {
        // The comparison is independent of the ts_kind, so we just choose any.
        Ord::cmp(self, &anchor.start())
    }

    /// Generate all possible 34-jumps.
    ///
    /// The `end` coordinates are in primary form and limit the jump to the left of (or into) them.
    pub fn generate_34_jumps(
        &self,
        end: PrimaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (isize, PrimaryAlignmentCoordinates)> {
        let Self {
            ancestor,
            descendant,
        } = *self;

        (0..=match ts_kind.descendant {
            TsDescendant::Seq1 => end.b(),
            TsDescendant::Seq2 => end.a(),
        })
            .map(move |new_ancestor| {
                (
                    new_ancestor as isize - ancestor as isize,
                    match ts_kind.descendant {
                        TsDescendant::Seq1 => PrimaryAlignmentCoordinates {
                            a: descendant,
                            b: new_ancestor,
                        },
                        TsDescendant::Seq2 => PrimaryAlignmentCoordinates {
                            a: new_ancestor,
                            b: descendant,
                        },
                    },
                )
            })
    }
}

impl SpecificSecondaryAlignmentCoordinates {
    pub fn new(ancestor: usize, descendant: usize, ts_kind: TsKind) -> Self {
        Self {
            coordinates: AnySecondaryAlignmentCoordinates::new(ancestor, descendant),
            ts_kind,
        }
    }

    /// Ancestor right index in the forward sequence.
    pub fn ancestor(&self) -> usize {
        self.coordinates.ancestor()
    }

    /// Descendant left index in the forward sequence.
    pub fn descendant(&self) -> usize {
        self.coordinates.descendant()
    }

    pub fn ts_kind(&self) -> TsKind {
        self.ts_kind
    }

    pub fn increment_ancestor(self) -> Self {
        Self {
            coordinates: self.coordinates.increment_ancestor(),
            ts_kind: self.ts_kind,
        }
    }

    pub fn increment_descendant(self) -> Self {
        Self {
            coordinates: self.coordinates.increment_descendant(),
            ts_kind: self.ts_kind,
        }
    }

    pub fn increment_both(self, increment: usize) -> Self {
        Self {
            coordinates: self.coordinates.increment_both(increment),
            ts_kind: self.ts_kind,
        }
    }

    /// Checks if the ancestor can be incremented.
    ///
    /// If the `end` coordinates are secondary, then the check is performed normally.
    /// If the `end` coordinates are primary, then there is a jump before the `end` boundary.
    pub fn can_increment_ancestor(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self).can_increment_ancestor(
            end,
            sequences,
            self.ts_kind,
        )
    }

    /// Checks if the ancestor can be incremented against primary coordinates.
    pub fn can_increment_ancestor_primary(&self, sequences: Option<&AlignmentSequences>) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self).can_increment_ancestor_primary(sequences)
    }

    /// Checks if the ancestor can be incremented against secondary coordinates.
    pub fn can_increment_ancestor_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
    ) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self)
            .can_increment_ancestor_secondary(end, self.ts_kind)
    }

    /// Checks if the descendant can be incremented.
    ///
    /// If the `end` coordinates are secondary, then the check is performed normally.
    /// If the `end` coordinates are primary, then there is a jump before the `end` boundary.
    pub fn can_increment_descendant(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self).can_increment_descendant(
            end,
            sequences,
            self.ts_kind,
        )
    }

    /// Checks if the descendant can be incremented against primary coordinates.
    pub fn can_increment_descendant_primary(
        &self,
        end: PrimaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self).can_increment_descendant_primary(
            end,
            sequences,
            self.ts_kind,
        )
    }

    /// Checks if the descendant can be incremented against secondary coordinates.
    pub fn can_increment_descendant_secondary(
        &self,
        end: SpecificSecondaryAlignmentCoordinates,
    ) -> bool {
        AnySecondaryAlignmentCoordinates::from(*self)
            .can_increment_descendant_secondary(end, self.ts_kind)
    }

    pub fn can_increment_both(
        &self,
        end: AlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        self.can_increment_ancestor(end, sequences) && self.can_increment_descendant(end, sequences)
    }

    pub fn can_increment_both_primary(
        &self,
        end: PrimaryAlignmentCoordinates,
        sequences: Option<&AlignmentSequences>,
    ) -> bool {
        self.can_increment_ancestor_primary(sequences)
            && self.can_increment_descendant_primary(end, sequences)
    }

    pub fn can_increment_both_secondary(&self, end: SpecificSecondaryAlignmentCoordinates) -> bool {
        self.can_increment_ancestor_secondary(end) && self.can_increment_descendant_secondary(end)
    }

    pub fn cmp<Cost>(&self, anchor: &SecondaryAnchor<Cost>) -> std::cmp::Ordering {
        Ord::cmp(
            &AnySecondaryAlignmentCoordinates::from(*self),
            &anchor.start(),
        )
    }

    /// Generate all possible 34-jumps.
    ///
    /// The `end` coordinates are in primary form and limit the jump to the left of (or into) them.
    pub fn generate_34_jumps(
        &self,
        end: PrimaryAlignmentCoordinates,
    ) -> impl Iterator<Item = (isize, PrimaryAlignmentCoordinates)> {
        self.coordinates.generate_34_jumps(end, self.ts_kind)
    }
}

impl Display for AlignmentCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentCoordinates::Primary(primary) => write!(f, "{primary}"),
            AlignmentCoordinates::Secondary(secondary) => write!(f, "{secondary}"),
        }
    }
}

impl Display for PrimaryAlignmentCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "P({}, {})", self.a, self.b)
    }
}

impl Display for AnySecondaryAlignmentCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "S({}, {})", self.ancestor(), self.descendant())
    }
}

impl Display for SpecificSecondaryAlignmentCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "S({}, {}, {})",
            self.ancestor(),
            self.descendant(),
            self.ts_kind,
        )
    }
}

impl From<PrimaryAlignmentCoordinates> for AlignmentCoordinates {
    fn from(value: PrimaryAlignmentCoordinates) -> Self {
        Self::Primary(value)
    }
}

impl From<SpecificSecondaryAlignmentCoordinates> for AlignmentCoordinates {
    fn from(value: SpecificSecondaryAlignmentCoordinates) -> Self {
        Self::Secondary(value)
    }
}

impl From<SpecificSecondaryAlignmentCoordinates> for AnySecondaryAlignmentCoordinates {
    fn from(value: SpecificSecondaryAlignmentCoordinates) -> Self {
        value.coordinates
    }
}

impl From<&'_ PrimaryAlignmentCoordinates> for AlignmentCoordinates {
    fn from(value: &PrimaryAlignmentCoordinates) -> Self {
        Self::Primary(*value)
    }
}

impl From<&'_ SpecificSecondaryAlignmentCoordinates> for AlignmentCoordinates {
    fn from(value: &SpecificSecondaryAlignmentCoordinates) -> Self {
        Self::Secondary(*value)
    }
}

impl From<&'_ SpecificSecondaryAlignmentCoordinates> for AnySecondaryAlignmentCoordinates {
    fn from(value: &SpecificSecondaryAlignmentCoordinates) -> Self {
        value.coordinates
    }
}
