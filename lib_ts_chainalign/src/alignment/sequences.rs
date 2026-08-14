use crate::alignment::{
    coordinates::{
        AlignmentCoordinates, PrimaryAlignmentCoordinates, SpecificSecondaryAlignmentCoordinates,
    },
    ts_kind::{TsAncestor, TsDescendant, TsKind},
};

pub struct AlignmentSequences {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
    seq1_name: String,
    seq2_name: String,
    start: PrimaryAlignmentCoordinates,
    end: PrimaryAlignmentCoordinates,
}

impl AlignmentSequences {
    pub fn new(
        seq1: Vec<u8>,
        seq2: Vec<u8>,
        start: PrimaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
    ) -> Self {
        Self::new_named(
            seq1,
            seq2,
            "seq1".to_string(),
            "seq2".to_string(),
            start,
            end,
        )
    }

    pub fn new_complete(seq1: Vec<u8>, seq2: Vec<u8>) -> Self {
        let end = PrimaryAlignmentCoordinates::new(seq1.len(), seq2.len());
        Self::new(seq1, seq2, PrimaryAlignmentCoordinates::new(0, 0), end)
    }

    pub fn new_named(
        seq1: Vec<u8>,
        seq2: Vec<u8>,
        seq1_name: String,
        seq2_name: String,
        start: PrimaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
    ) -> Self {
        Self {
            seq1,
            seq2,
            seq1_name,
            seq2_name,
            start,
            end,
        }
    }

    pub fn characters(
        &self,
        coordinates: AlignmentCoordinates,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> (u8, u8) {
        match coordinates {
            AlignmentCoordinates::Primary(primary) => self.primary_characters(primary),
            AlignmentCoordinates::Secondary(secondary) => {
                self.secondary_characters(secondary, rc_fn)
            }
        }
    }

    pub fn primary_characters(&self, coordinates: PrimaryAlignmentCoordinates) -> (u8, u8) {
        (self.seq1[coordinates.a()], self.seq2[coordinates.b()])
    }

    pub fn secondary_characters(
        &self,
        coordinates: SpecificSecondaryAlignmentCoordinates,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> (u8, u8) {
        (
            match coordinates.ts_kind().ancestor {
                TsAncestor::Seq1 => self.seq1[coordinates.ancestor() - 1],
                TsAncestor::Seq2 => self.seq2[coordinates.ancestor() - 1],
            },
            rc_fn(match coordinates.ts_kind().descendant {
                TsDescendant::Seq1 => self.seq1[coordinates.descendant()],
                TsDescendant::Seq2 => self.seq2[coordinates.descendant()],
            }),
        )
    }

    pub fn primary_start(&self) -> PrimaryAlignmentCoordinates {
        self.start
    }

    pub fn primary_end(&self) -> PrimaryAlignmentCoordinates {
        self.end
    }

    pub fn secondary_end(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        match ts_kind {
            ts_kind @ (TsKind::TS11 | TsKind::TS21) => {
                AlignmentCoordinates::new_secondary(0, self.seq1.len(), ts_kind)
            }
            ts_kind @ (TsKind::TS12 | TsKind::TS22) => {
                AlignmentCoordinates::new_secondary(0, self.seq2.len(), ts_kind)
            }
        }
    }

    pub fn end(&self, ts_kind: Option<TsKind>) -> AlignmentCoordinates {
        match ts_kind {
            None => self.primary_end().into(),
            Some(ts_kind) => self.secondary_end(ts_kind),
        }
    }

    /// Returns the full sequence1, without restricting to the alignment range.
    pub fn seq1(&self) -> &[u8] {
        &self.seq1
    }

    /// Returns the full sequence2, without restricting to the alignment range.
    pub fn seq2(&self) -> &[u8] {
        &self.seq2
    }

    /// Returns the name of sequence1.
    pub fn seq1_name(&self) -> &str {
        &self.seq1_name
    }

    /// Returns the name of sequence2.
    pub fn seq2_name(&self) -> &str {
        &self.seq2_name
    }
}
