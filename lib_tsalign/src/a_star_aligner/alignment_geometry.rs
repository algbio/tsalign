use std::{fmt::Display, ops::Range};

#[derive(Debug, Clone, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub struct AlignmentRange {
    offset: AlignmentCoordinates,
    limit: AlignmentCoordinates,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub struct AlignmentCoordinates {
    reference: usize,
    query: usize,
}

impl AlignmentRange {
    pub fn new(
        reference_offset: usize,
        query_offset: usize,
        reference_limit: usize,
        query_limit: usize,
    ) -> Self {
        Self {
            offset: AlignmentCoordinates::new(reference_offset, query_offset),
            limit: AlignmentCoordinates::new(reference_limit, query_limit),
        }
    }

    pub fn new_complete(reference: usize, query: usize) -> Self {
        Self {
            offset: AlignmentCoordinates::new_zero(),
            limit: AlignmentCoordinates::new(reference, query),
        }
    }

    pub fn new_offset_limit(offset: AlignmentCoordinates, limit: AlignmentCoordinates) -> Self {
        Self { offset, limit }
    }

    pub fn reference_offset(&self) -> usize {
        self.offset.reference
    }

    pub fn query_offset(&self) -> usize {
        self.offset.query
    }

    pub fn reference_limit(&self) -> usize {
        self.limit.reference
    }

    pub fn query_limit(&self) -> usize {
        self.limit.query
    }

    pub fn reference_range(&self) -> Range<usize> {
        self.offset.reference..self.limit.reference
    }

    pub fn query_range(&self) -> Range<usize> {
        self.offset.query..self.limit.query
    }

    #[must_use]
    pub fn move_offsets_left(&self) -> Option<Self> {
        Some(Self::new_offset_limit(
            AlignmentCoordinates::new(
                self.offset.reference.checked_sub(1)?,
                self.offset.query.checked_sub(1)?,
            ),
            self.limit,
        ))
    }

    #[must_use]
    pub fn move_offsets_right(&self) -> Option<Self> {
        let new_ref_offset = self.offset.reference.checked_add(1)?;
        let new_qry_offset = self.offset.query.checked_add(1)?;

        if new_ref_offset > self.limit.reference || new_qry_offset > self.limit.query {
            return None;
        }

        Some(Self::new_offset_limit(
            AlignmentCoordinates::new(new_ref_offset, new_qry_offset),
            self.limit,
        ))
    }

    #[must_use]
    pub fn move_limits_left(&self) -> Option<Self> {
        let new_ref_limit = self.limit.reference.checked_sub(1)?;
        let new_qry_limit = self.limit.query.checked_sub(1)?;

        if new_ref_limit < self.offset.reference || new_qry_limit > self.offset.query {
            return None;
        }
        Some(Self::new_offset_limit(
            self.offset,
            AlignmentCoordinates::new(new_ref_limit, new_qry_limit),
        ))
    }

    #[must_use]
    pub fn move_limits_right(&self) -> Option<Self> {
        Some(Self::new_offset_limit(
            self.offset,
            AlignmentCoordinates::new(
                self.limit.reference.checked_add(1)?,
                self.limit.query.checked_add(1)?,
            ),
        ))
    }
}

impl AlignmentCoordinates {
    pub fn new(reference: usize, query: usize) -> Self {
        Self { reference, query }
    }

    pub fn new_zero() -> Self {
        Self::new(0, 0)
    }
}

impl Display for AlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "R: {}..{}; Q: {}..{}",
            self.offset.reference, self.limit.reference, self.offset.query, self.limit.query
        )
    }
}
