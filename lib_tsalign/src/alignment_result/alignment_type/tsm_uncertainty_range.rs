/// Determines how the TSM uncertainty ranges are extended.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TSMUncertaintyRangeExtensionMode {
    /// Keep the TSM uncertainty range empty.
    None,

    /// Extend the TSM uncertainty range as far as possible without increasing cost.
    EqualCost,

    /// Extend the TSM uncertainty range as far as possible without increasing cost, but ignore cost increases caused by different TSM geometry.
    EqualCostIgnoreGeometry,
}

/// A heuristic range within which the start and end of the TS can be shifted.
/// This can be without increasing cost, or without introducing mismatches, or possibly other criteria in the future.
///
/// The range may not be maximal, but is required to be correct.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TSMUncertaintyRange {
    /// How much the start of the TS can be shifted to the left.
    pub start_left_shift: u16,

    /// How much the start of the TS can be shifted to the right.
    pub start_right_shift: u16,

    /// How much the end of the TS can be shifted to the left.
    pub end_left_shift: u16,

    /// How much the end of the TS can be shifted to the right.
    pub end_right_shift: u16,
}

impl TSMUncertaintyRange {
    pub const fn new() -> Self {
        Self {
            start_left_shift: 0,
            start_right_shift: 0,
            end_left_shift: 0,
            end_right_shift: 0,
        }
    }
}
