use std::ops::Range;

use num_traits::{Zero, bounds::UpperBounded};
use serde::{Deserialize, Serialize};

use crate::alignment::ts_kind::TsKind;

mod compat;

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct GapAffineCosts<Cost> {
    pub substitution: Cost,
    pub gap_open: Cost,
    pub gap_extend: Cost,
}

#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct TsLimits {
    /// The maximum range of the 12-jump of an inter (12 or 21) template switch.
    /// This parameter is ignored for now.
    pub inter_jump_12: Range<isize>,
    /// The maximum range of the 12-jump of an intra (11 or 22) template switch.
    /// This parameter is ignored for now.
    pub intra_jump_12: Range<isize>,
    /// The maximum range of the 34-jump of a template switch.
    /// This parameter is ignored for now.
    pub jump_34: Range<isize>,
    /// The range for the length of the 23-alignment of a template switch.
    pub length_23: Range<usize>,
    /// The range for the ancestor gap of a template switch.
    /// This parameter is ignored for now.
    pub ancestor_gap: Range<isize>,
}

/// The base cost of a template switch.
/// This is applied whenever a template switch is started.
#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct TsBaseCost<Cost> {
    cost_by_kind: [Cost; 4],
}

/// The cost function for alignments.
///
/// For convenience, it implements [`TryFrom<TemplateSwitchConfig<Alphabet, Cost>>`](std::convert::TryFrom).
/// Note that the conversion is very strict and only allows to convert from a [`TemplateSwitchConfig`](lib_tsalign::config::TemplateSwitchConfig) if the conversion loses no information.
#[derive(Debug, Clone, Serialize, Deserialize, Eq, PartialEq)]
pub struct AlignmentCosts<Cost> {
    /// Costs for primary alignment outside of template switches.
    pub primary_costs: GapAffineCosts<Cost>,
    /// Costs for secondary alignment, i.e. for the 23-alignment of a template switch.
    pub secondary_costs: GapAffineCosts<Cost>,
    /// The base cost of a template switch.
    /// This is applied whenever a template switch is started.
    pub ts_base_cost: TsBaseCost<Cost>,
    /// Limits on the geometry of a template switch.
    pub ts_limits: TsLimits,
}

impl<Cost> GapAffineCosts<Cost> {
    pub fn new(substitution: Cost, gap_open: Cost, gap_extend: Cost) -> Self {
        Self {
            substitution,
            gap_open,
            gap_extend,
        }
    }
}

impl<Cost: Zero> GapAffineCosts<Cost> {
    pub fn has_zero_cost(&self) -> bool {
        self.substitution.is_zero() || self.gap_open.is_zero() || self.gap_extend.is_zero()
    }
}

impl TsLimits {
    pub fn new_unlimited() -> Self {
        Self {
            inter_jump_12: isize::MIN..isize::MAX,
            intra_jump_12: isize::MIN..isize::MAX,
            jump_34: isize::MIN..isize::MAX,
            length_23: usize::MIN..usize::MAX,
            ancestor_gap: isize::MIN..isize::MAX,
        }
    }
}

impl<Cost> TsBaseCost<Cost> {
    pub fn new(cost_by_kind: [Cost; 4]) -> Self {
        Self { cost_by_kind }
    }

    pub fn has_zero_cost(&self) -> bool
    where
        Cost: Zero,
    {
        self.cost_by_kind.iter().any(|cost| cost.is_zero())
    }

    pub fn min(&self) -> Cost
    where
        Cost: Copy + Ord,
    {
        *self.cost_by_kind.iter().min().unwrap()
    }

    pub fn get(&self, ts_kind: TsKind) -> Cost
    where
        Cost: Copy,
    {
        self.cost_by_kind[ts_kind.index()]
    }
}

impl<Cost: UpperBounded + Eq> FromIterator<(TsKind, Cost)> for TsBaseCost<Cost> {
    fn from_iter<T: IntoIterator<Item = (TsKind, Cost)>>(iter: T) -> Self {
        let mut cost_by_kind = [
            Cost::max_value(),
            Cost::max_value(),
            Cost::max_value(),
            Cost::max_value(),
        ];

        let mut count = 0;
        for (ts_kind, cost) in iter {
            assert!(
                cost_by_kind[ts_kind.index()] == Cost::max_value(),
                "Duplicate TsKind {ts_kind} in iterator",
            );
            cost_by_kind[ts_kind.index()] = cost;
            count += 1;
        }
        assert_eq!(
            count, 4,
            "Iterator must contain exactly 4 entries, one for each TsKind",
        );

        Self::new(cost_by_kind)
    }
}

impl<Cost: Copy> From<Cost> for TsBaseCost<Cost> {
    fn from(value: Cost) -> Self {
        Self::new([value; 4])
    }
}

impl<Cost> AlignmentCosts<Cost> {
    pub fn new(
        primary_costs: GapAffineCosts<Cost>,
        secondary_costs: GapAffineCosts<Cost>,
        ts_base_cost: TsBaseCost<Cost>,
        ts_limits: TsLimits,
    ) -> Self {
        Self {
            primary_costs,
            secondary_costs,
            ts_base_cost,
            ts_limits,
        }
    }
}

impl<Cost: Zero> AlignmentCosts<Cost> {
    pub fn has_zero_cost(&self) -> bool {
        self.primary_costs.has_zero_cost()
            || self.secondary_costs.has_zero_cost()
            || self.ts_base_cost.has_zero_cost()
    }
}
