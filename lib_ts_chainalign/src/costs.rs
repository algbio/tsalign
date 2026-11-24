use std::{fmt::Debug, ops::Range};

use generic_a_star::cost::AStarCost;
use num_traits::{Zero, bounds::UpperBounded};
use serde::{Deserialize, Serialize};

use crate::alignment::{AlignmentType, ts_kind::TsKind};

mod compat;

pub use compat::CostConversionError;

/// The default cost configuration of tschainalign in the plain `.tsa` format of
/// [`TemplateSwitchConfig::read_plain`](lib_tsalign::config::TemplateSwitchConfig::read_plain).
///
/// The classic default config ([`lib_tsalign::config::DEFAULT_CONFIG_TSA`]) is not supported by tschainalign,
/// so the CLI (without `-c`) and the Python bindings (without costs) use this one instead:
/// only reverse template switches, over the DNA alphabet without `N`.
/// It is `lib_ts_chainalign/default_config.tsa` of the repository, which equals `test_files/config/chainalign/config.tsa`.
pub const DEFAULT_CONFIG_TSA: &str = include_str!("../default_config.tsa");

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

    pub fn has_zero_cost(&self) -> bool
    where
        Cost: Zero,
    {
        self.substitution.is_zero() || self.gap_open.is_zero() || self.gap_extend.is_zero()
    }

    pub fn gap_cost(&self, gap_length: usize) -> Cost
    where
        Cost: AStarCost,
    {
        if gap_length == 0 {
            return Cost::zero();
        }
        self.gap_open + self.gap_extend * Cost::from_usize(gap_length - 1)
    }

    pub fn into_costs<ToCost: From<Cost>>(self) -> GapAffineCosts<ToCost> {
        GapAffineCosts::new(
            ToCost::from(self.substitution),
            ToCost::from(self.gap_open),
            ToCost::from(self.gap_extend),
        )
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
        let mut cost_by_kind = [None, None, None, None];

        for (ts_kind, cost) in iter {
            assert!(
                cost_by_kind[ts_kind.index()].is_none(),
                "Duplicate TsKind {ts_kind} in iterator",
            );
            cost_by_kind[ts_kind.index()] = Some(cost);
        }

        Self::new(cost_by_kind.map(|cost| cost.expect("Missing TsKind")))
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

impl<Cost: AStarCost> AlignmentCosts<Cost> {
    /// Returns true if merging two adjacent gaps into one never increases the cost,
    /// i.e. if the gap extend cost is at most the gap open cost, both for primary and secondary alignments.
    pub fn gap_merging_never_increases_cost(&self) -> bool {
        self.primary_costs.gap_extend <= self.primary_costs.gap_open
            && self.secondary_costs.gap_extend <= self.secondary_costs.gap_open
    }

    /// Computes the cost of a complete alignment, given as a sequence of single alignment operations.
    ///
    /// This is the cost of the alignment as a whole, as it is defined by the classic template switch aligner:
    /// consecutive gap characters of the same kind form a single gap,
    /// which pays the gap open cost for its first character and the gap extend cost for every further character.
    /// The gap state is reset at the start and at the end of a template switch.
    /// A template switch costs the base cost of its kind.
    ///
    /// The geometry of template switches is not charged: tschainalign uses only cost functions for the geometry
    /// that are zero inside of [`TsLimits`] and infinite outside, and checks the limits only based on the anchors.
    /// The jump aligners may make a template switch longer or shorter than these limits allow
    /// (see [issue #207](https://github.com/algbio/tsalign/issues/207)).
    ///
    /// Note that this may be lower than the cost of the chain that the alignment was built from,
    /// since the chain cost is the sum of the costs of its gaps and anchors, each of which is aligned separately.
    /// See [issue #221](https://github.com/algbio/tsalign/issues/221).
    pub fn alignment_cost(&self, alignment: impl IntoIterator<Item = AlignmentType>) -> Cost {
        let mut cost = Cost::zero();
        let mut is_primary = true;
        let mut previous = None;

        for alignment_type in alignment {
            let gap_affine_costs = if is_primary {
                &self.primary_costs
            } else {
                &self.secondary_costs
            };

            cost += match alignment_type {
                AlignmentType::Match => Cost::zero(),
                AlignmentType::Substitution => gap_affine_costs.substitution,
                AlignmentType::GapA | AlignmentType::GapB => {
                    if previous == Some(alignment_type) {
                        gap_affine_costs.gap_extend
                    } else {
                        gap_affine_costs.gap_open
                    }
                }
                AlignmentType::TsStart { ts_kind, .. } => {
                    assert!(
                        is_primary,
                        "Template switch started within a template switch"
                    );
                    is_primary = false;
                    self.ts_base_cost.get(ts_kind)
                }
                AlignmentType::TsEnd { .. } => {
                    assert!(
                        !is_primary,
                        "Template switch ended outside of a template switch"
                    );
                    is_primary = true;
                    Cost::zero()
                }
            };
            previous = Some(alignment_type);
        }

        assert!(is_primary, "Alignment ends within a template switch");
        cost
    }
}

impl<Character: Eq, Cost> inexact_alignment_anchors::costs::AlignmentCost<Character, Cost>
    for GapAffineCosts<Cost>
where
    Cost: AStarCost + Copy,
{
    fn cost(&self, a: &Character, b: &Character) -> Cost {
        if a == b {
            Cost::zero()
        } else {
            self.substitution
        }
    }

    fn gap_open_cost(&self) -> Cost {
        self.gap_open
    }

    fn gap_extend_cost(&self) -> Cost {
        self.gap_extend
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    use generic_a_star::cost::U32Cost;
    use lib_tsalign::config::TemplateSwitchConfig;

    use super::{AlignmentCosts, DEFAULT_CONFIG_TSA};

    #[test]
    fn default_config_is_supported() {
        let config =
            TemplateSwitchConfig::<DnaAlphabet, U32Cost>::read_plain(DEFAULT_CONFIG_TSA.as_bytes())
                .unwrap();
        AlignmentCosts::try_from(config).unwrap();
    }

    #[test]
    fn default_config_equals_test_files_config() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../test_files/config/chainalign/config.tsa");
        assert_eq!(
            DEFAULT_CONFIG_TSA,
            std::fs::read_to_string(&path).unwrap(),
            "the built-in default config of tschainalign differs from {}",
            path.display()
        );
    }
}
