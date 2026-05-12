//! Functions for making tschainalign compatible with tsalign types.

use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;
use lib_tsalign::{
    a_star_aligner::template_switch_distance::TemplateSwitchDirection,
    config::TemplateSwitchConfig, costs::gap_affine::GapAffineAlignmentCostTable,
};

use crate::{
    alignment::ts_kind::TsKind,
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
};

impl<AlphabetType: Alphabet, Cost: AStarCost> From<TemplateSwitchConfig<AlphabetType, Cost>>
    for AlignmentCosts<Cost>
{
    fn from(value: TemplateSwitchConfig<AlphabetType, Cost>) -> Self {
        let value = &value;
        assert!(value.left_flank_length == 0 && value.right_flank_length == 0);

        let ts_base_cost = TsBaseCost::from_iter([
            (TsKind::TS11, value.base_cost.rrr),
            (TsKind::TS12, value.base_cost.qrr), // Note that the ancestor is the secondary and the descendant is the primary, hence these are flipped.
            (TsKind::TS21, value.base_cost.rqr),
            (TsKind::TS22, value.base_cost.qqr),
        ]);

        Self {
            primary_costs: (&value.primary_edit_costs).into(),
            secondary_costs: value
                .secondary_edit_costs(TemplateSwitchDirection::Reverse)
                .into(),
            ts_base_cost,
            ts_limits: TsLimits {
                inter_jump_12: value.rq_qr_offset_costs.zero_range().unwrap(),
                intra_jump_12: value.rr_qq_offset_costs.zero_range().unwrap(),
                // tsalign costs do not support limiting this.
                jump_34: isize::MIN..isize::MAX,
                length_23: value.length_costs.zero_range().unwrap(),
                ancestor_gap: value.reverse_anti_primary_gap_costs.zero_range().unwrap(),
            },
        }
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost>
    From<&'_ GapAffineAlignmentCostTable<AlphabetType, Cost>> for GapAffineCosts<Cost>
{
    fn from(value: &GapAffineAlignmentCostTable<AlphabetType, Cost>) -> Self {
        assert_eq!(value.unique_match_cost(), Some(Cost::zero()));
        Self {
            substitution: value.unique_substitution_cost().unwrap(),
            gap_open: value.unique_gap_open_cost().unwrap(),
            gap_extend: value.unique_gap_extend_cost().unwrap(),
        }
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost> From<GapAffineAlignmentCostTable<AlphabetType, Cost>>
    for GapAffineCosts<Cost>
{
    fn from(value: GapAffineAlignmentCostTable<AlphabetType, Cost>) -> Self {
        (&value).into()
    }
}
