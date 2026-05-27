use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;
use num_traits::{Bounded, bounds::UpperBounded};

use crate::{
    a_star_aligner::template_switch_distance::{
        TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection,
    },
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    error::{Error, Result},
};

pub mod io;

#[derive(Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        deserialize = "Cost: serde::Deserialize<'de>", // omit AlphabetType
        serialize = "Cost: serde::Serialize" // omit AlphabetType
    ))
)]
pub struct TemplateSwitchConfig<AlphabetType, Cost> {
    // Limits
    pub left_flank_length: isize,
    pub right_flank_length: isize,
    /// The minimum length of a template switch is the first length for which
    /// the length costs are finite.
    pub template_switch_min_length: usize,

    // Base cost
    pub base_cost: BaseCost<Cost>,

    // Edit costs
    pub primary_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_forward_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_reverse_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub left_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub right_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,

    // Jump costs
    pub rq_qr_offset_costs: CostFunction<isize, Cost>,
    pub rr_qq_offset_costs: CostFunction<isize, Cost>,
    pub length_costs: CostFunction<usize, Cost>,
    pub length_difference_costs: CostFunction<isize, Cost>,
    pub forward_anti_descendant_gap_costs: CostFunction<isize, Cost>,
    pub reverse_anti_descendant_gap_costs: CostFunction<isize, Cost>,
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCost<Cost> {
    /// Descendant: reference; ancestor: reference; direction: forward.
    pub rrf: Cost,
    /// Descendant: reference; ancestor: query; direction: forward.
    pub rqf: Cost,
    /// Descendant: query; ancestor: reference; direction: forward.
    pub qrf: Cost,
    /// Descendant: query; ancestor: query; direction: forward.
    pub qqf: Cost,
    /// Descendant: reference; ancestor: reference; direction: reverse.
    pub rrr: Cost,
    /// Descendant: reference; ancestor: query; direction: reverse.
    pub rqr: Cost,
    /// Descendant: query; ancestor: reference; direction: reverse.
    pub qrr: Cost,
    /// Descendant: query; ancestor: query; direction: reverse.
    pub qqr: Cost,
}

impl<AlphabetType, Cost: Bounded + Ord> TemplateSwitchConfig<AlphabetType, Cost> {
    /// Returns an error if any cost function is malformed.
    pub fn verify(&self) -> Result<()> {
        if !self.rq_qr_offset_costs.is_v_shaped() {
            Err(Error::RQQROffsetCostsNotVShaped)
        } else if !self.rr_qq_offset_costs.is_v_shaped() {
            Err(Error::RRQQOffsetCostsNotVShaped)
        } else if !self.length_difference_costs.is_v_shaped() {
            Err(Error::LengthDifferenceCostsNotVShaped)
        } else {
            Ok(())
        }
    }
}

impl<AlphabetType, Cost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn secondary_edit_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &GapAffineAlignmentCostTable<AlphabetType, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.secondary_forward_edit_costs,
            TemplateSwitchDirection::Reverse => &self.secondary_reverse_edit_costs,
        }
    }

    pub fn anti_descendant_gap_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &CostFunction<isize, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.forward_anti_descendant_gap_costs,
            TemplateSwitchDirection::Reverse => &self.reverse_anti_descendant_gap_costs,
        }
    }

    pub fn offset_costs(
        &self,
        descendant: TemplateSwitchDescendant,
        ancestor: TemplateSwitchAncestor,
    ) -> &CostFunction<isize, Cost> {
        match (descendant, ancestor) {
            (TemplateSwitchDescendant::Reference, TemplateSwitchAncestor::Reference) => {
                &self.rr_qq_offset_costs
            }
            (TemplateSwitchDescendant::Reference, TemplateSwitchAncestor::Query) => {
                &self.rq_qr_offset_costs
            }
            (TemplateSwitchDescendant::Query, TemplateSwitchAncestor::Reference) => {
                &self.rq_qr_offset_costs
            }
            (TemplateSwitchDescendant::Query, TemplateSwitchAncestor::Query) => {
                &self.rr_qq_offset_costs
            }
        }
    }
}

impl<Cost: UpperBounded> BaseCost<Cost> {
    pub fn new_max() -> Self {
        Self {
            rrf: Cost::max_value(),
            rqf: Cost::max_value(),
            qrf: Cost::max_value(),
            qqf: Cost::max_value(),
            rrr: Cost::max_value(),
            rqr: Cost::max_value(),
            qrr: Cost::max_value(),
            qqr: Cost::max_value(),
        }
    }
}

impl<Cost: Clone> BaseCost<Cost> {
    pub fn get(
        &self,
        descendant: TemplateSwitchDescendant,
        ancestor: TemplateSwitchAncestor,
        direction: TemplateSwitchDirection,
    ) -> Cost {
        match (descendant, ancestor, direction) {
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Forward,
            ) => self.rrf.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Reverse,
            ) => self.rrr.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Forward,
            ) => self.rqf.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Reverse,
            ) => self.rqr.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Forward,
            ) => self.qrf.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Reverse,
            ) => self.qrr.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Forward,
            ) => self.qqf.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Reverse,
            ) => self.qqr.clone(),
        }
    }
}

impl<AlphabetType: Alphabet, Cost: Clone> Clone for TemplateSwitchConfig<AlphabetType, Cost> {
    fn clone(&self) -> Self {
        Self {
            left_flank_length: self.left_flank_length,
            right_flank_length: self.right_flank_length,
            template_switch_min_length: self.template_switch_min_length,
            base_cost: self.base_cost.clone(),
            primary_edit_costs: self.primary_edit_costs.clone(),
            secondary_forward_edit_costs: self.secondary_forward_edit_costs.clone(),
            secondary_reverse_edit_costs: self.secondary_reverse_edit_costs.clone(),
            left_flank_edit_costs: self.left_flank_edit_costs.clone(),
            right_flank_edit_costs: self.right_flank_edit_costs.clone(),
            rq_qr_offset_costs: self.rq_qr_offset_costs.clone(),
            rr_qq_offset_costs: self.rr_qq_offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
            forward_anti_descendant_gap_costs: self.forward_anti_descendant_gap_costs.clone(),
            reverse_anti_descendant_gap_costs: self.reverse_anti_descendant_gap_costs.clone(),
        }
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost> Default for TemplateSwitchConfig<AlphabetType, Cost> {
    fn default() -> Self {
        Self {
            left_flank_length: 0,
            right_flank_length: 0,
            template_switch_min_length: 5,
            base_cost: BaseCost {
                rrf: 4.into(),
                rqf: 4.into(),
                qrf: 4.into(),
                qqf: 4.into(),
                rrr: 3.into(),
                rqr: 2.into(),
                qrr: 2.into(),
                qqr: 3.into(),
            },
            primary_edit_costs: GapAffineAlignmentCostTable::new_base_agnostic(
                "Primary Edit Costs",
                0.into(),
                2.into(),
                3.into(),
                1.into(),
            ),
            secondary_forward_edit_costs: GapAffineAlignmentCostTable::new_base_agnostic(
                "Secondary Forward Edit Costs",
                0.into(),
                2.into(),
                3.into(),
                1.into(),
            ),
            secondary_reverse_edit_costs: GapAffineAlignmentCostTable::new_base_agnostic(
                "Secondary Reverse Edit Costs",
                0.into(),
                2.into(),
                3.into(),
                1.into(),
            ),
            left_flank_edit_costs: GapAffineAlignmentCostTable::new_base_agnostic(
                "Left Flank Edit Costs",
                0.into(),
                2.into(),
                3.into(),
                1.into(),
            ),
            right_flank_edit_costs: GapAffineAlignmentCostTable::new_base_agnostic(
                "Right Flank Edit Costs",
                0.into(),
                2.into(),
                3.into(),
                1.into(),
            ),
            rq_qr_offset_costs: CostFunction::try_from(vec![
                (isize::MIN, Cost::max_value()),
                (-100, 0.into()),
                (101, Cost::max_value()),
            ])
            .unwrap(),
            rr_qq_offset_costs: CostFunction::try_from(vec![
                (isize::MIN, Cost::max_value()),
                (-100, 0.into()),
                (1, Cost::max_value()),
            ])
            .unwrap(),
            length_costs: CostFunction::try_from(vec![(0, Cost::max_value()), (5, 0.into())])
                .unwrap(),
            length_difference_costs: CostFunction::try_from(vec![
                (isize::MIN, Cost::max_value()),
                (-100, 0.into()),
                (101, Cost::max_value()),
            ])
            .unwrap(),
            forward_anti_descendant_gap_costs: CostFunction::try_from(vec![
                (isize::MIN, Cost::max_value()),
                (-100, 0.into()),
                (101, Cost::max_value()),
            ])
            .unwrap(),
            reverse_anti_descendant_gap_costs: CostFunction::try_from(vec![
                (isize::MIN, Cost::max_value()),
                (-100, 0.into()),
                (101, Cost::max_value()),
            ])
            .unwrap(),
        }
    }
}
