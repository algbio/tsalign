//! Functions for making tschainalign compatible with tsalign types.

use compact_genome::{
    implementation::alphabets::dna_alphabet::DnaAlphabet,
    interface::alphabet::{Alphabet, AlphabetError},
};
use generic_a_star::cost::AStarCost;
use lib_tsalign::{
    aligner::TemplateSwitchDirection, config::TemplateSwitchConfig,
    costs::gap_affine::GapAffineAlignmentCostTable,
};
use thiserror::Error;

use crate::{
    alignment::ts_kind::TsKind,
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
};

/// A part of a tsalign cost configuration that cannot be converted without losing information.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum CostConversionError {
    #[error("the alphabet must contain only characters A, C, G and T: {0}")]
    UnsupportedAlphabet(AlphabetError),
    #[error("flank lengths must be zero")]
    NonZeroFlankLength,
    #[error("all match costs must be zero")]
    NonZeroMatchCost,
    #[error("all substitution costs must be equal")]
    NonUniqueSubstitutionCost,
    #[error("all gap open costs must be equal")]
    NonUniqueGapOpenCost,
    #[error("all gap extend costs must be equal")]
    NonUniqueGapExtendCost,
    #[error("the {0} costs must be zero within a single range and infinite outside of it")]
    NoZeroRange(&'static str),
}

impl<AlphabetType: Alphabet, Cost: AStarCost> TryFrom<TemplateSwitchConfig<AlphabetType, Cost>>
    for AlignmentCosts<Cost>
{
    type Error = CostConversionError;

    fn try_from(value: TemplateSwitchConfig<AlphabetType, Cost>) -> Result<Self, Self::Error> {
        if value.left_flank_length != 0 || value.right_flank_length != 0 {
            return Err(CostConversionError::NonZeroFlankLength);
        }

        let ts_base_cost = TsBaseCost::from_iter([
            (TsKind::TS11, value.base_cost.rrr),
            (TsKind::TS12, value.base_cost.qrr), // Note that the ancestor is the secondary and the descendant is the primary, hence these are flipped.
            (TsKind::TS21, value.base_cost.rqr),
            (TsKind::TS22, value.base_cost.qqr),
        ]);

        Ok(Self {
            primary_costs: restrict_to_dna(&value.primary_edit_costs)?.try_into()?,
            secondary_costs: restrict_to_dna(
                value.secondary_edit_costs(TemplateSwitchDirection::Reverse),
            )?
            .try_into()?,
            ts_base_cost,
            ts_limits: TsLimits {
                inter_jump_12: value
                    .rq_qr_offset_costs
                    .zero_range()
                    .ok_or(CostConversionError::NoZeroRange("RQQROffset"))?,
                intra_jump_12: value
                    .rr_qq_offset_costs
                    .zero_range()
                    .ok_or(CostConversionError::NoZeroRange("RRQQOffset"))?,
                // tsalign costs do not support limiting this.
                jump_34: isize::MIN..isize::MAX,
                length_23: value
                    .length_costs
                    .zero_range()
                    .ok_or(CostConversionError::NoZeroRange("Length"))?,
                ancestor_gap: value
                    .reverse_anti_descendant_gap_costs
                    .zero_range()
                    .ok_or(CostConversionError::NoZeroRange("ReverseAntiDescendantGap"))?,
            },
        })
    }
}

fn restrict_to_dna<AlphabetType: Alphabet, Cost: AStarCost>(
    table: &GapAffineAlignmentCostTable<AlphabetType, Cost>,
) -> Result<GapAffineAlignmentCostTable<DnaAlphabet, Cost>, CostConversionError> {
    table
        .restrict_alphabet()
        .map_err(CostConversionError::UnsupportedAlphabet)
}

impl<Cost: AStarCost> TryFrom<&'_ GapAffineAlignmentCostTable<DnaAlphabet, Cost>>
    for GapAffineCosts<Cost>
{
    type Error = CostConversionError;

    fn try_from(
        value: &GapAffineAlignmentCostTable<DnaAlphabet, Cost>,
    ) -> Result<Self, Self::Error> {
        if value.unique_match_cost() != Some(Cost::zero()) {
            return Err(CostConversionError::NonZeroMatchCost);
        }

        Ok(Self {
            substitution: value
                .unique_substitution_cost()
                .ok_or(CostConversionError::NonUniqueSubstitutionCost)?,
            gap_open: value
                .unique_gap_open_cost()
                .ok_or(CostConversionError::NonUniqueGapOpenCost)?,
            gap_extend: value
                .unique_gap_extend_cost()
                .ok_or(CostConversionError::NonUniqueGapExtendCost)?,
        })
    }
}

impl<Cost: AStarCost> TryFrom<GapAffineAlignmentCostTable<DnaAlphabet, Cost>>
    for GapAffineCosts<Cost>
{
    type Error = CostConversionError;

    fn try_from(
        value: GapAffineAlignmentCostTable<DnaAlphabet, Cost>,
    ) -> Result<Self, Self::Error> {
        (&value).try_into()
    }
}
