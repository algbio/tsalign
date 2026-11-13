//! Compute lower bounds for chaining anchors with gaps.

use std::io::{Read, Write};

use generic_a_star::cost::AStarCost;
use serde::{Serialize, de::DeserializeOwned};

use crate::{
    chaining_lower_bounds::{gap_affine::GapAffineLowerBounds, ts_jump::TsJumpLowerBounds},
    config::ChainingLowerBoundConfig,
    costs::AlignmentCosts,
};

mod cost_array;
pub mod gap_affine;
#[cfg(test)]
mod tests;
pub mod ts_jump;

pub struct ChainingLowerBounds<Cost> {
    primary: GapAffineLowerBounds<Cost>,
    secondary: GapAffineLowerBounds<Cost>,
    jump: TsJumpLowerBounds<Cost>,
    alignment_costs: AlignmentCosts<Cost>,
    config: ChainingLowerBoundConfig,
}

impl<Cost: AStarCost> ChainingLowerBounds<Cost> {
    /// Compute chaining lower bounds.
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `anchor_k` is the size of the anchors.
    /// * `max_anchor_mutations` is the maximum number of mutations allowed in an anchor.
    /// * `alignment_costs` is the cost function for the alignment.
    pub fn new(
        max_n: usize,
        config: ChainingLowerBoundConfig,
        alignment_costs: AlignmentCosts<Cost>,
    ) -> Self {
        if config.max_anchor_mutations == 0 {
            Self {
                primary: GapAffineLowerBounds::new_exact_anchors(
                    max_n,
                    config.anchor_k - 1,
                    &alignment_costs.primary_costs,
                ),
                secondary: GapAffineLowerBounds::new_exact_anchors(
                    max_n,
                    config.anchor_k - 1,
                    &alignment_costs.secondary_costs,
                ),
                jump: TsJumpLowerBounds::new_exact(max_n, config.anchor_k - 1, &alignment_costs),
                alignment_costs,
                config,
            }
        } else {
            let primary = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                &config,
                &alignment_costs.primary_costs,
            );
            let secondary = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                &config,
                &alignment_costs.secondary_costs,
            );

            Self {
                jump: TsJumpLowerBounds::new_inexact(max_n, &alignment_costs, &primary, &secondary),
                primary,
                secondary,

                alignment_costs,
                config,
            }
        }
    }

    pub fn write(&self, mut write: impl Write) -> std::io::Result<()>
    where
        Cost: Copy + Serialize,
    {
        self.primary.write(&mut write)?;
        self.secondary.write(&mut write)?;
        self.jump.write(&mut write)?;
        bincode::serde::encode_into_std_write(
            &self.alignment_costs,
            &mut write,
            bincode::config::standard(),
        )
        .map_err(|error| match error {
            bincode::error::EncodeError::Io { inner, .. } => inner,
            error => panic!("I/O error: {error}"),
        })?;
        bincode::serde::encode_into_std_write(
            &self.config,
            &mut write,
            bincode::config::standard(),
        )
        .map_err(|error| match error {
            bincode::error::EncodeError::Io { inner, .. } => inner,
            error => panic!("I/O error: {error}"),
        })?;

        Ok(())
    }

    pub fn read(mut read: impl Read) -> std::io::Result<Self>
    where
        Cost: Copy + DeserializeOwned,
    {
        let primary = GapAffineLowerBounds::read(&mut read)?;
        let secondary = GapAffineLowerBounds::read(&mut read)?;
        let jump = TsJumpLowerBounds::read(&mut read)?;
        let alignment_costs =
            bincode::serde::decode_from_std_read(&mut read, bincode::config::standard()).map_err(
                |error| match error {
                    bincode::error::DecodeError::Io { inner, .. } => inner,
                    error => panic!("I/O error: {error}"),
                },
            )?;
        let config = bincode::serde::decode_from_std_read(&mut read, bincode::config::standard())
            .map_err(|error| match error {
            bincode::error::DecodeError::Io { inner, .. } => inner,
            error => panic!("I/O error: {error}"),
        })?;

        Ok(Self {
            primary,
            secondary,
            jump,
            alignment_costs,
            config,
        })
    }
}

impl<Cost: Copy> ChainingLowerBounds<Cost> {
    pub fn primary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.primary.lower_bound(gap1, gap2)
    }

    pub fn secondary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.secondary.lower_bound(gap1, gap2)
    }

    pub fn jump_12_lower_bound(&self, descendant_gap: usize) -> Cost {
        self.jump.lower_bound_12(descendant_gap)
    }

    pub fn jump_34_lower_bound(&self, descendant_gap: usize) -> Cost {
        self.jump.lower_bound_34(descendant_gap)
    }
}

/// Lower bounds for gaps that start at the start of the alignment or end at the end of the alignment.
///
/// For exact anchors, the lower bound tables assume that a gap starts and ends with a non-match,
/// because a match right next to an anchor implies that the anchor shifted by one also exists,
/// such that the chain can use the shifted anchor instead.
/// The start and the end of the alignment are no anchors, so this assumption does not hold for gaps that touch them.
/// Moreover, the aligner that computes the exact costs of gaps during chaining does not enforce non-matches next to anchors either,
/// so the bounds for gaps that touch the start or end make no assumption about either end of the gap.
///
/// Since runs of matches are limited to `k - 1`, such a gap consists of at most `k - 1` matches at its start,
/// at most `k - 1` matches at its end, and a remainder that starts and ends with a non-match, which is bounded by the tables.
///
/// For inexact anchors, the lower bound tables make no assumptions about the start and end of a gap,
/// so these methods return the same values as the methods for gaps between anchors.
impl<Cost: AStarCost> ChainingLowerBounds<Cost> {
    /// The maximum number of matches at the start or end of a gap that touches the start or end of the alignment,
    /// or `None` if the lower bound tables do not need to be adjusted for these.
    fn max_boundary_match_run(&self) -> Option<usize> {
        (self.config.max_anchor_mutations == 0)
            .then(|| usize::try_from(self.config.anchor_k - 1).unwrap())
    }

    /// A lower bound for a primary gap that may start and end with up to `max_match_run` matches.
    ///
    /// Gaps that consist of matches only are not covered, and need to be handled by the caller.
    fn primary_lower_bound_with_boundary_matches(
        &self,
        gap1: usize,
        gap2: usize,
        max_match_run: usize,
    ) -> Cost {
        let mut result = Cost::max_value();
        for start_matches in 0..=gap1.min(gap2).min(max_match_run) {
            let remaining1 = gap1 - start_matches;
            let remaining2 = gap2 - start_matches;
            for end_matches in 0..=remaining1.min(remaining2).min(max_match_run) {
                if remaining1 == end_matches && remaining2 == end_matches {
                    // The gap consists of matches only.
                    continue;
                }
                result = result.min(
                    self.primary
                        .lower_bound(remaining1 - end_matches, remaining2 - end_matches),
                );
            }
        }
        result
    }

    /// A lower bound for a primary gap from the start of the alignment to an anchor.
    pub fn primary_lower_bound_from_start(&self, gap1: usize, gap2: usize) -> Cost {
        let Some(max_match_run) = self.max_boundary_match_run() else {
            return self.primary_lower_bound(gap1, gap2);
        };
        if gap1 == 0 && gap2 == 0 {
            // Direct chaining from the start.
            return Cost::zero();
        }

        self.primary_lower_bound_with_boundary_matches(gap1, gap2, max_match_run)
    }

    /// A lower bound for a primary gap from an anchor to the end of the alignment.
    pub fn primary_lower_bound_to_end(&self, gap1: usize, gap2: usize) -> Cost {
        // Gaps to the end are subject to the same conditions as gaps from the start.
        self.primary_lower_bound_from_start(gap1, gap2)
    }

    /// A lower bound for a primary gap from the start to the end of the alignment.
    pub fn primary_lower_bound_start_to_end(&self, gap1: usize, gap2: usize) -> Cost {
        let Some(max_match_run) = self.max_boundary_match_run() else {
            return self.primary_lower_bound(gap1, gap2);
        };
        if gap1 == gap2 && gap1 <= max_match_run {
            // The alignment may consist of matches only.
            return Cost::zero();
        }

        self.primary_lower_bound_with_boundary_matches(gap1, gap2, max_match_run)
    }

    /// A lower bound for a 12-jump from the start of the alignment to a secondary anchor.
    pub fn jump_12_lower_bound_from_start(&self, descendant_gap: usize) -> Cost {
        let Some(max_match_run) = self.max_boundary_match_run() else {
            return self.jump_12_lower_bound(descendant_gap);
        };

        (0..=descendant_gap.min(max_match_run))
            .map(|matches| self.jump.lower_bound_12(descendant_gap - matches))
            .min()
            .unwrap()
    }

    /// A lower bound for a 34-jump from a secondary anchor to the end of the alignment.
    pub fn jump_34_lower_bound_to_end(&self, descendant_gap: usize) -> Cost {
        let Some(max_match_run) = self.max_boundary_match_run() else {
            return self.jump_34_lower_bound(descendant_gap);
        };

        (0..=descendant_gap.min(max_match_run))
            .map(|matches| self.jump.lower_bound_34(descendant_gap - matches))
            .min()
            .unwrap()
    }
}

impl<Cost> ChainingLowerBounds<Cost> {
    pub fn primary(&self) -> &GapAffineLowerBounds<Cost> {
        &self.primary
    }

    pub fn secondary(&self) -> &GapAffineLowerBounds<Cost> {
        &self.secondary
    }

    pub fn jump(&self) -> &TsJumpLowerBounds<Cost> {
        &self.jump
    }

    pub fn alignment_costs(&self) -> &AlignmentCosts<Cost> {
        &self.alignment_costs
    }

    pub fn config(&self) -> &ChainingLowerBoundConfig {
        &self.config
    }

    pub fn anchor_k(&self) -> u32 {
        self.config.anchor_k
    }

    pub fn max_anchor_mutations(&self) -> u8 {
        self.config.max_anchor_mutations
    }
}
