//! Compute lower bounds for chaining anchors with gaps.

use std::io::{Read, Write};

use generic_a_star::cost::AStarCost;
use serde::{Serialize, de::DeserializeOwned};

use crate::{
    chaining_lower_bounds::{gap_affine::GapAffineLowerBounds, ts_jump::TsJumpLowerBounds},
    costs::AlignmentCosts,
};

mod cost_array;
pub mod gap_affine;
pub mod ts_jump;

pub struct ChainingLowerBounds<Cost> {
    primary: GapAffineLowerBounds<Cost>,
    secondary: GapAffineLowerBounds<Cost>,
    jump: TsJumpLowerBounds<Cost>,
    alignment_costs: AlignmentCosts<Cost>,
    anchor_k: u32,
    max_anchor_mutations: u8,
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
        anchor_k: u32,
        max_anchor_mutations: u8,
        alignment_costs: AlignmentCosts<Cost>,
    ) -> Self {
        if max_anchor_mutations == 0 {
            Self {
                primary: GapAffineLowerBounds::new_exact_anchors(
                    max_n,
                    anchor_k - 1,
                    &alignment_costs.primary_costs,
                ),
                secondary: GapAffineLowerBounds::new_exact_anchors(
                    max_n,
                    anchor_k - 1,
                    &alignment_costs.secondary_costs,
                ),
                jump: TsJumpLowerBounds::new_exact(max_n, anchor_k - 1, &alignment_costs),
                alignment_costs,
                anchor_k,
                max_anchor_mutations,
            }
        } else {
            let anchor_k = u8::try_from(anchor_k).expect("anchor_k must be <= 255.");
            let primary = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                anchor_k,
                max_anchor_mutations,
                &alignment_costs.primary_costs,
            );
            let secondary = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                anchor_k,
                max_anchor_mutations,
                &alignment_costs.secondary_costs,
            );

            Self {
                jump: TsJumpLowerBounds::new_inexact(max_n, &alignment_costs, &primary, &secondary),
                primary,
                secondary,

                alignment_costs,
                anchor_k: anchor_k.into(),
                max_anchor_mutations,
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
        write.write_all(&self.anchor_k.to_ne_bytes())?;
        write.write_all(&self.max_anchor_mutations.to_ne_bytes())
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

        let mut buffer = [0; std::mem::size_of::<u32>()];
        read.read_exact(&mut buffer)?;
        let anchor_k = u32::from_ne_bytes(buffer);

        let mut buffer = [0; std::mem::size_of::<u8>()];
        read.read_exact(&mut buffer)?;
        let max_anchor_mutations = u8::from_ne_bytes(buffer);

        Ok(Self {
            primary,
            secondary,
            jump,
            alignment_costs,
            anchor_k,
            max_anchor_mutations,
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

    pub fn anchor_k(&self) -> u32 {
        self.anchor_k
    }

    pub fn max_anchor_mutations(&self) -> u8 {
        self.max_anchor_mutations
    }
}
