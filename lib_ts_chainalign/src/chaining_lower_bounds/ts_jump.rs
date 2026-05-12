use std::io::{Read, Write};

use generic_a_star::cost::AStarCost;

use crate::{
    chaining_lower_bounds::{cost_array::LowerBoundCostArray, gap_affine::GapAffineLowerBounds},
    costs::AlignmentCosts,
};

#[cfg(test)]
mod tests;

pub struct TsJumpLowerBounds<Cost> {
    lower_bounds_12: LowerBoundCostArray<1, Cost>,
    lower_bounds_34: LowerBoundCostArray<1, Cost>,
}

impl<Cost: AStarCost> TsJumpLowerBounds<Cost> {
    pub fn new(max_n: usize, max_match_run: u32, cost_table: &AlignmentCosts<Cost>) -> Self {
        let primary_lower_bounds = GapAffineLowerBounds::new_allow_all_matches(
            max_n,
            max_match_run,
            &cost_table.primary_costs,
        );
        let secondary_lower_bounds = GapAffineLowerBounds::new_allow_all_matches(
            max_n,
            max_match_run,
            &cost_table.secondary_costs,
        );

        // This way of calculating the lower bound for the 12-jump does not take the shape limits of the template switch into account.
        // However, most of the time these limits are gonna be big, so they should not have a big impact on the lower bound.
        let mut lower_bounds_12 =
            LowerBoundCostArray::new_from_cost([max_n + 1], Cost::max_value());
        for primary_descendant_gap in 0..=max_n {
            for secondary_descendant_gap in 0..=max_n - primary_descendant_gap {
                let lower_bound = primary_lower_bounds
                    .variable_gap2_lower_bound(primary_descendant_gap)
                    + cost_table.ts_base_cost.min()
                    + secondary_lower_bounds.variable_gap2_lower_bound(secondary_descendant_gap);
                lower_bounds_12[[primary_descendant_gap + secondary_descendant_gap]] =
                    lower_bounds_12[[primary_descendant_gap + secondary_descendant_gap]]
                        .min(lower_bound);
            }
        }

        let mut lower_bounds_34 =
            LowerBoundCostArray::new_from_cost([max_n + 1], Cost::max_value());
        for secondary_descendant_gap in 0..=max_n {
            for primary_descendant_gap in 0..=max_n - secondary_descendant_gap {
                let lower_bound = secondary_lower_bounds
                    .variable_gap2_lower_bound(secondary_descendant_gap)
                    + primary_lower_bounds.variable_gap2_lower_bound(primary_descendant_gap);
                lower_bounds_34[[primary_descendant_gap + secondary_descendant_gap]] =
                    lower_bounds_34[[primary_descendant_gap + secondary_descendant_gap]]
                        .min(lower_bound);
            }
        }

        Self {
            lower_bounds_12,
            lower_bounds_34,
        }
    }

    pub fn write(&self, mut write: impl Write) -> std::io::Result<()>
    where
        Cost: Copy,
    {
        self.lower_bounds_12.write(&mut write)?;
        self.lower_bounds_34.write(write)
    }

    pub fn read(mut read: impl Read) -> std::io::Result<Self>
    where
        Cost: Copy,
    {
        let lower_bounds_12 = LowerBoundCostArray::read(&mut read)?;
        let lower_bounds_34 = LowerBoundCostArray::read(read)?;
        Ok(Self {
            lower_bounds_12,
            lower_bounds_34,
        })
    }
}

impl<Cost: Copy> TsJumpLowerBounds<Cost> {
    /// A lower bound of the cost for chaining a primary anchor with a secondary anchor.
    ///
    /// This lower bound takes the template switch base cost into account.
    pub fn lower_bound_12(&self, descendant_gap: usize) -> Cost {
        self.lower_bounds_12[[descendant_gap]]
    }

    /// A lower bound of the cost for chaining a secondary anchor with a primary anchor.
    ///
    /// This lower bound does **not** take the template switch base cost into account.
    pub fn lower_bound_34(&self, descendant_gap: usize) -> Cost {
        self.lower_bounds_34[[descendant_gap]]
    }
}
