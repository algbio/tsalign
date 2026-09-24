use std::io::{Read, Write};

use generic_a_star::{AStar, AStarNode, cost::AStarCost};

use crate::{
    alignment::coordinates::PrimaryAlignmentCoordinates,
    chaining_lower_bounds::{
        cost_array::LowerBoundCostArray,
        gap_affine::{
            exact_algo::Context,
            inexact_algo::{AlignmentHistory, UnsignedIntAlignmentHistoryVec},
        },
    },
    costs::GapAffineCosts,
};

mod exact_algo;
mod inexact_algo;
#[cfg(test)]
mod tests;

pub struct GapAffineLowerBounds<Cost> {
    lower_bounds: LowerBoundCostArray<2, Cost>,
    variable_gap2_lower_bounds: LowerBoundCostArray<1, Cost>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[expect(clippy::enum_variant_names)]
enum BoundaryCondition {
    BothNonMatchAndNoDirectChaining,
    StartNonMatchAndDirectChaining,
    EndNonMatchAndDirectChaining,
}

impl<Cost: AStarCost> GapAffineLowerBounds<Cost> {
    /// Compute the lower bounds for the case that the anchors are exact matches.
    ///
    /// Disallow direct chaining between anchors and enforce that the alignment starts and ends with non-matches.
    ///
    /// # Parameters
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
    ///   Set this to `k-1`, if the anchors are `k`-mers.
    /// * `cost_table` is the cost function for the alignment.
    pub fn new_exact_anchors(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        Self::compute_exact(
            max_n,
            max_match_run,
            cost_table,
            BoundaryCondition::BothNonMatchAndNoDirectChaining,
        )
    }

    /// Compute the lower bounds for the case that the anchors are exact matches.
    ///
    /// Allow direct chaining between anchors and enforce that the alignment starts with a non-match.
    ///
    /// This is useful for computing the jump lower bounds before the jump, as the jump itself is a non-match operation,
    /// so there is no need to enforce non-matches right before the jump.
    ///
    /// # Parameters
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
    ///   Set this to `k-1`, if the anchors are `k`-mers.
    /// * `cost_table` is the cost function for the alignment.
    pub(super) fn new_exact_anchors_allow_end_match(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        Self::compute_exact(
            max_n,
            max_match_run,
            cost_table,
            BoundaryCondition::StartNonMatchAndDirectChaining,
        )
    }

    /// Compute the lower bounds for the case that the anchors are exact matches.
    ///
    /// Allow direct chaining between anchors and enforce that the alignment ends with a non-match.
    ///
    /// This is useful for computing the jump lower bounds after the jump, as the jump itself is a non-match operation,
    /// so there is no need to enforce non-matches right after the jump.
    ///
    /// # Parameters
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
    ///   Set this to `k-1`, if the anchors are `k`-mers.
    /// * `cost_table` is the cost function for the alignment.
    pub(super) fn new_exact_anchors_allow_start_match(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        Self::compute_exact(
            max_n,
            max_match_run,
            cost_table,
            BoundaryCondition::EndNonMatchAndDirectChaining,
        )
    }

    /// Compute the lower bounds for the case that the anchors are inexact matches.
    ///
    /// Inexact anchors don't have any nice overlap-extension property where a match right after or right before an anchor `A` means that there is another anchor `B` that can be chained with `A`.
    /// Therefore, we cannot assume that a gap starts or ends with a non-match, or even contains a non-match at all if the gap is at most `max_match_run` long.
    /// Hence, the lower bound is computed only with the condition that no subalignment of the alignment can be an anchor,
    /// so every subalignment where one of the sequences has length `k` must have at least `max_anchor_mutations + 1` mutations.
    ///
    /// # Parameters
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `anchor_k` is the length of the anchors.
    /// * `max_anchor_mutations` is the maximum number of mutations that an anchor can have.
    /// * `cost_table` is the cost function for the alignment.
    pub fn new_inexact_anchors(
        max_n: usize,
        anchor_k: u8,
        max_anchor_mutations: u8,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        if UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(anchor_k, max_anchor_mutations)
        {
            Self::new_inexact_anchors_with_history_type::<UnsignedIntAlignmentHistoryVec<u8>>(
                max_n,
                anchor_k,
                max_anchor_mutations,
                cost_table,
            )
        } else if UnsignedIntAlignmentHistoryVec::<u16>::has_enough_capacity(
            anchor_k,
            max_anchor_mutations,
        ) {
            Self::new_inexact_anchors_with_history_type::<UnsignedIntAlignmentHistoryVec<u16>>(
                max_n,
                anchor_k,
                max_anchor_mutations,
                cost_table,
            )
        } else if UnsignedIntAlignmentHistoryVec::<u32>::has_enough_capacity(
            anchor_k,
            max_anchor_mutations,
        ) {
            Self::new_inexact_anchors_with_history_type::<UnsignedIntAlignmentHistoryVec<u32>>(
                max_n,
                anchor_k,
                max_anchor_mutations,
                cost_table,
            )
        } else if UnsignedIntAlignmentHistoryVec::<u64>::has_enough_capacity(
            anchor_k,
            max_anchor_mutations,
        ) {
            Self::new_inexact_anchors_with_history_type::<UnsignedIntAlignmentHistoryVec<u64>>(
                max_n,
                anchor_k,
                max_anchor_mutations,
                cost_table,
            )
        } else if UnsignedIntAlignmentHistoryVec::<u128>::has_enough_capacity(
            anchor_k,
            max_anchor_mutations,
        ) {
            Self::new_inexact_anchors_with_history_type::<UnsignedIntAlignmentHistoryVec<u128>>(
                max_n,
                anchor_k,
                max_anchor_mutations,
                cost_table,
            )
        } else {
            panic!(
                "This combination of k and max_mismatches exceeds the maximum length of the alignment history vector that can be stored in a u128.",
            );
        }
    }

    fn new_inexact_anchors_with_history_type<AlignmentHistoryVec: AlignmentHistory>(
        max_n: usize,
        anchor_k: u8,
        max_anchor_mutations: u8,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        let mut lower_bounds =
            LowerBoundCostArray::new_from_cost([max_n + 1, max_n + 1], Cost::max_value());

        let context = inexact_algo::Context::<_, AlignmentHistoryVec>::new(
            cost_table,
            anchor_k,
            max_anchor_mutations,
            max_n,
        );
        let mut a_star = AStar::<_>::new(context);
        a_star.initialise();
        a_star.search_until(|_, node| {
            let lower_bound = &mut lower_bounds[[
                node.identifier.coordinates.primary_ordinate_a().unwrap(),
                node.identifier.coordinates.primary_ordinate_b().unwrap(),
            ]];
            *lower_bound = (*lower_bound).min(node.cost());

            false
        });
        let variable_gap2_lower_bounds = LowerBoundCostArray::from_iter((0..=max_n).map(|gap1| {
            (0..=max_n)
                .map(|gap2| lower_bounds[[gap1, gap2]])
                .min()
                .unwrap()
        }));

        Self {
            lower_bounds,
            variable_gap2_lower_bounds,
        }
    }

    /// Compute the lower bounds for the case that the anchors are exact matches.
    ///
    /// # Parameters
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
    ///   Set this to `k-1`, if the anchors are `k`-mers.
    /// * `cost_table` is the cost function for the alignment.
    /// * `boundary_condition` determines if there must be non-matches at the start or end of the alignment, and if direct chaining (i.e. a length-zero alignment) between anchors is allowed.
    fn compute_exact(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
        boundary_condition: BoundaryCondition,
    ) -> Self {
        let mut lower_bounds =
            LowerBoundCostArray::new_from_cost([max_n + 1, max_n + 1], Cost::max_value());
        if boundary_condition.allow_direct_chaining() {
            lower_bounds[[0, 0]] = Cost::zero();
        }

        let context = Context::new(
            cost_table,
            max_match_run,
            max_n,
            boundary_condition.allow_start_match(),
        );
        let mut a_star = AStar::<_>::new(context);
        a_star.initialise();
        a_star.search_until(|_, node| {
            if node.identifier.match_run == 0
                || (boundary_condition.allow_end_match()
                    && node.identifier.coordinates.into_primary().unwrap()
                        != PrimaryAlignmentCoordinates::new(0, 0))
            {
                let lower_bound = &mut lower_bounds[[
                    node.identifier.coordinates.primary_ordinate_a().unwrap(),
                    node.identifier.coordinates.primary_ordinate_b().unwrap(),
                ]];
                *lower_bound = (*lower_bound).min(node.cost());
            }
            false
        });
        let variable_gap2_lower_bounds = LowerBoundCostArray::from_iter((0..=max_n).map(|gap1| {
            (0..=max_n)
                .map(|gap2| lower_bounds[[gap1, gap2]])
                .min()
                .unwrap()
        }));

        Self {
            lower_bounds,
            variable_gap2_lower_bounds,
        }
    }

    pub fn write(&self, mut write: impl Write) -> std::io::Result<()>
    where
        Cost: Copy,
    {
        self.lower_bounds.write(&mut write)?;
        self.variable_gap2_lower_bounds.write(write)
    }

    pub fn read(mut read: impl Read) -> std::io::Result<Self>
    where
        Cost: Copy,
    {
        let lower_bounds = LowerBoundCostArray::read(&mut read)?;
        let variable_gap2_lower_bounds = LowerBoundCostArray::read(read)?;
        Ok(Self {
            lower_bounds,
            variable_gap2_lower_bounds,
        })
    }
}

impl<Cost: Copy> GapAffineLowerBounds<Cost> {
    /// A lower bound of the cost for chaining two anchors with the given gaps.
    /// The lower bound is symmetric, so the order of the gaps does not matter.
    pub fn lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.lower_bounds[[gap1, gap2]]
    }

    /// A lower bound of the cost for chaining two anchors with only one specified gap length.
    pub fn variable_gap2_lower_bound(&self, gap: usize) -> Cost {
        self.variable_gap2_lower_bounds[[gap]]
    }
}

impl BoundaryCondition {
    fn allow_direct_chaining(&self) -> bool {
        match self {
            BoundaryCondition::BothNonMatchAndNoDirectChaining => false,
            BoundaryCondition::StartNonMatchAndDirectChaining => true,
            BoundaryCondition::EndNonMatchAndDirectChaining => true,
        }
    }

    fn allow_start_match(&self) -> bool {
        match self {
            BoundaryCondition::BothNonMatchAndNoDirectChaining => false,
            BoundaryCondition::StartNonMatchAndDirectChaining => false,
            BoundaryCondition::EndNonMatchAndDirectChaining => true,
        }
    }

    fn allow_end_match(&self) -> bool {
        match self {
            BoundaryCondition::BothNonMatchAndNoDirectChaining => false,
            BoundaryCondition::StartNonMatchAndDirectChaining => true,
            BoundaryCondition::EndNonMatchAndDirectChaining => false,
        }
    }
}
