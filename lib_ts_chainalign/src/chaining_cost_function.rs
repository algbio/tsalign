use std::time::{Duration, Instant};

use generic_a_star::cost::AStarCost;
use itertools::Itertools;
use log::{debug, info, trace};
use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates, sequences::AlignmentSequences, ts_kind::TsKind,
    },
    anchors::{Anchors, index::AnchorIndex, primary::PrimaryAnchor, secondary::SecondaryAnchor},
    chaining_cost_function::cost_array::ChainingCostArray,
    chaining_lower_bounds::ChainingLowerBounds,
    exact_chaining::{
        gap_affine::GapAffineAligner, ts_12_jump::Ts12JumpAligner, ts_34_jump::Ts34JumpAligner,
    },
    panic_on_extend::PanicOnExtend,
};

mod cost_array;
#[cfg(test)]
mod tests;

pub struct ChainingCostFunction<Cost> {
    primary: ChainingCostArray<Cost>,
    secondaries: [ChainingCostArray<Cost>; 4],
    jump_12s: [ChainingCostArray<Cost>; 4],
    jump_34s: [ChainingCostArray<Cost>; 4],
}

impl<Cost: AStarCost> ChainingCostFunction<Cost> {
    pub fn new_from_lower_bounds(
        chaining_lower_bounds: &ChainingLowerBounds<Cost>,
        anchors: &Anchors,
        sequences: &AlignmentSequences,
        max_exact_cost_function_cost: Cost,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self {
        info!("Initialising chaining cost function...");
        let start_time = Instant::now();

        let start = sequences.primary_start();
        let end = sequences.primary_end();
        let k = usize::try_from(chaining_lower_bounds.max_match_run() + 1).unwrap();
        let primary_anchor_amount = anchors.primary_len() + 2;
        let primary_start_anchor_index = AnchorIndex::zero();
        let primary_end_anchor_index = primary_anchor_amount - 1;

        let mut primary_aligner = GapAffineAligner::new(
            sequences,
            &chaining_lower_bounds.alignment_costs().primary_costs,
            rc_fn,
            chaining_lower_bounds.max_match_run(),
        );
        let mut secondary_aligner = GapAffineAligner::new(
            sequences,
            &chaining_lower_bounds.alignment_costs().secondary_costs,
            rc_fn,
            chaining_lower_bounds.max_match_run(),
        );
        let mut ts_12_jump_aligner = Ts12JumpAligner::new(
            sequences,
            chaining_lower_bounds.alignment_costs(),
            rc_fn,
            chaining_lower_bounds.max_match_run(),
        );
        let mut ts_34_jump_aligner = Ts34JumpAligner::new(
            sequences,
            chaining_lower_bounds.alignment_costs(),
            rc_fn,
            chaining_lower_bounds.max_match_run(),
        );
        let mut additional_primary_targets_output = Vec::new();
        let mut additional_secondary_targets_output = Vec::new();

        // Initialise primary with infinity.
        trace!("Initialise primary with infinity");
        let mut primary = ChainingCostArray::new_from_cost(
            [primary_anchor_amount, primary_anchor_amount],
            Cost::max_value(),
            true,
        );

        // Fill primary from start with exact values.
        trace!("Fill primary");
        additional_primary_targets_output.clear();
        primary_aligner.align_until_cost_limit(
            start,
            max_exact_cost_function_cost,
            &mut additional_primary_targets_output,
            &mut PanicOnExtend,
        );
        // Allow start to chain without gaps to first primary anchor if it exists.
        additional_primary_targets_output
            .push((PrimaryAnchor::new_from_start(&start), Cost::zero()));
        additional_primary_targets_output.sort_unstable();
        for (to_index, cost) in
            anchors.primary_anchor_to_index_iter(additional_primary_targets_output.iter().copied())
        {
            let to_index = to_index + 1;
            primary[[primary_start_anchor_index, to_index]] =
                primary[[primary_start_anchor_index, to_index]].min(cost);
            if cost <= max_exact_cost_function_cost {
                primary.set_exact(primary_start_anchor_index, to_index);
            }
        }
        if let Some(start_end_cost) = additional_primary_targets_output
            .iter()
            .rev()
            .map_while(|(anchor, cost)| (anchor.start() == end).then_some(cost))
            .copied()
            .last()
        {
            primary[[primary_start_anchor_index, primary_end_anchor_index]] = start_end_cost;
            if start_end_cost <= max_exact_cost_function_cost {
                primary.set_exact(primary_start_anchor_index, primary_end_anchor_index);
            }
        }

        // Fill remaining primary with lower bound.
        let gap1 = end.primary_ordinate_a().unwrap() - start.primary_ordinate_a().unwrap();
        let gap2 = end.primary_ordinate_b().unwrap() - start.primary_ordinate_b().unwrap();
        primary[[primary_start_anchor_index, primary_end_anchor_index]] = chaining_lower_bounds
            .primary_lower_bound(gap1, gap2)
            .max(max_exact_cost_function_cost + Cost::from_usize(1))
            .min(primary[[primary_start_anchor_index, primary_end_anchor_index]]);
        for (from_index, from_anchor) in anchors.enumerate_primaries() {
            let from_index = from_index + 1;
            let (gap1, gap2) = from_anchor.chaining_gaps_from_start(start);
            primary[[primary_start_anchor_index, from_index]] = chaining_lower_bounds
                .primary_lower_bound(gap1, gap2)
                .max(max_exact_cost_function_cost + Cost::from_usize(1))
                .min(primary[[primary_start_anchor_index, from_index]]);

            // Fill primary from from_index with exact values.
            additional_primary_targets_output.clear();
            primary_aligner.align_until_cost_limit(
                from_anchor.end(k),
                max_exact_cost_function_cost,
                &mut additional_primary_targets_output,
                &mut PanicOnExtend,
            );
            additional_primary_targets_output.sort_unstable();
            for (to_index, cost) in anchors
                .primary_anchor_to_index_iter(additional_primary_targets_output.iter().copied())
            {
                let to_index = to_index + 1;
                primary[[from_index, to_index]] = primary[[from_index, to_index]].min(cost);
                if cost <= max_exact_cost_function_cost {
                    primary.set_exact(from_index, to_index);
                }
            }
            if let Some(to_end_cost) = additional_primary_targets_output
                .iter()
                .rev()
                .map_while(|(anchor, cost)| (anchor.start() == end).then_some(cost))
                .copied()
                .last()
            {
                primary[[from_index, primary_end_anchor_index]] = to_end_cost;
                if to_end_cost <= max_exact_cost_function_cost {
                    primary.set_exact(from_index, primary_end_anchor_index);
                }
            }

            // Fill remaining primary with lower bound.
            let (gap1, gap2) = from_anchor.chaining_gaps_to_end(end, k);
            if gap1 == 0 && gap2 == 0 {
                // Allow last primary anchor to chain without gaps to end anchor.
                primary[[from_index, primary_end_anchor_index]] = Cost::zero();
            } else {
                primary[[from_index, primary_end_anchor_index]] = chaining_lower_bounds
                    .primary_lower_bound(gap1, gap2)
                    .max(max_exact_cost_function_cost + Cost::from_usize(1))
                    .min(primary[[from_index, primary_end_anchor_index]]);
            }

            for (to_index, to_anchor) in anchors.enumerate_primaries() {
                let to_index = to_index + 1;
                if let Some((gap1, gap2)) = from_anchor.chaining_gaps(&to_anchor, k) {
                    primary[[from_index, to_index]] = chaining_lower_bounds
                        .primary_lower_bound(gap1, gap2)
                        .max(max_exact_cost_function_cost + Cost::from_usize(1))
                        .min(primary[[from_index, to_index]]);
                }
                if from_anchor.is_direct_predecessor_of(&to_anchor) {
                    debug_assert!(
                        primary[[from_index, to_index]].is_zero()
                            || primary[[from_index, to_index]] == Cost::max_value()
                    );
                    primary[[from_index, to_index]] = Cost::zero();
                }
            }
        }

        // Initialise secondaries with infinity.
        trace!("Initialise secondaries with infinity");
        let mut secondaries = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [
                        anchors.secondary_len(ts_kind),
                        anchors.secondary_len(ts_kind),
                    ],
                    Cost::max_value(),
                    true,
                )
            })
            .collect_array()
            .unwrap();
        trace!("Fill secondaries");
        for (ts_kind, secondary) in TsKind::iter().zip(&mut secondaries) {
            trace!("Fill secondaries S{}", ts_kind.digits());
            for (from_index, from_anchor) in anchors.enumerate_secondaries(ts_kind) {
                // Fill secondary from from_index with exact values.
                additional_secondary_targets_output.clear();
                secondary_aligner.align_until_cost_limit(
                    anchors.secondary(from_index, ts_kind).end(ts_kind, k),
                    max_exact_cost_function_cost,
                    &mut PanicOnExtend,
                    &mut additional_secondary_targets_output,
                );
                additional_secondary_targets_output.sort_unstable();
                for (to_index, cost) in anchors.secondary_anchor_to_index_iter(
                    additional_secondary_targets_output.iter().copied(),
                    ts_kind,
                ) {
                    secondary[[from_index, to_index]] = secondary[[from_index, to_index]].min(cost);
                    if cost <= max_exact_cost_function_cost {
                        secondary.set_exact(from_index, to_index);
                    }
                }

                // Fill remaining secondary with lower bound.
                for (to_index, to_anchor) in anchors.enumerate_secondaries(ts_kind) {
                    if let Some((gap1, gap2)) = from_anchor.chaining_gaps(&to_anchor, ts_kind, k) {
                        secondary[[from_index, to_index]] = chaining_lower_bounds
                            .secondary_lower_bound(gap1, gap2)
                            .max(max_exact_cost_function_cost + Cost::from_usize(1))
                            .min(secondary[[from_index, to_index]]);
                    }
                    if from_anchor.is_direct_predecessor_of(&to_anchor) {
                        debug_assert!(
                            secondary[[from_index, to_index]].is_zero()
                                || secondary[[from_index, to_index]] == Cost::max_value(),
                            "Direct predecessor relationship from S{}{} to S{}{} has cost {}",
                            ts_kind.digits(),
                            anchors.secondary(from_index, ts_kind),
                            ts_kind.digits(),
                            anchors.secondary(to_index, ts_kind),
                            secondary[[from_index, to_index]],
                        );
                        secondary[[from_index, to_index]] = Cost::zero();
                    }
                }
            }
        }

        // Initialise 12-jumps with infinity.
        trace!("Initialise 12-jumps with infinity");
        let mut jump_12s = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [primary_anchor_amount, anchors.secondary_len(ts_kind)],
                    Cost::max_value(),
                    false,
                )
            })
            .collect_array()
            .unwrap();

        let mut total_12_jump_exact_align_time = Duration::default();
        let mut total_12_jump_exact_evaluation_time = Duration::default();
        let mut total_12_jump_exact_evaluations = 0;
        let mut total_12_jump_exact_evaluation_opened_nodes = 0usize;
        for (ts_kind, jump_12) in TsKind::iter().zip(&mut jump_12s) {
            trace!("Filling 12-jumps to S{}", ts_kind.digits());
            for (from_index, from_anchor) in anchors.enumerate_primaries() {
                let from_index = from_index + 1;

                // Fill 12-jumps with lower bound.
                let mut eligible_anchors = Vec::new();
                for (to_index, to_anchor) in anchors.enumerate_secondaries(ts_kind) {
                    if let Some(gap) = from_anchor.chaining_jump_gap(&to_anchor, ts_kind, k) {
                        let lower_bound = chaining_lower_bounds.jump_12_lower_bound(gap);
                        jump_12[[from_index, to_index]] = lower_bound.max(
                            max_exact_cost_function_cost
                                + chaining_lower_bounds
                                    .alignment_costs()
                                    .ts_base_cost
                                    .get(ts_kind)
                                + Cost::from_usize(1),
                        );
                        if lower_bound
                            <= max_exact_cost_function_cost
                                + chaining_lower_bounds
                                    .alignment_costs()
                                    .ts_base_cost
                                    .get(ts_kind)
                        {
                            eligible_anchors.push(to_anchor);
                        }
                    }
                }

                if !eligible_anchors.is_empty() {
                    total_12_jump_exact_evaluations += 1;
                    eligible_anchors.sort_unstable_by_key(|anchor| {
                        anchor.start(ts_kind).secondary_ordinate_ancestor().unwrap()
                    });
                    // TODO: use eligible anchors for much more detailed filtering or for an A* lower bound in alignment.
                    let eligible_anchors = eligible_anchors;
                    let min_eligible_ancestor = eligible_anchors
                        .first()
                        .unwrap()
                        .start(ts_kind)
                        .secondary_ordinate_ancestor()
                        .unwrap();
                    let align_end = AlignmentCoordinates::new_secondary(
                        min_eligible_ancestor,
                        sequences
                            .secondary_end(ts_kind)
                            .secondary_ordinate_descendant()
                            .unwrap(),
                        ts_kind,
                    );

                    // Correct 12-jumps from from_index to exact values.
                    additional_secondary_targets_output.clear();
                    // FIXME: this finds all positions in the ancestor, but almost none of them belong to any anchor.
                    let start = Instant::now();
                    total_12_jump_exact_evaluation_opened_nodes += ts_12_jump_aligner
                        .align_until_cost_limit(
                            anchors.primary(from_index - 1).end(k),
                            align_end,
                            ts_kind,
                            max_exact_cost_function_cost
                                + chaining_lower_bounds
                                    .alignment_costs()
                                    .ts_base_cost
                                    .get(ts_kind),
                            &mut additional_secondary_targets_output,
                        );
                    let after_align = Instant::now();
                    additional_secondary_targets_output.sort_unstable();
                    let mut additional_secondary_anchor_count = 0;
                    for (to_index, cost) in anchors.secondary_anchor_to_index_iter(
                        additional_secondary_targets_output.iter().copied(),
                        ts_kind,
                    ) {
                        additional_secondary_anchor_count += 1;
                        jump_12[[from_index, to_index]] = jump_12[[from_index, to_index]].min(cost);
                        if cost
                            <= max_exact_cost_function_cost
                                + chaining_lower_bounds
                                    .alignment_costs()
                                    .ts_base_cost
                                    .get(ts_kind)
                        {
                            jump_12.set_exact(from_index, to_index);
                        }
                    }
                    let end = Instant::now();
                    total_12_jump_exact_align_time += after_align - start;
                    total_12_jump_exact_evaluation_time += end - after_align;
                    trace!(
                        "Found {additional_secondary_anchor_count}/{} additional anchors from P[{}] to S{}",
                        additional_secondary_targets_output.len(),
                        from_index - 1,
                        ts_kind.digits(),
                    );
                }
            }
        }

        debug!(
            "Exact cost evaluation for 12-jumps took {:.0}s to align and {:.0}s to evaluate",
            total_12_jump_exact_align_time.as_secs_f64(),
            total_12_jump_exact_evaluation_time.as_secs_f64(),
        );
        debug!(
            "Exact cost evaluation for 12-jumps opened on average {} nodes (without skipped evaluations)",
            total_12_jump_exact_evaluation_opened_nodes
                .checked_div(total_12_jump_exact_evaluations)
                .unwrap_or(0),
        );

        // Initialise 34-jumps with infinity.
        trace!("Initialise 34-jumps with infinity");
        let mut jump_34s = TsKind::iter()
            .map(|ts_kind| {
                ChainingCostArray::new_from_cost(
                    [anchors.secondary_len(ts_kind), primary_anchor_amount],
                    Cost::max_value(),
                    false,
                )
            })
            .collect_array()
            .unwrap();
        for (ts_kind, jump_34) in TsKind::iter().zip(&mut jump_34s) {
            trace!("Filling 34-jumps from S{}", ts_kind.digits());
            for (from_index, from_anchor) in anchors.enumerate_secondaries(ts_kind) {
                // Fill 34-jumps from from_index with exact values.
                additional_primary_targets_output.clear();
                ts_34_jump_aligner.align_until_cost_limit(
                    anchors.secondary(from_index, ts_kind).end(ts_kind, k),
                    max_exact_cost_function_cost,
                    &mut additional_primary_targets_output,
                );
                additional_primary_targets_output.sort_unstable();
                for (to_index, cost) in anchors
                    .primary_anchor_to_index_iter(additional_primary_targets_output.iter().copied())
                {
                    let to_index = to_index + 1;
                    jump_34[[from_index, to_index]] = jump_34[[from_index, to_index]].min(cost);
                    if cost <= max_exact_cost_function_cost {
                        jump_34.set_exact(from_index, to_index);
                    }
                }
                if let Some(to_end_cost) = additional_primary_targets_output
                    .iter()
                    .rev()
                    .map_while(|(anchor, cost)| (anchor.start() == end).then_some(cost))
                    .copied()
                    .last()
                {
                    jump_34[[from_index, primary_end_anchor_index]] = to_end_cost;
                    if to_end_cost <= max_exact_cost_function_cost {
                        jump_34.set_exact(from_index, primary_end_anchor_index);
                    }
                }

                // Fill remaining 34-jumps with lower bound.
                for (to_index, to_anchor) in anchors.enumerate_primaries() {
                    let to_index = to_index + 1;
                    if let Some(gap) = from_anchor.chaining_jump_gap(&to_anchor, ts_kind, k) {
                        jump_34[[from_index, to_index]] = chaining_lower_bounds
                            .jump_34_lower_bound(gap)
                            .max(max_exact_cost_function_cost + Cost::from_usize(1))
                            .min(jump_34[[from_index, to_index]]);
                    }
                }
            }
        }

        trace!("Fill jumps from start and to end");
        for (ts_kind, (jump_12, jump_34)) in
            TsKind::iter().zip(jump_12s.iter_mut().zip(&mut jump_34s))
        {
            // Fill 12-jumps from start with exact values.
            additional_secondary_targets_output.clear();
            ts_12_jump_aligner.align_until_cost_limit(
                start,
                sequences.secondary_end(ts_kind),
                ts_kind,
                max_exact_cost_function_cost
                    + chaining_lower_bounds
                        .alignment_costs()
                        .ts_base_cost
                        .get(ts_kind),
                &mut additional_secondary_targets_output,
            );
            additional_secondary_targets_output.sort_unstable();
            for (to_index, cost) in anchors.secondary_anchor_to_index_iter(
                additional_secondary_targets_output.iter().copied(),
                ts_kind,
            ) {
                jump_12[[primary_start_anchor_index, to_index]] = cost;
                if cost
                    <= max_exact_cost_function_cost
                        + chaining_lower_bounds
                            .alignment_costs()
                            .ts_base_cost
                            .get(ts_kind)
                {
                    jump_12.set_exact(primary_start_anchor_index, to_index);
                }
            }

            for (index, anchor) in anchors.enumerate_secondaries(ts_kind) {
                // Fill remaining 12-jumps from start with lower bound.
                let gap = anchor.chaining_jump_gap_from_start(start, ts_kind);
                jump_12[[primary_start_anchor_index, index]] = chaining_lower_bounds
                    .jump_12_lower_bound(gap)
                    .max(
                        max_exact_cost_function_cost
                            + chaining_lower_bounds
                                .alignment_costs()
                                .ts_base_cost
                                .get(ts_kind)
                            + Cost::from_usize(1),
                    )
                    .min(jump_12[[primary_start_anchor_index, index]]);

                // Fill remaining 34-jumps to end with lower bound.
                let gap = anchor.chaining_jump_gap_to_end(end, ts_kind, k);
                jump_34[[index, primary_end_anchor_index]] = chaining_lower_bounds
                    .jump_34_lower_bound(gap)
                    .max(max_exact_cost_function_cost + Cost::from_usize(1))
                    .min(jump_34[[index, primary_end_anchor_index]]);
            }
        }

        let end_time = Instant::now();
        let duration = end_time - start_time;
        debug!(
            "Initialising chaining cost function took {:.0}ms",
            duration.as_secs_f64() * 1e3,
        );

        Self {
            primary,
            secondaries,
            jump_12s,
            jump_34s,
        }
    }

    fn primary_start_anchor_index() -> AnchorIndex {
        AnchorIndex::zero()
    }

    fn primary_end_anchor_index(&self) -> AnchorIndex {
        self.primary.dim().0 - 1
    }

    fn secondary_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.secondaries[ts_kind.index()]
    }

    fn secondary_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.secondaries[ts_kind.index()]
    }

    fn jump_12_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.jump_12s[ts_kind.index()]
    }

    fn jump_12_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.jump_12s[ts_kind.index()]
    }

    fn jump_34_array(&self, ts_kind: TsKind) -> &ChainingCostArray<Cost> {
        &self.jump_34s[ts_kind.index()]
    }

    fn jump_34_array_mut(&mut self, ts_kind: TsKind) -> &mut ChainingCostArray<Cost> {
        &mut self.jump_34s[ts_kind.index()]
    }

    pub fn primary(&self, from_primary_index: AnchorIndex, to_primary_index: AnchorIndex) -> Cost {
        self.primary[[from_primary_index + 1, to_primary_index + 1]]
    }

    pub fn primary_from_start(&self, primary_index: AnchorIndex) -> Cost {
        self.primary[[Self::primary_start_anchor_index(), primary_index + 1]]
    }

    pub fn primary_to_end(&self, primary_index: AnchorIndex) -> Cost {
        self.primary[[primary_index + 1, self.primary_end_anchor_index()]]
    }

    pub fn start_to_end(&self) -> Cost {
        self.primary[[
            Self::primary_start_anchor_index(),
            self.primary_end_anchor_index(),
        ]]
    }

    pub fn jump_12(
        &self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.jump_12_array(ts_kind)[[from_primary_index + 1, to_secondary_index]]
    }

    pub fn jump_12_from_start(&self, to_secondary_index: AnchorIndex, ts_kind: TsKind) -> Cost {
        self.jump_12_array(ts_kind)[[Self::primary_start_anchor_index(), to_secondary_index]]
    }

    pub fn secondary(
        &self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.secondary_array(ts_kind)[[from_secondary_index, to_secondary_index]]
    }

    pub fn jump_34(
        &self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> Cost {
        self.jump_34_array(ts_kind)[[from_secondary_index, to_primary_index + 1]]
    }

    pub fn jump_34_to_end(&self, from_secondary_index: AnchorIndex, ts_kind: TsKind) -> Cost {
        self.jump_34_array(ts_kind)[[from_secondary_index, self.primary_end_anchor_index()]]
    }

    pub fn is_primary_exact(
        &self,
        from_primary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
    ) -> bool {
        self.primary
            .is_exact(from_primary_index + 1, to_primary_index + 1)
    }

    pub fn is_primary_from_start_exact(&self, primary_index: AnchorIndex) -> bool {
        self.primary
            .is_exact(Self::primary_start_anchor_index(), primary_index + 1)
    }

    pub fn is_primary_to_end_exact(&self, primary_index: AnchorIndex) -> bool {
        self.primary
            .is_exact(primary_index + 1, self.primary_end_anchor_index())
    }

    pub fn is_start_to_end_exact(&self) -> bool {
        self.primary.is_exact(
            Self::primary_start_anchor_index(),
            self.primary_end_anchor_index(),
        )
    }

    pub fn is_jump_12_exact(
        &self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_12_array(ts_kind)
            .is_exact(from_primary_index + 1, to_secondary_index)
    }

    pub fn is_jump_12_from_start_exact(
        &self,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_12_array(ts_kind)
            .is_exact(Self::primary_start_anchor_index(), to_secondary_index)
    }

    pub fn is_secondary_exact(
        &self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.secondary_array(ts_kind)
            .is_exact(from_secondary_index, to_secondary_index)
    }

    pub fn is_jump_34_exact(
        &self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_34_array(ts_kind)
            .is_exact(from_secondary_index, to_primary_index + 1)
    }

    pub fn is_jump_34_to_end_exact(
        &self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
    ) -> bool {
        self.jump_34_array(ts_kind)
            .is_exact(from_secondary_index, self.primary_end_anchor_index())
    }

    pub fn update_primary(
        &mut self,
        from_primary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary
                .set_exact(from_primary_index + 1, to_primary_index + 1);
        }
        let target = &mut self.primary[[from_primary_index + 1, to_primary_index + 1]];
        assert!(
            *target <= cost,
            "Target is larger than cost.\ntarget: {target}; cost: {cost}; from_primary_index: {from_primary_index}; to_primary_index: {to_primary_index}; is_exact: {is_exact}",
        );
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_from_start(
        &mut self,
        primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.primary
                .set_exact(Self::primary_start_anchor_index(), primary_index + 1);
        }
        let target = &mut self.primary[[Self::primary_start_anchor_index(), primary_index + 1]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_primary_to_end(
        &mut self,
        primary_index: AnchorIndex,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let end_index = self.primary_end_anchor_index();
        if is_exact {
            self.primary.set_exact(primary_index + 1, end_index);
        }
        let target = &mut self.primary[[primary_index + 1, end_index]];
        assert!(
            *target <= cost,
            "Target is larger than cost.\ntarget: {target}; cost: {cost}; from_primary_index: {primary_index}; is_exact: {is_exact}",
        );
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_start_to_end(&mut self, cost: Cost, is_exact: bool) -> bool {
        debug_assert!(
            cost < Cost::max_value(),
            "Chaining start to end is updated to infinite costs, but this chaining should always be possible."
        );
        let end_index = self.primary.dim().1 - 1;
        if is_exact {
            self.primary
                .set_exact(Self::primary_start_anchor_index(), end_index);
        }
        let target = &mut self.primary[[Self::primary_start_anchor_index(), end_index]];
        assert!(
            *target <= cost,
            "Start-to-end target is larger than cost.\ntarget: {target}; cost: {cost}; is_exact: {is_exact}",
        );
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12(
        &mut self,
        from_primary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_12_array_mut(ts_kind)
                .set_exact(from_primary_index + 1, to_secondary_index);
        }
        let target =
            &mut self.jump_12_array_mut(ts_kind)[[from_primary_index + 1, to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_12_from_start(
        &mut self,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_12_array_mut(ts_kind)
                .set_exact(Self::primary_start_anchor_index(), to_secondary_index);
        }
        let target = &mut self.jump_12_array_mut(ts_kind)
            [[Self::primary_start_anchor_index(), to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_secondary(
        &mut self,
        from_secondary_index: AnchorIndex,
        to_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.secondary_array_mut(ts_kind)
                .set_exact(from_secondary_index, to_secondary_index);
        }
        let target =
            &mut self.secondary_array_mut(ts_kind)[[from_secondary_index, to_secondary_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34(
        &mut self,
        from_secondary_index: AnchorIndex,
        to_primary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        if is_exact {
            self.jump_34_array_mut(ts_kind)
                .set_exact(from_secondary_index, to_primary_index + 1);
        }
        let target =
            &mut self.jump_34_array_mut(ts_kind)[[from_secondary_index, to_primary_index + 1]];
        assert!(
            *target <= cost,
            "Reducing 34-jump from S{}-{from_secondary_index} to P-{to_primary_index}: cost from {target} to {cost}",
            ts_kind.digits()
        );
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn update_jump_34_to_end(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        cost: Cost,
        is_exact: bool,
    ) -> bool {
        let end_index = self.primary_end_anchor_index();
        if is_exact {
            self.jump_34_array_mut(ts_kind)
                .set_exact(from_secondary_index, end_index);
        }
        let target = &mut self.jump_34_array_mut(ts_kind)[[from_secondary_index, end_index]];
        assert!(*target <= cost);
        let result = *target < cost;
        *target = cost;
        result
    }

    pub fn iter_primary_in_cost_order(
        &mut self,
        from_primary_index: AnchorIndex,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary
            .iter_in_cost_order_from(from_primary_index + 1, offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn iter_primary_from_start_in_cost_order(
        &mut self,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.primary
            .iter_in_cost_order_from(AnchorIndex::zero(), offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn iter_jump_12_in_cost_order(
        &mut self,
        from_primary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_12_array_mut(ts_kind)
            .iter_in_cost_order_from(from_primary_index + 1, offset)
    }

    pub fn iter_jump_12_from_start_in_cost_order(
        &mut self,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_12_array_mut(ts_kind)
            .iter_in_cost_order_from(Self::primary_start_anchor_index(), offset)
    }

    pub fn iter_secondary_in_cost_order(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.secondary_array_mut(ts_kind)
            .iter_in_cost_order_from(from_secondary_index, offset)
    }

    pub fn iter_jump_34_in_cost_order(
        &mut self,
        from_secondary_index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    ) -> impl Iterator<Item = (AnchorIndex, Cost)>
    where
        Cost: Copy + Ord,
    {
        self.jump_34_array_mut(ts_kind)
            .iter_in_cost_order_from(from_secondary_index, offset)
            .filter(|(to_primary_index, _)| *to_primary_index != AnchorIndex::zero())
            .map(|(to_primary_index, chaining_cost)| (to_primary_index - 1, chaining_cost))
    }

    pub fn update_additional_primary_targets(
        &mut self,
        from_primary_index: AnchorIndex,
        additional_targets: &mut [(PrimaryAnchor, Cost)],
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_primary_index, cost) in
            anchors.primary_anchor_to_index_iter(additional_targets.iter().copied())
        {
            if self.is_primary_exact(from_primary_index, to_primary_index) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(self.primary(from_primary_index, to_primary_index), cost);
            } else {
                self.update_primary(from_primary_index, to_primary_index, cost, true);
            }
        }
    }

    pub fn update_additional_primary_targets_from_start(
        &mut self,
        additional_targets: &mut [(PrimaryAnchor, Cost)],
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_primary_index, cost) in
            anchors.primary_anchor_to_index_iter(additional_targets.iter().copied())
        {
            if self.is_primary_from_start_exact(to_primary_index) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(self.primary_from_start(to_primary_index), cost);
            } else {
                self.update_primary_from_start(to_primary_index, cost, true);
            }
        }
    }

    pub fn update_additional_secondary_targets(
        &mut self,
        from_secondary_index: AnchorIndex,
        additional_targets: &mut [(SecondaryAnchor, Cost)],
        ts_kind: TsKind,
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_secondary_index, cost) in
            anchors.secondary_anchor_to_index_iter(additional_targets.iter().copied(), ts_kind)
        {
            if self.is_secondary_exact(from_secondary_index, to_secondary_index, ts_kind) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(
                    self.secondary(from_secondary_index, to_secondary_index, ts_kind),
                    cost,
                );
            } else {
                self.update_secondary(
                    from_secondary_index,
                    to_secondary_index,
                    ts_kind,
                    cost,
                    true,
                );
            }
        }
    }

    pub fn update_additional_12_jump_targets(
        &mut self,
        from_primary_index: AnchorIndex,
        additional_targets: &mut [(SecondaryAnchor, Cost)],
        ts_kind: TsKind,
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_secondary_index, cost) in
            anchors.secondary_anchor_to_index_iter(additional_targets.iter().copied(), ts_kind)
        {
            if self.is_jump_12_exact(from_primary_index, to_secondary_index, ts_kind) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(
                    self.jump_12(from_primary_index, to_secondary_index, ts_kind),
                    cost,
                    "Jump12: Previous exact cost was {} but additional target cost is {cost}.\nFrom P-{from_primary_index} to S{}-{to_secondary_index}.",
                    self.jump_12(from_primary_index, to_secondary_index, ts_kind),
                    ts_kind.digits(),
                );
            } else {
                self.update_jump_12(from_primary_index, to_secondary_index, ts_kind, cost, true);
            }
        }
    }

    pub fn update_additional_12_jump_targets_from_start(
        &mut self,
        additional_targets: &mut [(SecondaryAnchor, Cost)],
        ts_kind: TsKind,
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_secondary_index, cost) in
            anchors.secondary_anchor_to_index_iter(additional_targets.iter().copied(), ts_kind)
        {
            if self.is_jump_12_from_start_exact(to_secondary_index, ts_kind) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(self.jump_12_from_start(to_secondary_index, ts_kind), cost,);
            } else {
                self.update_jump_12_from_start(to_secondary_index, ts_kind, cost, true);
            }
        }
    }

    pub fn update_additional_34_jump_targets(
        &mut self,
        from_secondary_index: AnchorIndex,
        additional_targets: &mut [(PrimaryAnchor, Cost)],
        ts_kind: TsKind,
        anchors: &Anchors,
        total_redundant_gap_fillings: &mut u64,
    ) {
        additional_targets.sort_unstable();
        for (to_primary_index, cost) in
            anchors.primary_anchor_to_index_iter(additional_targets.iter().copied())
        {
            if self.is_jump_34_exact(from_secondary_index, to_primary_index, ts_kind) {
                *total_redundant_gap_fillings += 1;
                debug_assert_eq!(
                    self.jump_34(from_secondary_index, to_primary_index, ts_kind),
                    cost,
                    "Jump34: Previous exact cost was {} but additional target cost is {cost}.\nFrom S{}-{from_secondary_index}{} to P-{to_primary_index}{}.",
                    self.jump_34(from_secondary_index, to_primary_index, ts_kind),
                    ts_kind.digits(),
                    anchors.secondary(from_secondary_index, ts_kind),
                    anchors.primary(to_primary_index),
                );
            } else {
                self.update_jump_34(from_secondary_index, to_primary_index, ts_kind, cost, true);
            }
        }
    }
}
