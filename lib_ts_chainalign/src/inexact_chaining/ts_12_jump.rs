use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{
        Alignment,
        coordinates::{
            AlignmentCoordinates, AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            SpecificSecondaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    alignment_history::{extension_graph::HistoryExtensionGraph, history_vec::AlignmentHistory},
    costs::AlignmentCosts,
    inexact_chaining::ts_12_jump::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts12JumpAligner<
    'sequences,
    'alignment_costs,
    'rc_fn,
    Cost: AStarCost,
    AlignmentHistoryVec,
> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    history_graph: HistoryExtensionGraph<AlignmentHistoryVec>,
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    anchor_k: u8,
    max_anchor_mutations: u8,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory>
    Ts12JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self {
        Self {
            a_star_buffers: Some(Default::default()),
            history_graph: HistoryExtensionGraph::new(),
            sequences,
            alignment_costs,
            rc_fn,
            anchor_k,
            max_anchor_mutations,
        }
    }

    /// Align from start to end.
    ///
    /// Additionally continue the alignment to all nodes with the same cost as the alignment cost from start to end.
    ///
    /// Collect all closed nodes into the given output list.
    /// Note that the output list may contain duplicate anchors with different cost.
    pub fn align(
        &mut self,
        start: PrimaryAlignmentCoordinates,
        end: SpecificSecondaryAlignmentCoordinates,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) -> (Cost, Alignment) {
        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            &mut self.history_graph,
            start,
            end,
            self.anchor_k,
            self.max_anchor_mutations,
        );
        let mut a_star = AStar::<_>::new_with_buffers(context, self.a_star_buffers.take().unwrap());

        a_star.initialise();
        let (cost, alignment) = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => (cost, a_star.reconstruct_path().into()),
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => (Cost::max_value(), Vec::new().into()),
        };

        a_star.search_until_with_target_policy(|_, node| node.cost > cost, true);
        Self::fill_additional_targets(&a_star, end.ts_kind(), additional_secondary_targets_output);
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output list.
    /// Note that the output list may contain duplicate anchors with different cost.
    pub fn align_until_cost_limit(
        &mut self,
        start: PrimaryAlignmentCoordinates,
        end: SpecificSecondaryAlignmentCoordinates,
        cost_limit: Cost,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) -> usize {
        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            &mut self.history_graph,
            start,
            end,
            self.anchor_k,
            self.max_anchor_mutations,
        );
        let mut a_star = AStar::new_with_buffers(context, self.a_star_buffers.take().unwrap());
        a_star.initialise();
        a_star.search_until_with_target_policy(|_, node| node.cost > cost_limit, true);

        Self::fill_additional_targets(&a_star, end.ts_kind(), additional_secondary_targets_output);

        let opened_node_amount = a_star.performance_counters().opened_nodes;
        self.a_star_buffers = Some(a_star.into_buffers());
        opened_node_amount
    }

    fn fill_additional_targets(
        a_star: &AStar<Context<Cost, AlignmentHistoryVec>>,
        ts_kind: TsKind,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) {
        additional_secondary_targets_output.extend(a_star.iter_closed_nodes().filter_map(|node| {
            if let AlignmentCoordinates::Secondary(secondary) = node.identifier.coordinates(ts_kind)
            {
                Some((secondary.into(), node.cost))
            } else {
                None
            }
        }));
    }
}
