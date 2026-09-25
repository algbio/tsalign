use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{
        Alignment,
        coordinates::{
            AlignmentCoordinates, AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
    },
    alignment_history::{extension_graph::HistoryExtensionGraph, history_vec::AlignmentHistory},
    costs::GapAffineCosts,
    inexact_chaining::gap_affine::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost, AlignmentHistoryVec> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    history_graph: HistoryExtensionGraph<AlignmentHistoryVec>,
    sequences: &'sequences AlignmentSequences,
    cost_table: &'cost_table GapAffineCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    anchor_k: u8,
    max_anchor_mutations: u8,
}

impl<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory>
    GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost, AlignmentHistoryVec>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        cost_table: &'cost_table GapAffineCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self {
        Self {
            a_star_buffers: Some(Default::default()),
            history_graph: HistoryExtensionGraph::new(),
            sequences,
            cost_table,
            rc_fn,
            anchor_k,
            max_anchor_mutations,
        }
    }

    /// Align from start to end.
    ///
    /// Additionally continue the alignment to all nodes with the same cost as the alignment cost from start to end.
    ///
    /// Collect all closed nodes into the given output lists.
    /// Note that the output lists may contain duplicate anchors with different cost.
    pub fn align(
        &mut self,
        start: impl Into<AlignmentCoordinates>,
        end: impl Into<AlignmentCoordinates>,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) -> (Cost, Alignment) {
        let start = start.into();
        let end = end.into();

        assert!(
            start.is_primary() && end.is_primary() || start.is_secondary() && end.is_secondary()
        );

        let context = Context::new(
            self.cost_table,
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
        let (cost, alignment) = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => {
                let alignment = a_star.reconstruct_path().into();
                (cost, alignment)
            }
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => (Cost::max_value(), Vec::new().into()),
        };
        a_star.search_until_with_target_policy(|_, node| node.cost > cost, true);

        Self::fill_additional_targets(
            &a_star,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output lists.
    /// Note that the output lists may contain duplicate coordinates with different cost.
    pub fn align_until_cost_limit(
        &mut self,
        start: impl Into<AlignmentCoordinates>,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) {
        let start = start.into();
        let end = self.sequences.end(start.ts_kind());
        debug_assert!(
            start.is_primary() && end.is_primary() || start.is_secondary() && end.is_secondary()
        );

        let context = Context::new(
            self.cost_table,
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

        Self::fill_additional_targets(
            &a_star,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());
    }

    fn fill_additional_targets(
        a_star: &AStar<Context<Cost, AlignmentHistoryVec>>,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) {
        additional_primary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter_map(|node| Some((node.identifier.coordinates.into_primary()?, node.cost))),
        );
        additional_secondary_targets_output.extend(a_star.iter_closed_nodes().filter_map(|node| {
            Some((
                node.identifier.coordinates.into_secondary()?.into_any(),
                node.cost,
            ))
        }));
    }
}
