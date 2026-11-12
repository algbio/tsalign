use generic_a_star::cost::AStarCost;
use multicriteria_dijkstra::{
    MulticriteriaDijkstra, MulticriteriaDijkstraDatastructures,
    prune_list::{MulticriteriaDijkstraPruneList, NodeWithOptionalPredecessor},
};

use crate::{
    alignment::{
        Alignment,
        coordinates::{
            AlignmentCoordinates, AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
    },
    costs::GapAffineCosts,
    max_match_run_chaining::{
        gap_affine::algo::{Context, Identifier, Node, NodeCost},
        multicriteria_dijkstra::PruneList,
    },
};

mod algo;
#[cfg(test)]
mod tests;

type Datastructures<Cost> = MulticriteriaDijkstraDatastructures<PruneList<Node<Cost>>>;
type Dijkstra<'context, Cost> =
    MulticriteriaDijkstra<PruneList<Node<Cost>>, Context<'context, 'context, 'context, Cost>>;

/// Aligns with a bicriteria Dijkstra, minimising the alignment cost and the length of the current match run.
///
/// Target nodes are not reported to the prune list, since target pruning would also prune nodes with
/// the same alignment cost as a target, but a larger match run.
/// Instead, nodes whose alignment cost exceeds the best target cost (or the cost limit) are not inserted.
pub struct GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost> {
    datastructures: Option<Datastructures<Cost>>,
    sequences: &'sequences AlignmentSequences,
    cost_table: &'cost_table GapAffineCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
    enforcement_offset: u8,
}

impl<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost>
    GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        cost_table: &'cost_table GapAffineCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
        enforcement_offset: u8,
    ) -> Self {
        Self {
            datastructures: Some(Default::default()),
            sequences,
            cost_table,
            rc_fn,
            max_match_run,
            enforcement_offset,
        }
    }

    /// Evaluate if `allow_direct_chaining` should be set in the alignment context.
    fn allow_direct_chaining(
        &self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> bool {
        start == self.sequences.primary_start().into()
            || end == self.sequences.primary_end().into()
            // If enforcement_offset > 0 then we are in inexact anchor chaining, where direct chaining is always allowed.
            || self.enforcement_offset > 0
    }

    /// Evaluate if `allow_all_matches` should be set in the alignment context.
    fn allow_all_matches(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        let minimum_primary_sequence_length = (self.sequences.primary_end().a()
            - self.sequences.primary_start().a())
        .min(self.sequences.primary_end().b() - self.sequences.primary_start().b());
        (start == self.sequences.primary_start().into()
            && end == self.sequences.primary_end().into()
            && u32::try_from(minimum_primary_sequence_length).unwrap() <= self.max_match_run)
            // If enforcement_offset > 0 then we are in inexact anchor chaining, where all matches is always allowed.
            || self.enforcement_offset > 0
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
            start,
            end,
            self.allow_direct_chaining(start, end),
            self.allow_all_matches(start, end),
            self.max_match_run,
            self.enforcement_offset,
        );
        let mut dijkstra = self.dijkstra(context);
        let (cost, alignment) = Self::search_best_target(&mut dijkstra);

        self.fill_additional_targets(
            dijkstra.prune_list(),
            start,
            cost,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.datastructures = Some(dijkstra.into_datastructures());

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
            start,
            end,
            true,
            true,
            self.max_match_run,
            self.enforcement_offset,
        );
        let mut dijkstra = self.dijkstra(context);
        // The match run is not limited, such that only the alignment cost is limited.
        dijkstra.search_until_cost_limit(&NodeCost {
            cost: cost_limit,
            match_run: u32::MAX,
        });

        self.fill_additional_targets(
            dijkstra.prune_list(),
            start,
            cost_limit,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.datastructures = Some(dijkstra.into_datastructures());
    }

    /// Creates a Dijkstra with the reused datastructures of this aligner.
    ///
    /// The datastructures must be returned with [`MulticriteriaDijkstra::into_datastructures`].
    fn dijkstra<'context>(
        &mut self,
        context: Context<'context, 'context, 'context, Cost>,
    ) -> Dijkstra<'context, Cost> {
        let mut dijkstra = MulticriteriaDijkstra::from_datastructures(
            self.datastructures.take().unwrap(),
            context,
        );
        dijkstra.reset(|_| {});
        dijkstra
    }

    /// Search all nodes with an alignment cost of at most the cost of the best target.
    ///
    /// Returns the cost of the best target and its alignment,
    /// or `Cost::max_value()` and an empty alignment if no target was found.
    fn search_best_target(dijkstra: &mut Dijkstra<'_, Cost>) -> (Cost, Alignment) {
        use multicriteria_dijkstra::MulticriteriaDijkstraContext;

        // Nodes with a cost exceeding the best target cost found so far are not inserted.
        // Some nodes exceeding the final best target cost may be inserted before the best target is found.
        let mut best_target_cost = Cost::max_value();
        dijkstra.search_with_prune_condition(|context, node| {
            if context.is_target(node) {
                best_target_cost = best_target_cost.min(node.cost.cost);
            }
            node.cost.cost > best_target_cost
        });

        let prune_list = dijkstra.prune_list();
        let Some(best_target) = dijkstra
            .context()
            .target_identifiers()
            .filter_map(|identifier| prune_list.min_first_node_at(&identifier))
            .min_by_key(|node| node.cost.cost)
        else {
            return (Cost::max_value(), Vec::new().into());
        };

        // Backtrack the alignment.
        let mut alignment = Vec::new();
        let mut current = prune_list
            .iter_nodes_with_predecessor_at_identifier(&best_target.identifier)
            .find(|node| node.node == best_target)
            .unwrap();
        loop {
            let NodeWithOptionalPredecessor {
                node,
                predecessor_index,
            } = current;
            alignment.extend(node.predecessor_alignment_type);
            let Some(predecessor_index) = predecessor_index else {
                break;
            };
            current = prune_list.get(predecessor_index).unwrap();
        }
        alignment.reverse();

        (best_target.cost.cost, alignment.into())
    }

    /// Output the minimum cost of each node that is within the cost limit.
    fn fill_additional_targets(
        &self,
        prune_list: &PruneList<Node<Cost>>,
        start: AlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) {
        for (identifier, mut front) in prune_list.iter_fronts() {
            // The front is ordered by cost, so the first node has minimum cost.
            let Some(node) = front.next() else {
                continue;
            };
            let coordinates = identifier.coordinates;
            let cost = node.cost.cost;

            if cost > cost_limit
                || (start == coordinates && !self.allow_direct_chaining(start, coordinates))
                || (start != coordinates
                    && !identifier.has_non_match
                    && !self.allow_all_matches(start, coordinates))
            {
                continue;
            }

            // Filter duplicates.
            let pair_identifier = Identifier {
                has_non_match: !identifier.has_non_match,
                ..*identifier
            };
            if prune_list
                .min_first_node_at(&pair_identifier)
                .is_some_and(|pair| {
                    pair.cost.cost < cost || (pair.cost.cost == cost && !identifier.has_non_match)
                })
            {
                continue;
            }

            if let Some(primary) = coordinates.into_primary() {
                additional_primary_targets_output.extend(std::iter::once((primary, cost)));
            } else {
                additional_secondary_targets_output.extend(std::iter::once((
                    coordinates.into_secondary().unwrap().into(),
                    cost,
                )));
            }
        }
    }

    /// Align an anchor from start to end.
    ///
    /// None of the restrictions typically imposed by this aligner imply, e.g. `max_match_run` is ignored.
    pub fn align_anchor(
        &mut self,
        start: impl Into<AlignmentCoordinates>,
        end: impl Into<AlignmentCoordinates>,
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
            start,
            end,
            true,
            true,
            u32::MAX,
            0,
        );
        let mut dijkstra = self.dijkstra(context);
        let (cost, alignment) = Self::search_best_target(&mut dijkstra);
        self.datastructures = Some(dijkstra.into_datastructures());

        (cost, alignment)
    }
}
