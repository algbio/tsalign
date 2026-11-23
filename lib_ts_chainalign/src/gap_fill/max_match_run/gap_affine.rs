use generic_a_star::cost::AStarCost;

use crate::{
    alignment::{
        Alignment,
        coordinates::{
            AlignmentCoordinates, AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
    },
    costs::GapAffineCosts,
    gap_fill::max_match_run::{
        alignment_search::{
            Datastructures, Dijkstra, iter_min_cost_nodes_matching, new_dijkstra,
            search_best_target, search_until_cost_limit,
        },
        gap_affine::algo::{Context, Identifier, Node},
    },
};

mod algo;
#[cfg(test)]
mod tests;

/// Aligns with a bicriteria Dijkstra, see [`alignment_search`](crate::gap_fill::max_match_run::alignment_search).
pub struct GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost> {
    datastructures: Option<Datastructures<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    cost_table: &'cost_table GapAffineCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
    enforcement_offset: u8,
    enforce_non_match_next_to_anchors: bool,
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
            enforce_non_match_next_to_anchors: enforcement_offset == 0,
        }
    }

    /// Do not enforce non-matches next to anchors.
    ///
    /// By default, non-matches are enforced next to anchors for exact anchors (`enforcement_offset == 0`),
    /// where the lower bounds assume it.
    /// This allows to use the aligner for alignments between arbitrary coordinates.
    #[cfg(test)]
    pub fn without_enforced_non_match_next_to_anchors(mut self) -> Self {
        self.enforce_non_match_next_to_anchors = false;
        self
    }

    /// Returns true if non-matches are enforced next to anchors.
    fn enforce_non_match_next_to_anchors(&self) -> bool {
        self.enforce_non_match_next_to_anchors
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
            self.enforce_non_match_next_to_anchors(),
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        let (cost, alignment) = search_best_target(&mut dijkstra);

        self.fill_additional_targets(
            &dijkstra,
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
            self.enforce_non_match_next_to_anchors(),
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        search_until_cost_limit(&mut dijkstra, cost_limit);

        self.fill_additional_targets(
            &dijkstra,
            start,
            cost_limit,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.datastructures = Some(dijkstra.into_datastructures());
    }

    /// Output the minimum cost of each node that is within the cost limit.
    fn fill_additional_targets(
        &self,
        dijkstra: &Dijkstra<Context<'_, '_, '_, Cost>>,
        start: AlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(AnySecondaryAlignmentCoordinates, Cost)>,
    ) {
        let prune_list = dijkstra.prune_list();
        let context = dijkstra.context();
        for node in iter_min_cost_nodes_matching(prune_list, cost_limit, |node| {
            context.fulfils_end_condition(node)
        }) {
            let identifier = &node.identifier;
            let coordinates = identifier.coordinates;
            let cost = node.cost.cost;

            if (start == coordinates && !self.allow_direct_chaining(start, coordinates))
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
                .min_first_node_at_matching(&pair_identifier, |node| {
                    context.fulfils_end_condition(node)
                })
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
            false,
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        let (cost, alignment) = search_best_target(&mut dijkstra);
        self.datastructures = Some(dijkstra.into_datastructures());

        (cost, alignment)
    }
}
