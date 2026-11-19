use generic_a_star::cost::AStarCost;

use crate::{
    alignment::{
        Alignment,
        coordinates::{
            AlignmentCoordinates, PrimaryAlignmentCoordinates,
            SpecificSecondaryAlignmentCoordinates,
        },
        sequences::AlignmentSequences,
        ts_kind::{TsDescendant, TsKind},
    },
    costs::AlignmentCosts,
    max_match_run_chaining::{
        alignment_search::{
            Datastructures, iter_min_cost_nodes_matching, new_dijkstra, search_best_target,
            search_until_cost_limit,
        },
        multicriteria_dijkstra::PruneList,
        ts_34_jump::algo::{Context, Node},
    },
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    datastructures: Option<Datastructures<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
    enforcement_offset: u8,
    enforce_non_match_next_to_anchors: bool,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost>
    Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
        enforcement_offset: u8,
    ) -> Self {
        Self {
            datastructures: Some(Default::default()),
            sequences,
            alignment_costs,
            rc_fn,
            max_match_run,
            enforcement_offset,
            enforce_non_match_next_to_anchors: enforcement_offset == 0,
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
        start: SpecificSecondaryAlignmentCoordinates,
        end: PrimaryAlignmentCoordinates,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
    ) -> (Cost, Alignment) {
        let descendant_start = start.descendant();
        let enforce_non_match = self.enforce_non_match(
            descendant_start,
            end,
            Self::descendant_end(end, start.ts_kind()),
        );

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
            self.max_match_run,
            self.enforcement_offset,
            self.enforce_non_match_next_to_anchors(),
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        let (cost, alignment) = search_best_target(&mut dijkstra);

        Self::fill_additional_targets(
            dijkstra.prune_list(),
            self.enforce_non_match_next_to_anchors(),
            self.enforcement_offset,
            self.sequences.primary_end(),
            cost,
            descendant_start,
            start.ts_kind(),
            additional_primary_targets_output,
        );
        self.datastructures = Some(dijkstra.into_datastructures());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output list.
    /// Note that the output list may contain duplicate anchors with different cost.
    pub fn align_until_cost_limit(
        &mut self,
        start: SpecificSecondaryAlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
    ) {
        let end = self.sequences.primary_end();
        let descendant_start = start.descendant();
        let ts_kind = start.ts_kind();
        // Targets are not used when searching until the cost limit, so the end only bounds the search space.
        // The non-match condition for the reported coordinates is applied in `fill_additional_targets`,
        // this only keeps the context consistent with `align`.
        let enforce_non_match =
            self.enforce_non_match(descendant_start, end, Self::descendant_end(end, ts_kind));

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
            self.max_match_run,
            self.enforcement_offset,
            self.enforce_non_match_next_to_anchors(),
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        search_until_cost_limit(&mut dijkstra, cost_limit);

        Self::fill_additional_targets(
            dijkstra.prune_list(),
            self.enforce_non_match_next_to_anchors(),
            self.enforcement_offset,
            self.sequences.primary_end(),
            cost_limit,
            descendant_start,
            ts_kind,
            additional_primary_targets_output,
        );
        self.datastructures = Some(dijkstra.into_datastructures());
    }

    /// Do not enforce non-matches next to anchors.
    ///
    /// By default, non-matches are enforced next to anchors for exact anchors (`enforcement_offset == 0`),
    /// where the lower bounds assume it.
    /// This allows to use the aligner for alignments between arbitrary coordinates.
    pub fn without_enforced_non_match_next_to_anchors(mut self) -> Self {
        self.enforce_non_match_next_to_anchors = false;
        self
    }

    /// Returns true if non-matches are enforced next to anchors.
    fn enforce_non_match_next_to_anchors(&self) -> bool {
        self.enforce_non_match_next_to_anchors
    }

    /// Returns the descendant ordinate of the given end.
    fn descendant_end(end: PrimaryAlignmentCoordinates, ts_kind: TsKind) -> usize {
        match ts_kind.descendant {
            TsDescendant::Seq1 => end.a(),
            TsDescendant::Seq2 => end.b(),
        }
    }

    /// Returns true if a 34-jump from a secondary anchor ending at `descendant_start` in the descendant to `end`
    /// must contain a non-match.
    ///
    /// See [`Self::requires_non_match`].
    fn enforce_non_match(
        &self,
        descendant_start: usize,
        end: PrimaryAlignmentCoordinates,
        descendant_end: usize,
    ) -> bool {
        Self::requires_non_match(
            self.enforcement_offset,
            end == self.sequences.primary_end(),
            descendant_start,
            descendant_end,
        )
    }

    /// Returns true if a 34-jump whose descendant gap is `descendant_start..descendant_end` must contain a non-match.
    ///
    /// A non-match is enforced if there is a gap between the anchors in the descendant
    /// and the anchors are exact (`enforcement_offset == 0`).
    /// For exact anchors, an all-match jump part next to the secondary or primary anchor can be represented by a
    /// shifted anchor that is chained for free (see `is_direct_free_predecessor_of`).
    /// This matches the lower-bound computation of 34-jumps between exact anchors,
    /// and it discourages chains to deviate from the alignment geometry boundaries.
    ///
    /// Inexact anchors (`enforcement_offset > 0`) do not have such free shifted twins,
    /// so the jump may consist of matches only, and there is no other way to represent such a jump.
    /// The lower bounds of 34-jumps between inexact anchors allow this as well
    /// (they are built from the inexact gap-affine lower bounds, which allow all matches),
    /// and the gap-affine aligner does the same for `enforcement_offset > 0`.
    ///
    /// The end of the alignment is no anchor, so a jump to there may consist of matches only as well.
    /// The lower bounds of 34-jumps to the end allow this.
    ///
    /// If the descendant gap is empty, no non-match is enforced.
    /// Then the secondary and primary parts can only contain gaps in the descendant (deletions),
    /// and every such alignment is dominated by the bare jump directly to the target ancestor position,
    /// which costs nothing.
    ///
    /// This is independent of [`Self::enforce_non_match_next_to_anchors`].
    /// It is used for both the target in [`Self::align`] and the additional targets
    /// in [`Self::fill_additional_targets`], so that both agree on the exact costs.
    fn requires_non_match(
        enforcement_offset: u8,
        is_alignment_end: bool,
        descendant_start: usize,
        descendant_end: usize,
    ) -> bool {
        enforcement_offset == 0 && !is_alignment_end && descendant_start != descendant_end
    }

    /// Fills the primary coordinates reached with at most the given cost limit into the output.
    ///
    /// A node is reported if it would be a target of [`Self::align`] with its coordinates as `end`,
    /// i.e. if it has a non-match or if no non-match is required (see [`Self::requires_non_match`]).
    /// In particular, if the descendant gap is empty, all nodes are reported,
    /// including the ones reached only via deletions (which are dominated by the bare jump to the same coordinates).
    #[expect(clippy::too_many_arguments)]
    fn fill_additional_targets(
        prune_list: &PruneList<Node<Cost>>,
        enforce_non_match_next_to_anchors: bool,
        enforcement_offset: u8,
        primary_end: PrimaryAlignmentCoordinates,
        cost_limit: Cost,
        descendant_start: usize,
        ts_kind: TsKind,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAlignmentCoordinates, Cost)>,
    ) {
        additional_primary_targets_output.extend(
            // Additional targets are primary anchors or the end of the alignment,
            // so the same condition as in `Context::fulfils_end_condition` applies.
            iter_min_cost_nodes_matching(prune_list, cost_limit, |node| {
                !enforce_non_match_next_to_anchors
                    || node.identifier.coordinates(ts_kind) == primary_end.into()
                    || node.cost.match_run == 0
            })
            .filter_map(|node| {
                if let AlignmentCoordinates::Primary(primary) = node.identifier.coordinates(ts_kind)
                {
                    Some((primary, node.identifier.has_non_match(), node.cost.cost))
                } else {
                    None
                }
            })
            // The same condition as in `Context::is_target_identifier`.
            .filter(|(coordinates, has_non_match, _)| {
                *has_non_match
                    || !Self::requires_non_match(
                        enforcement_offset,
                        *coordinates == primary_end,
                        descendant_start,
                        Self::descendant_end(*coordinates, ts_kind),
                    )
            })
            .map(|(coordinates, _, cost)| (coordinates, cost)),
        );
    }
}
