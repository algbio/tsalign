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
        // Enfore non-match if there is a gap between the anchors in the descendant.
        // This is to match the lower-bound computation.
        // It also discourages chains to deviate from the alignment geometry boundaries.
        let descendant_end = match start.ts_kind().descendant {
            TsDescendant::Seq1 => end.a(),
            TsDescendant::Seq2 => end.b(),
        };
        let descendant_start = start.descendant();
        let enforce_non_match = descendant_end != start.descendant();

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

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            true,
            self.max_match_run,
            self.enforcement_offset,
            self.enforce_non_match_next_to_anchors(),
        );
        let mut dijkstra = new_dijkstra(&mut self.datastructures, context);
        search_until_cost_limit(&mut dijkstra, cost_limit);

        let descendant_start = start.descendant();
        let ts_kind = start.ts_kind();
        Self::fill_additional_targets(
            dijkstra.prune_list(),
            self.enforce_non_match_next_to_anchors(),
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

    fn fill_additional_targets(
        prune_list: &PruneList<Node<Cost>>,
        enforce_non_match_next_to_anchors: bool,
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
            .filter(|(coordinates, has_non_match, _)| {
                *has_non_match
                    != (descendant_start
                        == match ts_kind.descendant {
                            TsDescendant::Seq1 => coordinates.a(),
                            TsDescendant::Seq2 => coordinates.b(),
                        })
            })
            .map(|(coordinates, _, cost)| (coordinates, cost)),
        );
    }
}
