use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{
        Alignment,
        coordinates::AlignmentCoordinates,
        sequences::AlignmentSequences,
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::primary::PrimaryAnchor,
    costs::AlignmentCosts,
    exact_chaining::ts_34_jump::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost>
    Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
    ) -> Self {
        Self {
            a_star_buffers: Some(Default::default()),
            sequences,
            alignment_costs,
            rc_fn,
            max_match_run,
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
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
    ) -> (Cost, Alignment) {
        assert!(start.is_secondary());
        assert!(end.is_primary());

        // Enfore non-match if there is a gap between the anchors in the descendant.
        // This is to match the lower-bound computation.
        // It also discourages chains to deviate from the alignment geometry boundaries.
        let descendant_end = match start.ts_kind().unwrap().descendant {
            TsDescendant::Seq1 => end.primary_ordinate_a(),
            TsDescendant::Seq2 => end.primary_ordinate_b(),
        }
        .unwrap();
        let descendant_start = start.secondary_ordinate_descendant().unwrap();
        // Towards the alignment end, a match run after the jump is not covered by any anchor
        let enforce_non_match = end != self.sequences.primary_end()
            && descendant_end != start.secondary_ordinate_descendant().unwrap();

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
            self.max_match_run,
        );
        let mut a_star = AStar::<_>::new_with_buffers(context, self.a_star_buffers.take().unwrap());

        a_star.initialise();
        let (cost, alignment) = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => {
                let cost = cost.0;
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
            descendant_start,
            start.ts_kind().unwrap(),
            self.sequences.primary_end(),
            additional_primary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output list.
    /// Note that the output list may contain duplicate anchors with different cost.
    pub fn align_until_cost_limit(
        &mut self,
        start: AlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
    ) {
        assert!(start.is_secondary());
        let end = self.sequences.primary_end();

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            true,
            self.max_match_run,
        );
        let mut a_star = AStar::new_with_buffers(context, self.a_star_buffers.take().unwrap());
        a_star.initialise();
        a_star.search_until_with_target_policy(|_, node| node.cost > cost_limit, true);

        let descendant_start = start.secondary_ordinate_descendant().unwrap();
        let ts_kind = start.ts_kind().unwrap();
        Self::fill_additional_targets(
            &a_star,
            descendant_start,
            ts_kind,
            end,
            additional_primary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());
    }

    fn fill_additional_targets(
        a_star: &AStar<Context<Cost>>,
        descendant_start: usize,
        ts_kind: TsKind,
        primary_end: AlignmentCoordinates,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
    ) {
        additional_primary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| {
                    node.identifier.coordinates().is_primary()
                        && (node.identifier.has_non_match()
                            != (descendant_start
                                == match ts_kind.descendant {
                                    TsDescendant::Seq1 => {
                                        node.identifier.coordinates().primary_ordinate_a().unwrap()
                                    }
                                    TsDescendant::Seq2 => {
                                        node.identifier.coordinates().primary_ordinate_b().unwrap()
                                    }
                                })
                            || (node.identifier.coordinates() == primary_end
                                && !node.identifier.has_non_match()))
                })
                .map(|node| {
                    (
                        PrimaryAnchor::new_from_start(&node.identifier.coordinates()),
                        node.cost,
                    )
                }),
        );
    }
}
