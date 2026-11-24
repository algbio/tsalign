use std::fmt::Display;

use generic_a_star::cost::AStarCost;

use crate::{
    alignment::{Alignment, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::{
        chainer::Identifier,
        evaluation::{
            max_match_run::MaxMatchRunChainEvaluator,
            windowed_min_mutations::WindowedMinMutationsChainEvaluator,
        },
    },
    chaining_cost_function::ChainingCostFunction,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::AlignmentCosts,
};

pub mod max_match_run;
pub mod windowed_min_mutations;

pub trait ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    fn evaluate_chain(
        &mut self,
        anchors: &Anchors<Cost>,
        chain: &[Identifier],
        config: &ChainingLowerBoundConfig,
        chaining_cost_function: &mut ChainingCostFunction<Cost>,
        final_evaluation: bool,
    ) -> (Cost, Vec<Alignment>);

    fn total_gaps(&self) -> u64;

    fn total_gap_fillings(&self) -> u64;

    fn total_redundant_gap_fillings(&self) -> u64;

    fn gap_fill_alignments_per_chain(&self) -> &[u32];
}

/// Checks the cost of the re-alignment of an anchor in the final evaluation of a chain against the cost of the anchor.
///
/// The final evaluation aligns each anchor of the chain again, without any restrictions.
/// The cost of an exact anchor is zero, and so is its re-alignment.
/// But the cost of an inexact anchor is the minimum cost of its alignments with at most `max_anchor_mutations` mutations,
/// while the re-alignment may use more mutations, and hence may be cheaper.
/// The chain cost uses the anchor cost, while the cost of the reported alignment uses the re-alignment.
/// This is one of the differences between the chain cost and the cost of the reported alignment,
/// see [issue #221](https://github.com/algbio/tsalign/issues/221).
fn check_anchor_alignment_cost<Cost: AStarCost>(
    anchor_alignment_cost: Cost,
    anchor_cost: Cost,
    anchor: impl Display,
) {
    assert!(
        anchor_alignment_cost <= anchor_cost,
        "The re-alignment of anchor {anchor} costs {anchor_alignment_cost}, which is more than the anchor cost {anchor_cost}",
    );
    if anchor_cost.is_zero() {
        assert!(anchor_alignment_cost.is_zero());
    }
}

pub fn build_chain_evaluator<
    'sequences: 'result,
    'alignment_costs: 'result,
    'rc_fn: 'result,
    'result,
    Cost: AStarCost,
>(
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    config: &ChainingLowerBoundConfig,
) -> Box<dyn 'result + ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost>> {
    if config.max_anchor_mutations == 0 {
        Box::new(MaxMatchRunChainEvaluator::new(
            sequences,
            alignment_costs,
            rc_fn,
            config.anchor_k - 1,
            0,
        ))
    } else {
        let anchor_k = config
            .anchor_k
            .try_into()
            .expect("anchor_k must be <= 255.");

        match config.inexact_lower_bound {
            InexactLowerBoundKind::WindowedMinMutations => {
                Box::new(WindowedMinMutationsChainEvaluator::new(
                    sequences,
                    alignment_costs,
                    rc_fn,
                    anchor_k,
                    config.max_anchor_mutations,
                ))
            }

            InexactLowerBoundKind::MaxMatchRun => Box::new(MaxMatchRunChainEvaluator::new(
                sequences,
                alignment_costs,
                rc_fn,
                (anchor_k - 1 - config.max_anchor_mutations).into(),
                anchor_k,
            )),
        }
    }
}
