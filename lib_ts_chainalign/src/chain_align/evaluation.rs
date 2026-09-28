use generic_a_star::cost::AStarCost;

use crate::{
    alignment::{Alignment, sequences::AlignmentSequences},
    alignment_history::history_vec::{AlignmentHistory, UnsignedIntAlignmentHistoryVec},
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
                if UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
                    anchor_k,
                    config.max_anchor_mutations,
                ) {
                    Box::new(WindowedMinMutationsChainEvaluator::<
                        _,
                        UnsignedIntAlignmentHistoryVec<u8>,
                    >::new(
                        sequences,
                        alignment_costs,
                        rc_fn,
                        anchor_k,
                        config.max_anchor_mutations,
                    ))
                } else if UnsignedIntAlignmentHistoryVec::<u16>::has_enough_capacity(
                    anchor_k,
                    config.max_anchor_mutations,
                ) {
                    Box::new(WindowedMinMutationsChainEvaluator::<
                        _,
                        UnsignedIntAlignmentHistoryVec<u16>,
                    >::new(
                        sequences,
                        alignment_costs,
                        rc_fn,
                        anchor_k,
                        config.max_anchor_mutations,
                    ))
                } else if UnsignedIntAlignmentHistoryVec::<u32>::has_enough_capacity(
                    anchor_k,
                    config.max_anchor_mutations,
                ) {
                    Box::new(WindowedMinMutationsChainEvaluator::<
                        _,
                        UnsignedIntAlignmentHistoryVec<u32>,
                    >::new(
                        sequences,
                        alignment_costs,
                        rc_fn,
                        anchor_k,
                        config.max_anchor_mutations,
                    ))
                } else if UnsignedIntAlignmentHistoryVec::<u64>::has_enough_capacity(
                    anchor_k,
                    config.max_anchor_mutations,
                ) {
                    Box::new(WindowedMinMutationsChainEvaluator::<
                        _,
                        UnsignedIntAlignmentHistoryVec<u64>,
                    >::new(
                        sequences,
                        alignment_costs,
                        rc_fn,
                        anchor_k,
                        config.max_anchor_mutations,
                    ))
                } else if UnsignedIntAlignmentHistoryVec::<u128>::has_enough_capacity(
                    anchor_k,
                    config.max_anchor_mutations,
                ) {
                    Box::new(WindowedMinMutationsChainEvaluator::<
                        _,
                        UnsignedIntAlignmentHistoryVec<u128>,
                    >::new(
                        sequences,
                        alignment_costs,
                        rc_fn,
                        anchor_k,
                        config.max_anchor_mutations,
                    ))
                } else {
                    panic!(
                        "This combination of k and max_mismatches exceeds the maximum length of the alignment history vector that can be stored in a u128.",
                    );
                }
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
