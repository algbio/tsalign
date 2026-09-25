use generic_a_star::cost::AStarCost;

use crate::{
    alignment::Alignment, anchors::Anchors, chain_align::chainer::Identifier,
    chaining_cost_function::ChainingCostFunction,
};

pub mod exact;
pub mod inexact;

pub trait ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    fn evaluate_chain(
        &mut self,
        anchors: &Anchors<Cost>,
        chain: &[Identifier],
        anchor_k: u32,
        max_anchor_mutations: u8,
        chaining_cost_function: &mut ChainingCostFunction<Cost>,
        final_evaluation: bool,
    ) -> (Cost, Vec<Alignment>);

    fn total_gaps(&self) -> u64;

    fn total_gap_fillings(&self) -> u64;

    fn total_redundant_gap_fillings(&self) -> u64;

    fn gap_fill_alignments_per_chain(&self) -> &[u32];
}
