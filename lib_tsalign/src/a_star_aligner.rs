use std::{
    fmt::{Debug, Display},
    time::Instant,
};

use alignment_geometry::AlignmentRange;
use alignment_result::{AlignmentResult, IAlignmentType};
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStar, AStarContext, AStarNode, AStarResult, cost::AStarCost};
use log::{debug, info};
use template_switch_distance::{
    context::Memory,
    strategies::{
        AlignmentStrategySelector, chaining::ChainingStrategy,
        primary_match::AllowPrimaryMatchStrategy, shortcut::NoShortcutStrategy,
        template_switch_count::TemplateSwitchCountStrategy,
        template_switch_min_length::TemplateSwitchMinLengthStrategy,
    },
};
use traitsequence::interface::Sequence;

use crate::{a_star_aligner::alignment_result::alignment::Alignment, config};

pub mod alignment_geometry;
pub mod alignment_result;
pub mod configurable_a_star_align;
pub mod gap_affine_edit_distance;
pub mod template_switch_distance;
#[cfg(test)]
mod tests;

pub trait AlignmentContext: AStarContext {
    type AlphabetType: Alphabet;

    type AlignmentType: IAlignmentType
        + Debug
        + From<<<Self as AStarContext>::Node as AStarNode>::EdgeType>;

    type SubsequenceType: GenomeSequence<Self::AlphabetType, Self::SubsequenceType> + ?Sized;

    fn reference(&self) -> &Self::SubsequenceType;

    fn query(&self) -> &Self::SubsequenceType;

    fn reference_name(&self) -> &str;

    fn query_name(&self) -> &str;

    fn range(&self) -> &AlignmentRange;
}

fn a_star_align<Context: AStarContext + AlignmentContext>(
    context: Context,
    debug_backtracking_identifiers: Vec<<Context::Node as AStarNode>::Identifier>,
) -> AlignmentResult<Context::AlignmentType, <<Context as AStarContext>::Node as AStarNode>::Cost>
where
    <Context::Node as AStarNode>::Identifier: Display,
    <Context::Node as AStarNode>::EdgeType: IAlignmentType + Display,
{
    info!("Aligning on subsequence {}", context.range());
    debug!("Is label setting: {}", context.is_label_setting());

    let start_time = Instant::now();

    // Perform forwards search.
    let mut a_star = AStar::<_>::new(context);
    a_star.initialise();
    let result = a_star.search();
    let has_target = matches!(result, AStarResult::FoundTarget { .. });

    for identifier in &debug_backtracking_identifiers {
        let mut debug_alignment = Vec::new();
        if let Some(backtracking_iterator) = a_star.backtrack_from(identifier) {
            for alignment_type in backtracking_iterator {
                if !alignment_type.is_internal() {
                    if let Some((count, previous_alignment_type)) = debug_alignment.last_mut() {
                        if alignment_type.is_repeated(previous_alignment_type) {
                            *count += 1;
                        } else {
                            debug_alignment.push((1, alignment_type));
                        }
                    } else {
                        debug_alignment.push((1, alignment_type));
                    }
                }
            }

            debug_alignment.reverse();
            let debug_alignment = Alignment::from(debug_alignment).cigar();
            println!("Alignment for identifier {identifier}: {debug_alignment}");
        }
    }

    let mut alignment = Vec::new();

    if has_target {
        // Backtrack.
        for alignment_type in a_star
            .backtrack()
            .map(<Context as AlignmentContext>::AlignmentType::from)
        {
            if !alignment_type.is_internal() {
                if let Some((count, previous_alignment_type)) = alignment.last_mut() {
                    if alignment_type.is_repeated(previous_alignment_type) {
                        *count += 1;
                    } else {
                        alignment.push((1, alignment_type));
                    }
                } else {
                    alignment.push((1, alignment_type));
                }
            }
        }

        alignment.reverse();
    }

    let end_time = Instant::now();
    let duration = (end_time - start_time).as_secs_f64();

    if has_target {
        AlignmentResult::new_with_target(
            alignment,
            a_star.context().reference(),
            a_star.context().query(),
            a_star.context().reference_name(),
            a_star.context().query_name(),
            a_star.context().range().reference_offset(),
            a_star.context().range().query_offset(),
            result.without_node_identifier(),
            duration,
            a_star.performance_counters().opened_nodes,
            a_star.performance_counters().closed_nodes,
            a_star.performance_counters().suboptimal_opened_nodes,
            a_star.context().reference().len(),
            a_star.context().query().len(),
        )
    } else {
        AlignmentResult::new_without_target(
            result.without_node_identifier(),
            a_star.context().reference(),
            a_star.context().query(),
            a_star.context().reference_name(),
            a_star.context().query_name(),
            a_star.context().range().reference_offset(),
            a_star.context().range().query_offset(),
            duration,
            a_star.performance_counters().opened_nodes,
            a_star.performance_counters().closed_nodes,
            a_star.performance_counters().suboptimal_opened_nodes,
            a_star.context().reference().len(),
            a_star.context().query().len(),
        )
    }
}

pub fn gap_affine_edit_distance_a_star_align<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    scoring_table: gap_affine_edit_distance::ScoringTable<Cost>,
) -> AlignmentResult<gap_affine_edit_distance::AlignmentType, Cost> {
    a_star_align(
        gap_affine_edit_distance::Context::new(reference, query, scoring_table),
        vec![],
    )
}

#[expect(clippy::too_many_arguments)]
pub fn template_switch_distance_a_star_align<
    Strategies: AlignmentStrategySelector<
            Shortcut = NoShortcutStrategy<<Strategies as AlignmentStrategySelector>::Cost>,
            PrimaryMatch = AllowPrimaryMatchStrategy,
        >,
    SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    reference_name: &str,
    query_name: &str,
    mut range: AlignmentRange,
    config: &config::TemplateSwitchConfig<
        Strategies::Alphabet,
        <Strategies as AlignmentStrategySelector>::Cost,
    >,
    cost_limit: Option<Strategies::Cost>,
    memory_limit: Option<usize>,
    force_label_correcting: bool,
    extend_beyond_range: bool,
    template_switch_count_memory: <Strategies::TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) -> AlignmentResult<template_switch_distance::AlignmentType, Strategies::Cost>
where
    Strategies::Cost: From<u64>,
{
    debug!("Reference sequence: {}", reference.as_string());
    debug!("Query sequence: {}", query.as_string());

    let memory = Memory {
        template_switch_min_length: <<Strategies as AlignmentStrategySelector>::TemplateSwitchMinLength as
            TemplateSwitchMinLengthStrategy<<Strategies as AlignmentStrategySelector>::Cost,
        >>::initialise_memory(reference, query, config),
        chaining: <<Strategies as AlignmentStrategySelector>::Chaining as ChainingStrategy<
            <Strategies as AlignmentStrategySelector>::Cost,
        >>::initialise_memory(reference, query, config, 20),
        template_switch_count: template_switch_count_memory,
        shortcut: (),
        primary_match: (),
    };

    info!("Calling aligner...");
    let mut result = a_star_align(
        template_switch_distance::Context::<SubsequenceType, Strategies>::new(
            reference,
            query,
            reference_name,
            query_name,
            range.clone(),
            config.clone(),
            memory,
            cost_limit,
            memory_limit,
            force_label_correcting,
        ),
        vec![],
    );
    info!("Main alignment finished");

    if extend_beyond_range {
        info!("Extending range");
        debug!("CIGAR before extending: {}", result.cigar());
        let extension_steps = result
            .extend_beyond_range_without_increasing_cost(reference, query, &mut range, config);
        info!(
            "Extended alignment {extension_steps} steps beyond the alignment range without increasing alignment costs"
        );
        info!("Alignment ranges after extension {range}");
    }

    info!("Extending template switches");
    result.compute_ts_equal_cost_ranges(reference, query, &range, config);
    result
}
