use std::time::Instant;

use self::{
    context::Memory,
    strategies::{
        AlignmentStrategySelector, chaining::ChainingStrategy,
        primary_match::AllowPrimaryMatchStrategy, shortcut::NoShortcutStrategy,
        template_switch_count::TemplateSwitchCountStrategy,
        template_switch_min_length::TemplateSwitchMinLengthStrategy,
    },
};
use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::{AStar, AStarContext, AStarResult, reset::Reset};
use log::{debug, info};

use crate::{
    aligner::{context::DynamicStrategies, ts_14_out_of_range::AdditionalExplicitTSMStartsAndEnds},
    alignment_geometry::AlignmentRange,
    alignment_result::{AlignmentResult, TSMUncertaintyRangeExtensionMode},
    config,
};

pub mod configurable;
pub mod context;
pub mod display;
mod identifier;
pub mod lower_bounds;
mod node;
pub mod strategies;
pub mod ts_14_out_of_range;

pub use context::Context;
pub use identifier::{
    GapType, Identifier, TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection,
};
pub use node::{Node, NodeData};

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
    range: AlignmentRange,
    additional_tsm_starts_and_ends: AdditionalExplicitTSMStartsAndEnds,
    config: &config::TemplateSwitchConfig<
        Strategies::Alphabet,
        <Strategies as AlignmentStrategySelector>::Cost,
    >,
    mut dynamic_strategies: DynamicStrategies,
    cost_limit: Option<Strategies::Cost>,
    memory_limit: Option<usize>,
    force_label_correcting: bool,
    extend_beyond_range: bool,
    tsm_uncertainty_range_extension_mode: TSMUncertaintyRangeExtensionMode,
    template_switch_count_memory: <Strategies::TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) -> AlignmentResult<Strategies::Cost>
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
    dynamic_strategies.reset();

    info!("Calling aligner...");
    let context = Context::<SubsequenceType, Strategies>::new(
        reference,
        query,
        reference_name,
        query_name,
        range,
        additional_tsm_starts_and_ends,
        config.clone(),
        memory,
        dynamic_strategies,
        cost_limit,
        memory_limit,
        force_label_correcting,
    );
    info!("Aligning on subsequence {}", context.range());
    debug!("Is label setting: {}", context.is_label_setting());

    let start_time = Instant::now();

    // Perform forwards search.
    let mut a_star = AStar::<_>::new(context);
    a_star.initialise();
    let search_result = a_star.search();
    let has_target = matches!(search_result, AStarResult::FoundTarget { .. });

    let mut alignment = Vec::new();
    let mut first_alignment = None;
    let mut last_alignment = None;

    if has_target {
        // Backtrack.
        for alignment_type in a_star.backtrack() {
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

            if last_alignment.is_none() {
                last_alignment = Some(alignment_type);
            }
            first_alignment = Some(alignment_type);
        }

        alignment.reverse();
    }

    let mut alignment_range = a_star.context().range().clone();
    if let Some(alignment_type) = first_alignment {
        if let Some(alternative_start) = alignment_type.alternative_start() {
            assert!(alternative_start.reference() <= alignment_range.reference_offset());
            assert!(alternative_start.query() <= alignment_range.query_offset());
            alignment_range = alignment_range.with_offset(alternative_start);
        }
    }
    if let Some(alignment_type) = last_alignment {
        if let Some(alternative_end) = alignment_type.alternative_end() {
            assert!(alternative_end.reference() >= alignment_range.reference_limit());
            assert!(alternative_end.query() >= alignment_range.query_limit());
            alignment_range = alignment_range.with_limit(alternative_end);
        }
    }

    let end_time = Instant::now();
    let duration = (end_time - start_time).as_secs_f64();

    let mut result = if has_target {
        AlignmentResult::new_with_target(
            alignment,
            a_star.context().reference(),
            a_star.context().query(),
            a_star.context().reference_name(),
            a_star.context().query_name(),
            alignment_range,
            search_result.without_node_identifier(),
            duration,
            a_star.performance_counters().opened_nodes,
            a_star.performance_counters().closed_nodes,
            a_star.performance_counters().suboptimal_opened_nodes,
            a_star.context().reference().len(),
            a_star.context().query().len(),
        )
    } else {
        AlignmentResult::new_without_target(
            search_result.without_node_identifier(),
            a_star.context().reference(),
            a_star.context().query(),
            a_star.context().reference_name(),
            a_star.context().query_name(),
            alignment_range,
            duration,
            a_star.performance_counters().opened_nodes,
            a_star.performance_counters().closed_nodes,
            a_star.performance_counters().suboptimal_opened_nodes,
            a_star.context().reference().len(),
            a_star.context().query().len(),
        )
    };
    info!("Main alignment finished");

    let mut range = result.alignment_range();

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

    match tsm_uncertainty_range_extension_mode {
        TSMUncertaintyRangeExtensionMode::None => {
            debug!("TSM uncertainty range extension is disabled.")
        }

        TSMUncertaintyRangeExtensionMode::EqualCost => {
            info!("Extending TSM uncertainty ranges without increasing cost");
            result.compute_ts_uncertainty_ranges(
                reference,
                query,
                &range,
                config,
                tsm_uncertainty_range_extension_mode,
            );
        }

        TSMUncertaintyRangeExtensionMode::EqualCostIgnoreGeometry => {
            info!(
                "Extending TSM uncertainty ranges without increasing cost but ignoring cost increases due to TSM geometry",
            );
            result.compute_ts_uncertainty_ranges(
                reference,
                query,
                &range,
                config,
                tsm_uncertainty_range_extension_mode,
            );
        }
    }

    result
}
