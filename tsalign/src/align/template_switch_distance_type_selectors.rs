use std::fmt::Debug;

use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::{
    a_star_aligner::{
        alignment_geometry::AlignmentRange,
        template_switch_distance::{
            context::DynamicStrategies,
            strategies::{
                AlignmentStrategySelection,
                chaining::{ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy},
                descendant::{
                    AnyTemplateSwitchDescendantStrategy, OnlyEqualTemplateSwitchDescendantStrategy,
                    TemplateSwitchDescendantStrategy,
                },
                node_ord::{AntiDiagonalNodeOrdStrategy, NodeOrdStrategy},
                primary_match::AllowPrimaryMatchStrategy,
                primary_range::NoPrunePrimaryRangeStrategy,
                secondary_deletion::AllowSecondaryDeletionStrategy,
                shortcut::NoShortcutStrategy,
                template_switch_count::{
                    MaxTemplateSwitchCountStrategy, NoTemplateSwitchCountStrategy,
                    TemplateSwitchCountStrategy,
                },
                template_switch_min_length::{
                    LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                    PreprocessedLookaheadTemplateSwitchMinLengthStrategy,
                    PreprocessedTemplateSwitchMinLengthStrategy, TemplateSwitchMinLengthStrategy,
                },
                template_switch_total_length::{
                    MaxTemplateSwitchTotalLengthStrategy, NoTemplateSwitchTotalLengthStrategy,
                    TemplateSwitchTotalLengthStrategy,
                },
            },
        },
        template_switch_distance_a_star_align,
    },
    costs::U64Cost,
};
use log::info;

use super::Cli;
use crate::util::load_tsa_config;

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchNodeOrdStrategySelector {
    // CostOnly,
    AntiDiagonal,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchMinLengthStrategySelector {
    None,
    Lookahead,
    /// Check for perfect matches of minimum length and compute a lower bound for imperfect matches.
    PreprocessPrice,
    /// Check for perfect matches of minimum length and filter out imperfect matches.
    PreprocessFilter,
    /// Check for perfect matches of minimum length and cache the costs of imperfect matches.
    PreprocessLookahead,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchChainingStrategySelector {
    None,
    // PrecomputeOnly,
    LowerBound,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchTotalLengthStrategySelector {
    None,
    Maximise,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchDescendantStrategySelector {
    AllowAny,
    AllowOnlyAllEqual,
}

pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    align_a_star_template_switch_distance_select_node_ord_strategy(
        cli,
        reference,
        query,
        range,
        reference_name,
        query_name,
    );
}

fn align_a_star_template_switch_distance_select_node_ord_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_node_ord_strategy {
        /*TemplateSwitchNodeOrdStrategySelector::CostOnly => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                CostOnlyNodeOrdStrategy,
            >(cli, reference, query, reference_name, query_name)
            unimplemented!("The other option appears to always be better.");
        }*/
        TemplateSwitchNodeOrdStrategySelector::AntiDiagonal => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                AntiDiagonalNodeOrdStrategy,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_distance_select_template_switch_min_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_min_length_strategy {
        TemplateSwitchMinLengthStrategySelector::None => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                NoTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::Lookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::PreprocessPrice => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                PreprocessedTemplateSwitchMinLengthStrategy<false, U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::PreprocessFilter => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                PreprocessedTemplateSwitchMinLengthStrategy<true, U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::PreprocessLookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                PreprocessedLookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_select_chaining_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_chaining_strategy {
        TemplateSwitchChainingStrategySelector::None => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                NoChainingStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        /*TemplateSwitchChainingStrategySelector::PrecomputeOnly => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
            unimplemented!("No reason to precompute without using the information.");
        }*/
        TemplateSwitchChainingStrategySelector::LowerBound => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                LowerBoundChainingStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_select_no_ts_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) {
    if cli.no_ts {
        align_a_star_template_switch_select_template_switch_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            MaxTemplateSwitchCountStrategy,
        >(cli, reference, query, range, reference_name, query_name, 0)
    } else {
        align_a_star_template_switch_select_template_switch_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
        >(cli, reference, query, range, reference_name, query_name, ())
    }
}

fn align_a_star_template_switch_select_template_switch_total_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) {
    match cli.ts_total_length_strategy {
        TemplateSwitchTotalLengthStrategySelector::None => {
            align_a_star_template_switch_select_template_switch_descendant_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                NoTemplateSwitchTotalLengthStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
        TemplateSwitchTotalLengthStrategySelector::Maximise => {
            align_a_star_template_switch_select_template_switch_descendant_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                MaxTemplateSwitchTotalLengthStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
    }
}

fn align_a_star_template_switch_select_template_switch_descendant_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    TemplateSwitchTotalLength: TemplateSwitchTotalLengthStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) {
    match cli.ts_descendant_strategy {
        TemplateSwitchDescendantStrategySelector::AllowAny => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                TemplateSwitchTotalLength,
                AnyTemplateSwitchDescendantStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
        TemplateSwitchDescendantStrategySelector::AllowOnlyAllEqual => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                TemplateSwitchTotalLength,
                OnlyEqualTemplateSwitchDescendantStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
    }
}

fn align_a_star_template_switch_distance_call<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    TemplateSwitchTotalLength: TemplateSwitchTotalLengthStrategy,
    TemplateSwitchDescendant: TemplateSwitchDescendantStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) {
    let costs = load_tsa_config(&cli.configuration_directory).unwrap();
    let dynamic_strategies = DynamicStrategies {};

    let alignment = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            AlphabetType,
            U64Cost,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            TemplateSwitchCount,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            TemplateSwitchTotalLength,
            TemplateSwitchDescendant,
        >,
        _,
    >(
        reference,
        query,
        reference_name,
        query_name,
        range,
        &costs,
        dynamic_strategies,
        cli.cost_limit,
        cli.memory_limit,
        cli.force_label_correcting,
        !cli.dont_extend_beyond_range,
        template_switch_count_memory,
    );
    info!("Finished aligning");

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{alignment}");
}
