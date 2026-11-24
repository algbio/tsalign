use std::{fmt::Debug, ops::Range, path::PathBuf};

use anyhow::{Context, Result, anyhow, bail, ensure};
use astar::{
    TemplateSwitchNodeOrdStrategySelector, align_a_star_template_switch_distance, configure_aligner,
};
use clap::{
    ArgAction, ArgMatches, Args, Parser, ValueEnum,
    builder::{BoolishValueParser, TypedValueParser},
};
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_ts_chainalign::chain_align::performance_parameters::{
    AnchorPruning, ChainingClosedList, ChainingOpenList, InexactPruningRules,
};
use lib_tsalign::{
    aligner::{
        configurable::{
            ChainingStrategySelector, DescendantStrategySelector, MinLengthStrategySelector,
            TotalLengthStrategySelector,
        },
        ts_14_out_of_range::{AdditionalExplicitTSMStartsAndEnds, Ts14OutOfRangeStrategy},
    },
    alignment_geometry::{AlignmentCoordinates, AlignmentRange},
    alignment_result::TSMUncertaintyRangeExtensionMode,
    costs::U64Cost,
    error::{SequenceRole, sequence_from_ascii},
};
use log::{LevelFilter, debug, info};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

use crate::{
    align::{
        chain::align_a_star_chain_ts,
        fasta_parser::{parse_pair_fasta_file, parse_single_fasta_file},
    },
    alphabet::{InputAlphabet, WithAlphabet},
    util::{is_explicit, reject_inexact_anchor_arguments_for_exact_anchors},
};

mod astar;
mod chain;
mod fasta_parser;

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    #[command(flatten)]
    input: CliInput,

    /// The file to store the alignment statistics in toml format.
    #[clap(long, short = 'o')]
    output: Option<PathBuf>,

    #[command(flatten)]
    common: CommonArgs,

    #[command(
        flatten,
        next_help_heading = "Classic aligner options (a-star-template-switch)"
    )]
    a_star: AStarArgs,

    #[command(flatten, next_help_heading = "tschainalign options (a-star-chain-ts)")]
    chain: ChainArgs,
}

// Options that apply to all alignment methods.
// The argument groups have no doc comments, because clap would use them as the description of the command.
#[derive(Args)]
struct CommonArgs {
    /// The alphabet present in the input files.
    ///
    /// This must also match the alphabet used in the config.
    #[clap(long, short = 'a', default_value = "dna-n")]
    alphabet: InputAlphabet,

    /// A string of (ASCII) characters that should be skipped in the input fasta.
    ///
    /// For example, `-` characters caused by alignment hints can be skipped this way.
    #[clap(long, default_value = "")]
    skip_characters: String,

    /// A directory containing the configuration files.
    ///
    /// See the README for its layout.
    /// Default: `sample_tsa_config` for `a-star-template-switch`;
    /// for `a-star-chain-ts`, the built-in default config of tschainalign
    /// (equal to `test_files/config/chainalign`), since tschainalign does not support `sample_tsa_config`.
    #[clap(long, short = 'c')]
    configuration_directory: Option<PathBuf>,

    #[clap(long, default_value = "a-star-template-switch")]
    alignment_method: AlignmentMethod,

    /// First character in the reference to start the alignment from.
    ///
    /// Skipped characters are ignored for computing this index.
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    reference_offset: Option<usize>,

    /// First character in the query to start the alignment from.
    ///
    /// Skipped characters are ignored for computing this index.
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    query_offset: Option<usize>,

    /// First character after the last character in the reference to end the alignment at.
    ///
    /// Skipped characters are ignored for computing this index.
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    reference_limit: Option<usize>,

    /// First character after the last character in the query to end the alignment at.
    ///
    /// Skipped characters are ignored for computing this index.
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    query_limit: Option<usize>,

    /// Ranges for computing the alignment.
    ///
    /// Skipped characters are ignored for the range indices.
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    rq_ranges: Option<String>,

    /// Use the ranges for alignment that are embedded in the fasta file(s).
    ///
    /// The start and end of the range must be delimited by '|' characters in the fasta file.
    /// Both reference and query must contain exactly two delimiters each.
    ///
    /// Template switch inners can still align to the full sequence.
    #[clap(long)]
    use_embedded_rq_ranges: bool,

    /// An inclusive cost limit for the alignment.
    ///
    /// If there is no alignment with at most that cost, the aligner will abort without result.
    ///
    /// For tschainalign (a-star-chain-ts), the limit applies to the chain cost, which tschainalign minimises.
    /// For inexact anchors, the cost of the alignment may be lower than its chain cost.
    #[clap(long)]
    cost_limit: Option<U64Cost>,
}

// Options of the classic aligner (`a-star-template-switch`).
#[derive(Args)]
struct AStarArgs {
    /// Has no effect, the only implemented node ordering is anti-diagonal.
    #[clap(long, default_value = "anti-diagonal", hide = true)]
    #[expect(
        dead_code,
        reason = "kept as a no-op so that existing command lines keep working"
    )]
    ts_node_ord_strategy: TemplateSwitchNodeOrdStrategySelector,

    #[clap(long, default_value = "lookahead")]
    ts_min_length_strategy: MinLengthStrategySelector,

    /// The chaining strategy of the classic aligner.
    ///
    /// The `lower-bound` strategy is currently disabled (issue #248).
    #[clap(long, default_value = "none")]
    ts_chaining_strategy: ChainingStrategySelector,

    /// If set to maximise, the total length of template switches becomes a secondary search criterion.
    ///
    /// As a result, instead of reporting an arbitrary alignment with optimal cost,
    /// tsalign selects an alignment with maximal total TS length out of all the alignments with optimal cost.
    #[clap(long, default_value = "maximise")]
    ts_total_length_strategy: TotalLengthStrategySelector,

    /// If set to allow-only-all-equal, then only solutions are considered where all descendants of all TSMs are equal.
    #[clap(long, default_value = "allow-any")]
    ts_descendant_strategy: DescendantStrategySelector,

    /// If set, then TSMs can start or end outside of the specified alignment ranges.
    /// In this case the alignment would start or end outside of the specified alignment ranges,
    /// but no other alignment operations outside of the ranges are allowed.
    ///
    /// The TSMs that start or end outside of the specified alignment ranges can only start or end at coordinates that correspond to an indel-free alignment of the sequences outside of the specified ranges.
    /// **Every skip character is treated as an indel.**
    /// This feature is intended for sequences that stem from an existing alignment, where gaps are marked as skip characters.
    ///
    /// In case the alignment ranges are not set or set to the whole sequences, this option has no effect.
    #[clap(
        long,
        action = ArgAction::SetTrue,
        value_parser = BoolishValueParser::new().map(Ts14OutOfRangeStrategy::from),
    )]
    allow_ts_14_out_of_range: Ts14OutOfRangeStrategy,

    /// Determines how the TSM uncertainty ranges are extended.
    #[clap(long, default_value = "equal-cost-ignore-geometry")]
    ts_uncertainty_range_extension_mode: TSMUncertaintyRangeExtensionMode,

    /// If set, template switches are not allowed.
    ///
    /// Use this to compare a template switch alignment against an alignment with out template switches.
    #[clap(long)]
    no_ts: bool,

    /// An approximate memory limit in bytes for the aligner.
    ///
    /// If it is exceeded, then the aligner will abort without result.
    /// Not supported by tschainalign (a-star-chain-ts).
    #[clap(long)]
    memory_limit: Option<usize>,

    /// Force the search to run in label-correcting mode.
    ///
    /// This is for debug purposes only.
    #[clap(long)]
    force_label_correcting: bool,

    /// If set, the aligner will not attempt to extend the alignment beyond the specified range
    /// without increasing cost.
    #[clap(long)]
    dont_extend_beyond_range: bool,
}

// Options of tschainalign (`a-star-chain-ts`).
#[derive(Args)]
struct ChainArgs {
    /// The directory in which preprocessed data is stored.
    #[clap(long)]
    cache_directory: Option<PathBuf>,

    /// If set, then the aligner performs no preprocessing, but instead crashes if the cached file is missing.
    #[clap(long)]
    force_no_preprocessing: bool,

    /// k-mer size for tschainalign.
    ///
    /// If it is not specified, it is inferred from the sequence lengths.
    #[clap(short)]
    k: Option<u32>,

    /// The maximum number of mutations allowed in an anchor.
    ///
    /// Higher values produce a more accurate alignment at the cost of increased runtime.
    #[clap(long, default_value = "0")]
    max_anchor_mutations: u8,

    /// The lower bound to use for inexact chaining.
    #[clap(long, default_value = "max-match-run")]
    inexact_chaining_lower_bound: InexactLowerBound,

    /// The maximum amount of successors to generate while processing a node during chaining.
    /// Can be tuned to optimise performance.
    #[clap(long, default_value = "10")]
    max_chaining_successors: usize,

    /// The maximum cost until which the cost function is initialised exactly.
    ///
    /// Setting this to a higher value will increase the time required to initialise the cost function, but decrease the amount of chains computed.
    #[clap(long, default_value = "0")]
    max_exact_cost_function_cost: u32,

    /// The closed list type to use for chaining.
    #[clap(long, default_value = "special")]
    chaining_closed_list: ChainingClosedList,

    /// The open list type to use for chaining.
    #[clap(long, default_value = "linear-heap")]
    chaining_open_list: ChainingOpenList,

    /// Which inexact anchors to remove before chaining (only for --max-anchor-mutations >= 1).
    ///
    /// `exact` removes dominated primary anchors only where this provably does not change the cost of the optimal chain
    /// (one anchor mutation). `inexact` removes primary and secondary anchors with heuristic rules for any number of
    /// anchor mutations, which may increase the cost and lose template switches (see --inexact-pruning-rules).
    /// `none` keeps all anchors, which is useful to compare against the chaining with all anchors.
    #[clap(long, default_value = "exact")]
    anchor_pruning: AnchorPruning,

    /// Deprecated alias of `--anchor-pruning none`.
    #[clap(long, conflicts_with = "anchor_pruning")]
    no_anchor_domination_pruning: bool,

    /// The rules of `--anchor-pruning inexact`, as a comma-separated list [default: all].
    ///
    /// Rules: primary-drop-guards, primary-prefix-domination, primary-suffix-domination, primary-edge-substitution,
    /// primary-exact-overlap, and the same with the prefix `secondary-`.
    /// Keywords: `all`, `none`, `primary` (all primary rules), `secondary` (all secondary rules).
    /// An item prefixed with `-` removes rules again, e.g. `all,-secondary-exact-overlap`.
    /// Requires `--anchor-pruning inexact`.
    #[clap(long, allow_hyphen_values = true)]
    inexact_pruning_rules: Option<InexactPruningRules>,

    /// The maximum distance between the diagonal of an exact anchor and the diagonal range of an inexact anchor
    /// for the `*-exact-overlap` rules of `--anchor-pruning inexact` [default: 1].
    ///
    /// Requires `--anchor-pruning inexact`.
    #[clap(long)]
    exact_overlap_diagonal_distance: Option<usize>,

    /// The maximum number of chaining rounds.
    ///
    /// Each round searches the cheapest chain with regard to the current lower bounds of the gap costs,
    /// and then aligns its gaps exactly, which improves the lower bounds.
    /// If the optimal chain is not found within this many rounds, the aligner will abort without result,
    /// but reports a lower bound of the chain cost.
    #[clap(long)]
    iteration_limit: Option<u64>,
}

#[derive(Args)]
struct CliInput {
    #[clap(flatten)]
    separate_input: Option<CliSeparateInput>,

    #[clap(flatten)]
    pair_input: Option<CliPairInput>,
}

#[derive(Args)]
#[group(multiple = true)]
struct CliSeparateInput {
    /// The path to the reference fasta file.
    #[clap(long, short = 'r', required = false, requires = "query")]
    reference: PathBuf,

    /// The path to the query fasta file.
    #[clap(long, short = 'q', required = false, requires = "reference")]
    query: PathBuf,
}

#[derive(Args)]
struct CliPairInput {
    /// The path to a fasta file containing both the reference and the query.
    #[clap(long, short = 'p', required = false, conflicts_with_all = ["reference", "query"])]
    pair_fasta: PathBuf,
}

#[derive(Clone, PartialEq, Eq, ValueEnum)]
enum AlignmentMethod {
    AStarTemplateSwitch,
    AStarChainTS,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum InexactLowerBound {
    /// Compute a lower bound assuming that in every window of size `anchor_k` there have to be at least `max_anchor_mutations + 1` mutations.
    ///
    /// Currently disabled in the CLI because it is broken.
    #[value(skip)]
    WindowedMinMutations,
    /// Compute a lower bound assuming that there is no match run longer than `anchor_k - 1 - max_anchor_mutations`,
    /// except for alignments that are shorter than `anchor_k` on both sequences.
    MaxMatchRun,
}

/// Arguments that only the classic aligner (`a-star-template-switch`) uses: pairs of argument id and argument.
///
/// `--cost-limit` is a common argument, since both methods support it.
/// `--ts-node-ord-strategy` is not listed, since it is a hidden no-op for compatibility.
const A_STAR_ONLY_ARGUMENTS: &[(&str, &str)] = &[
    ("ts_min_length_strategy", "--ts-min-length-strategy"),
    ("ts_chaining_strategy", "--ts-chaining-strategy"),
    ("ts_total_length_strategy", "--ts-total-length-strategy"),
    ("ts_descendant_strategy", "--ts-descendant-strategy"),
    ("allow_ts_14_out_of_range", "--allow-ts-14-out-of-range"),
    (
        "ts_uncertainty_range_extension_mode",
        "--ts-uncertainty-range-extension-mode",
    ),
    ("no_ts", "--no-ts"),
    ("force_label_correcting", "--force-label-correcting"),
    ("dont_extend_beyond_range", "--dont-extend-beyond-range"),
    ("memory_limit", "--memory-limit"),
];

/// Arguments that only tschainalign (`a-star-chain-ts`) uses: pairs of argument id and argument.
const CHAIN_ONLY_ARGUMENTS: &[(&str, &str)] = &[
    ("cache_directory", "--cache-directory"),
    ("force_no_preprocessing", "--force-no-preprocessing"),
    ("k", "-k"),
    ("max_anchor_mutations", "--max-anchor-mutations"),
    (
        "inexact_chaining_lower_bound",
        "--inexact-chaining-lower-bound",
    ),
    ("max_chaining_successors", "--max-chaining-successors"),
    (
        "max_exact_cost_function_cost",
        "--max-exact-cost-function-cost",
    ),
    ("chaining_closed_list", "--chaining-closed-list"),
    ("chaining_open_list", "--chaining-open-list"),
    ("iteration_limit", "--iteration-limit"),
    (
        "no_anchor_domination_pruning",
        "--no-anchor-domination-pruning",
    ),
    ("anchor_pruning", "--anchor-pruning"),
    ("inexact_pruning_rules", "--inexact-pruning-rules"),
    (
        "exact_overlap_diagonal_distance",
        "--exact-overlap-diagonal-distance",
    ),
];

/// Arguments of the classic aligner that have no effect if template switches are disabled with `--no-ts`.
const TS_ONLY_ARGUMENTS: &[(&str, &str)] = &[
    ("ts_min_length_strategy", "--ts-min-length-strategy"),
    ("ts_chaining_strategy", "--ts-chaining-strategy"),
    ("ts_total_length_strategy", "--ts-total-length-strategy"),
    ("ts_descendant_strategy", "--ts-descendant-strategy"),
    ("allow_ts_14_out_of_range", "--allow-ts-14-out-of-range"),
    (
        "ts_uncertainty_range_extension_mode",
        "--ts-uncertainty-range-extension-mode",
    ),
];

/// Arguments of tschainalign that only affect inexact anchors (`--max-anchor-mutations` > 0).
const INEXACT_ANCHOR_ONLY_ARGUMENTS: &[(&str, &str)] = &[
    (
        "inexact_chaining_lower_bound",
        "--inexact-chaining-lower-bound",
    ),
    (
        "no_anchor_domination_pruning",
        "--no-anchor-domination-pruning",
    ),
    ("anchor_pruning", "--anchor-pruning"),
    ("inexact_pruning_rules", "--inexact-pruning-rules"),
    (
        "exact_overlap_diagonal_distance",
        "--exact-overlap-diagonal-distance",
    ),
];

/// Arguments that set an alignment range.
const RANGE_ARGUMENTS: &[&str] = &[
    "reference_offset",
    "query_offset",
    "reference_limit",
    "query_limit",
    "rq_ranges",
    "use_embedded_rq_ranges",
];

/// Reject arguments that the chosen alignment method ignores, if they were given explicitly.
///
/// Arguments that keep their default value are never rejected.
fn reject_ignored_arguments(cli: &Cli, matches: &ArgMatches) -> Result<()> {
    let (ignored, method, other_method) = match cli.common.alignment_method {
        AlignmentMethod::AStarTemplateSwitch => (
            CHAIN_ONLY_ARGUMENTS,
            "a-star-template-switch",
            "a-star-chain-ts",
        ),
        AlignmentMethod::AStarChainTS => (
            A_STAR_ONLY_ARGUMENTS,
            "a-star-chain-ts",
            "a-star-template-switch",
        ),
    };
    for (id, flag) in ignored {
        if is_explicit(matches, id) {
            bail!(
                "{flag} is only used by --alignment-method {other_method}, \
                 but the alignment method is {method}; remove {flag}"
            );
        }
    }

    match cli.common.alignment_method {
        AlignmentMethod::AStarTemplateSwitch => {
            if cli.a_star.no_ts {
                for (id, flag) in TS_ONLY_ARGUMENTS {
                    if is_explicit(matches, id) {
                        bail!(
                            "{flag} has no effect with --no-ts, because template switches are disabled; remove {flag}"
                        );
                    }
                }
            }

            if cli.a_star.dont_extend_beyond_range
                && !RANGE_ARGUMENTS.iter().any(|id| is_explicit(matches, id))
            {
                bail!(
                    "--dont-extend-beyond-range has no effect without an alignment range \
                     (--rq-ranges, --use-embedded-rq-ranges, or --reference-offset etc.); remove it"
                );
            }
        }
        AlignmentMethod::AStarChainTS => reject_inexact_anchor_arguments_for_exact_anchors(
            matches,
            cli.chain.max_anchor_mutations,
            INEXACT_ANCHOR_ONLY_ARGUMENTS,
        )?,
    }

    Ok(())
}

pub fn cli(cli: Cli, matches: &ArgMatches) -> Result<()> {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    reject_ignored_arguments(&cli, matches)?;

    cli.common.alphabet.dispatch(Execute(cli))
}

struct Execute(Cli);

impl WithAlphabet for Execute {
    type Output = Result<()>;

    fn run<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(self) -> Result<()> {
        execute_with_alphabet::<AlphabetType>(self.0)
    }
}

fn execute_with_alphabet<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(
    cli: Cli,
) -> Result<()> {
    // Check the settings of the classic aligner before loading the inputs.
    let aligner = if cli.common.alignment_method == AlignmentMethod::AStarTemplateSwitch {
        let aligner = configure_aligner::<AlphabetType>(&cli);
        aligner.validate_settings().map_err(|error| {
            if cli.a_star.allow_ts_14_out_of_range == Ts14OutOfRangeStrategy::Allow {
                anyhow!("--allow-ts-14-out-of-range: {error}")
            } else {
                anyhow!(error)
            }
        })?;
        Some(aligner)
    } else {
        None
    };

    // Load input sequences.
    let (mut reference_record, mut query_record) =
        if let Some(CliPairInput { pair_fasta }) = &cli.input.pair_input {
            info!("Loading pair file {pair_fasta:?}");
            parse_pair_fasta_file(pair_fasta)?
        } else if let Some(CliSeparateInput { reference, query }) = &cli.input.separate_input {
            info!("Loading reference file {reference:?}");
            let reference = parse_single_fasta_file(reference)?;

            info!("Loading query file {query:?}");
            let query = parse_single_fasta_file(query)?;

            (reference, query)
        } else {
            return Err(anyhow!("No fasta input file given"));
        };

    // Save a copy of the original sequences for the computation of additional TSM starts and ends, if requested.
    let original_reference = reference_record.sequence_handle.clone();
    let original_query = query_record.sequence_handle.clone();

    // Remove skip characters.
    let skip_characters = cli.common.skip_characters.chars().collect::<Vec<_>>();
    ensure!(
        !cli.common.use_embedded_rq_ranges || !skip_characters.contains(&'|'),
        "Using embedded RQ ranges, but '|' is part of the skip characters"
    );
    reference_record
        .sequence_handle
        .retain(|c| !skip_characters.contains(&c));
    query_record
        .sequence_handle
        .retain(|c| !skip_characters.contains(&c));

    // Convert sequences to upper case.
    reference_record.sequence_handle.make_ascii_uppercase();
    query_record.sequence_handle.make_ascii_uppercase();

    // Parse RQ ranges.
    let range = if cli.common.use_embedded_rq_ranges {
        ensure!(
            cli.common.rq_ranges.is_none()
                && cli.common.reference_offset.is_none()
                && cli.common.reference_limit.is_none()
                && cli.common.query_offset.is_none()
                && cli.common.query_limit.is_none(),
            "Redundant specification of RQ ranges"
        );

        let reference_offset = reference_record.sequence_handle.find('|').ok_or_else(|| {
            anyhow!("Using embedded RQ ranges, but reference sequence contains no '|' character.")
        })?;
        let reference_limit = reference_offset + reference_record.sequence_handle[reference_offset+1..].find('|').ok_or_else(|| {
            anyhow!("Using embedded RQ ranges, but reference sequence contains only one '|' character.")
        })?;
        ensure!(
            reference_record.sequence_handle[reference_limit + 2..]
                .find('|')
                .is_none(),
            "Using embedded RQ ranges, but reference sequence contains more than two '|' characters"
        );
        reference_record.sequence_handle = reference_record.sequence_handle.replace('|', "");

        let query_offset = query_record.sequence_handle.find('|').ok_or_else(|| {
            anyhow!("Using embedded RQ ranges, but query sequence contains no '|' character.")
        })?;
        let query_limit = query_offset + query_record.sequence_handle[query_offset+1..].find('|').ok_or_else(|| {
            anyhow!("Using embedded RQ ranges, but query sequence contains only one '|' character.")
        })?;
        ensure!(
            query_record.sequence_handle[query_limit + 2..]
                .find('|')
                .is_none(),
            "Using embedded RQ ranges, but query sequence contains more than two '|' characters"
        );
        query_record.sequence_handle = query_record.sequence_handle.replace('|', "");

        AlignmentRange::new_offset_limit(
            AlignmentCoordinates::new(reference_offset, query_offset),
            AlignmentCoordinates::new(reference_limit, query_limit),
        )
    } else {
        parse_range(
            &cli,
            reference_record.sequence_handle.len(),
            query_record.sequence_handle.len(),
        )?
    };
    range.validate(
        reference_record.sequence_handle.len(),
        query_record.sequence_handle.len(),
    )?;

    // If requested, compute the additional starts and ends of TSMs outside of the given ranges.
    let additional_tsm_starts_and_ends =
        if cli.a_star.allow_ts_14_out_of_range == Ts14OutOfRangeStrategy::Allow {
            AdditionalExplicitTSMStartsAndEnds::new(
                &original_reference,
                &original_query,
                &range,
                &skip_characters,
                cli.common.use_embedded_rq_ranges,
            )?
        } else {
            AdditionalExplicitTSMStartsAndEnds::default()
        };

    // Convert sequences into the alphabet.
    let reference_sequence = sequence_from_ascii::<AlphabetType>(
        SequenceRole::Reference,
        reference_record.sequence_handle.as_bytes(),
    )?;
    let query_sequence = sequence_from_ascii::<AlphabetType>(
        SequenceRole::Query,
        query_record.sequence_handle.as_bytes(),
    )?;
    let reference_name = format!("{} {}", reference_record.id, reference_record.comment);
    let query_name = format!("{} {}", query_record.id, query_record.comment);

    debug!("Choosing alignment method...");
    match aligner {
        Some(aligner) => align_a_star_template_switch_distance(
            cli,
            aligner,
            reference_sequence.as_genome_subsequence(),
            query_sequence.as_genome_subsequence(),
            range,
            additional_tsm_starts_and_ends,
            &reference_name,
            &query_name,
        ),
        None => align_a_star_chain_ts(
            cli,
            reference_sequence.as_genome_subsequence(),
            query_sequence.as_genome_subsequence(),
            range,
            &reference_name,
            &query_name,
        ),
    }
}

fn parse_range(cli: &Cli, reference_length: usize, query_length: usize) -> Result<AlignmentRange> {
    let complete_reference_range = 0..reference_length;
    let complete_query_range = 0..query_length;

    let (reference_range, query_range) = if let Some(rq_ranges) = cli.common.rq_ranges.as_ref() {
        let (reference_range, query_range) = parse_rq_ranges(rq_ranges).map_err(|error| {
            anyhow!("Invalid --rq-ranges {rq_ranges:?}: {error}; expected e.g. R10..20Q12..25")
        })?;

        ensure!(
            reference_range.is_none()
                || (cli.common.reference_offset.is_none() && cli.common.reference_limit.is_none()),
            "--rq-ranges already specifies the reference range, so --reference-offset and --reference-limit must not be given"
        );
        ensure!(
            query_range.is_none()
                || (cli.common.query_offset.is_none() && cli.common.query_limit.is_none()),
            "--rq-ranges already specifies the query range, so --query-offset and --query-limit must not be given"
        );

        (
            reference_range.unwrap_or(complete_reference_range),
            query_range.unwrap_or(complete_query_range),
        )
    } else {
        (complete_reference_range, complete_query_range)
    };

    Ok(AlignmentRange::new_offset_limit(
        AlignmentCoordinates::new(
            cli.common.reference_offset.unwrap_or(reference_range.start),
            cli.common.query_offset.unwrap_or(query_range.start),
        ),
        AlignmentCoordinates::new(
            cli.common.reference_limit.unwrap_or(reference_range.end),
            cli.common.query_limit.unwrap_or(query_range.end),
        ),
    ))
}

/// The optional reference and query ranges given with `--rq-ranges`.
type RqRanges = (Option<Range<usize>>, Option<Range<usize>>);

/// Parse ranges of the form `R<offset>..<limit>Q<offset>..<limit>`, where each of the two parts is optional.
/// Whitespace is allowed after `R` and `Q`.
fn parse_rq_ranges(rq_ranges: &str) -> Result<RqRanges> {
    let mut rq_ranges = rq_ranges.chars().peekable();
    let mut reference_range = None;
    let mut query_range = None;

    let parse_number = |rq_ranges: &mut std::iter::Peekable<std::str::Chars>, what| {
        let mut number = String::new();
        while let Some(c) = rq_ranges.next_if(char::is_ascii_digit) {
            number.push(c);
        }
        number
            .parse::<usize>()
            .with_context(|| format!("expected a number as {what}"))
    };

    while let Some(rq) = rq_ranges.next() {
        while rq_ranges.next_if(|c| c.is_whitespace()).is_some() {}

        let offset = parse_number(&mut rq_ranges, "offset")?;
        ensure!(
            rq_ranges.next() == Some('.') && rq_ranges.next() == Some('.'),
            "expected '..' after the offset"
        );
        let limit = parse_number(&mut rq_ranges, "limit")?;

        let range = match rq {
            'R' => &mut reference_range,
            'Q' => &mut query_range,
            other => bail!("expected 'R' or 'Q', but found {other:?}"),
        };
        ensure!(range.is_none(), "the range for {rq} is given twice");
        *range = Some(offset..limit);
    }

    Ok((reference_range, query_range))
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::{
        A_STAR_ONLY_ARGUMENTS, CHAIN_ONLY_ARGUMENTS, Cli, INEXACT_ANCHOR_ONLY_ARGUMENTS,
        RANGE_ARGUMENTS, TS_ONLY_ARGUMENTS,
    };

    /// The argument tables name existing arguments, and the flags match the argument ids.
    #[test]
    fn argument_tables_name_existing_arguments() {
        let command = Cli::command();
        let find = |id: &str| {
            command
                .get_arguments()
                .find(|argument| argument.get_id() == id)
                .unwrap_or_else(|| panic!("no argument with id {id}"))
        };

        for (id, flag) in A_STAR_ONLY_ARGUMENTS
            .iter()
            .chain(CHAIN_ONLY_ARGUMENTS)
            .chain(TS_ONLY_ARGUMENTS)
            .chain(INEXACT_ANCHOR_ONLY_ARGUMENTS)
        {
            let argument = find(id);
            let actual_flag = match (argument.get_long(), argument.get_short()) {
                (Some(long), _) => format!("--{long}"),
                (None, Some(short)) => format!("-{short}"),
                (None, None) => panic!("argument {id} is positional"),
            };
            assert_eq!(&actual_flag, flag, "argument {id}");
        }
        for id in RANGE_ARGUMENTS {
            find(id);
        }

        // Every argument of a method-specific group is listed for its method.
        let mut method_specific_arguments = 0;
        for argument in command.get_arguments() {
            let id = argument.get_id().as_str();
            let heading = argument.get_help_heading().unwrap_or_default();
            if heading.starts_with("Classic aligner options") {
                method_specific_arguments += 1;
                assert!(
                    A_STAR_ONLY_ARGUMENTS.iter().any(|(other, _)| *other == id)
                        || id == "ts_node_ord_strategy",
                    "classic aligner argument {id} is not listed"
                );
            } else if heading.starts_with("tschainalign options") {
                method_specific_arguments += 1;
                assert!(
                    CHAIN_ONLY_ARGUMENTS.iter().any(|(other, _)| *other == id),
                    "tschainalign argument {id} is not listed"
                );
            }
        }
        assert_eq!(
            method_specific_arguments,
            A_STAR_ONLY_ARGUMENTS.len() + 1 + CHAIN_ONLY_ARGUMENTS.len()
        );
    }
}
