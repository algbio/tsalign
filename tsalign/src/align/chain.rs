use anyhow::{Result, ensure};
use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use generic_a_star::cost::AStarCost;
use lib_ts_chainalign::{
    chain_align::{
        limits::AlignmentLimits,
        performance_parameters::{
            AlignmentPerformanceParameters, AnchorPruning, InexactPruningParameters,
        },
    },
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::AlignmentCosts,
};
use lib_tsalign::alignment_geometry::AlignmentRange;
use log::{debug, info, warn};

use std::{fmt::Debug, fs::File, path::PathBuf};

use crate::{
    align::{Cli, InexactLowerBound},
    util::{
        infer_tschain_k, infer_tschain_max_n, load_tsa_config, tschain_preprocess_cache_file,
        validate_tschain_config, write_cache_file_atomically,
    },
};

pub fn align_a_star_chain_ts<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    reference_name: &str,
    query_name: &str,
) -> Result<()> {
    // tschainalign has several searches and data structures whose memory is not accounted for by a node count,
    // so an approximate memory limit like the one of the classic aligner would be misleading.
    ensure!(
        cli.a_star.memory_limit.is_none(),
        "--memory-limit is not supported by tschainalign (a-star-chain-ts); use --cost-limit or --iteration-limit"
    );
    let alignment_costs = load_tsa_config::<AlphabetType, _>(&cli.common.configuration_directory)?;
    let alignment_costs = AlignmentCosts::try_from(alignment_costs)?;

    let cache_directory = cli.chain.cache_directory.clone().unwrap_or_else(|| {
        warn!("No cache directory specified, dropping files into current working directory.");
        PathBuf::new()
    });
    let max_n = infer_tschain_max_n(reference.len(), query.len());
    let k = if let Some(k) = cli.chain.k {
        k
    } else {
        infer_tschain_k(reference.len(), query.len(), cli.chain.max_anchor_mutations).unwrap()
    };
    debug!("Using max_n = {max_n}");
    info!("Using k = {k}");
    let config = ChainingLowerBoundConfig {
        anchor_k: k,
        max_anchor_mutations: cli.chain.max_anchor_mutations,
        inexact_lower_bound: InexactLowerBoundKind::from(cli.chain.inexact_chaining_lower_bound),
    };
    validate_tschain_config(&config, cli.chain.k.is_some())?;

    let cache_file = tschain_preprocess_cache_file(
        &alignment_costs,
        &cache_directory,
        k,
        cli.chain.max_anchor_mutations,
        cli.chain.inexact_chaining_lower_bound.into(),
        max_n,
    )
    .unwrap();

    let chaining_lower_bounds = if let Ok(mut file) = File::open(&cache_file) {
        info!("Loading preprocessed data from cache at {cache_file:?}");
        let chaining_lower_bounds = ChainingLowerBounds::read(&mut file).unwrap();
        assert_eq!(chaining_lower_bounds.alignment_costs(), &alignment_costs);
        assert_eq!(chaining_lower_bounds.anchor_k(), k);
        assert_eq!(
            chaining_lower_bounds.max_anchor_mutations(),
            cli.chain.max_anchor_mutations
        );
        chaining_lower_bounds
    } else {
        assert!(
            !cli.chain.force_no_preprocessing,
            "Forcing no preprocessing, but no cache file found at {cache_file:?}."
        );

        info!("Preprocessing...");
        let chaining_lower_bounds = lib_ts_chainalign::preprocess(max_n, config, alignment_costs);

        info!("Storing preprocessed data into cache at {cache_file:?}");
        write_cache_file_atomically(&cache_file, |file| chaining_lower_bounds.write(file)).unwrap();
        chaining_lower_bounds
    };

    let anchor_pruning = if cli.chain.no_anchor_domination_pruning {
        warn!("--no-anchor-domination-pruning is deprecated, use --anchor-pruning none");
        AnchorPruning::None
    } else {
        cli.chain.anchor_pruning
    };
    // The inexact pruning options are rejected for the other modes, so that they are not silently ignored.
    ensure!(
        anchor_pruning == AnchorPruning::Inexact
            || (cli.chain.inexact_pruning_rules.is_none()
                && cli.chain.exact_overlap_diagonal_distance.is_none()),
        "--inexact-pruning-rules and --exact-overlap-diagonal-distance require --anchor-pruning inexact"
    );
    let default_inexact_pruning = InexactPruningParameters::default();
    let inexact_pruning = InexactPruningParameters {
        rules: cli
            .chain
            .inexact_pruning_rules
            .unwrap_or(default_inexact_pruning.rules),
        exact_overlap_diagonal_distance: cli
            .chain
            .exact_overlap_diagonal_distance
            .unwrap_or(default_inexact_pruning.exact_overlap_diagonal_distance),
    };

    let reference = reference.clone_as_vec();
    let query = query.clone_as_vec();
    let performance_parameters = AlignmentPerformanceParameters {
        max_successors: cli.chain.max_chaining_successors,
        max_exact_cost_function_cost: cli.chain.max_exact_cost_function_cost.into(),
        closed_list: cli.chain.chaining_closed_list,
        open_list: cli.chain.chaining_open_list,
        anchor_pruning,
        inexact_pruning,
    };

    let limits = AlignmentLimits {
        // Chain costs are at most `u32::MAX`, so a larger limit is no limit.
        cost_limit: cli.common.cost_limit.map(|cost_limit| {
            u32::try_from(cost_limit.as_u64())
                .unwrap_or(u32::MAX)
                .into()
        }),
        iteration_limit: cli.chain.iteration_limit,
    };

    let alignment = lib_ts_chainalign::align_with_limits(
        reference,
        query,
        range,
        &performance_parameters,
        &limits,
        &|c| {
            AlphabetType::character_to_ascii(
                AlphabetType::ascii_to_character(c).unwrap().complement(),
            )
        },
        reference_name,
        query_name,
        &chaining_lower_bounds,
    );
    info!("Finished aligning");

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{alignment}");
    Ok(())
}

impl From<InexactLowerBound> for InexactLowerBoundKind {
    fn from(value: InexactLowerBound) -> Self {
        match value {
            InexactLowerBound::WindowedMinMutations => InexactLowerBoundKind::WindowedMinMutations,
            InexactLowerBound::MaxMatchRun => InexactLowerBoundKind::MaxMatchRun,
        }
    }
}
