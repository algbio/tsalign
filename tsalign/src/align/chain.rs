use anyhow::{Context, Result, bail, ensure};
use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use lib_ts_chainalign::{
    chain_align::performance_parameters::AlignmentPerformanceParameters,
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
    let reference = reference.clone_as_vec();
    let query = query.clone_as_vec();
    // Check the inputs before preprocessing.
    lib_ts_chainalign::validate_sequences(&reference, &query, &range)?;
    ensure!(
        !reference.is_empty() || !query.is_empty(),
        "tschainalign needs at least one non-empty sequence"
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
        infer_tschain_k(reference.len(), query.len(), cli.chain.max_anchor_mutations)?
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
    )?;

    let chaining_lower_bounds = if let Ok(mut file) = File::open(&cache_file) {
        info!("Loading preprocessed data from cache at {cache_file:?}");
        let chaining_lower_bounds = ChainingLowerBounds::read(&mut file)
            .with_context(|| format!("Reading cache file {cache_file:?}"))?;
        if chaining_lower_bounds.alignment_costs() != &alignment_costs
            || chaining_lower_bounds.anchor_k() != k
            || chaining_lower_bounds.max_anchor_mutations() != cli.chain.max_anchor_mutations
        {
            bail!(
                "The cache file {cache_file:?} was preprocessed with different costs, k or max_anchor_mutations; \
                 remove it to preprocess again"
            );
        }
        chaining_lower_bounds
    } else {
        if cli.chain.force_no_preprocessing {
            bail!("--force-no-preprocessing is set, but no cache file found at {cache_file:?}");
        }

        info!("Preprocessing...");
        let chaining_lower_bounds = lib_ts_chainalign::preprocess(max_n, config, alignment_costs)?;

        info!("Storing preprocessed data into cache at {cache_file:?}");
        write_cache_file_atomically(&cache_file, |file| chaining_lower_bounds.write(file))?;
        chaining_lower_bounds
    };

    let performance_parameters = AlignmentPerformanceParameters {
        max_successors: cli.chain.max_chaining_successors,
        max_exact_cost_function_cost: cli.chain.max_exact_cost_function_cost.into(),
        closed_list: cli.chain.chaining_closed_list,
        open_list: cli.chain.chaining_open_list,
        prune_dominated_anchors: !cli.chain.no_anchor_domination_pruning,
    };

    let alignment = lib_ts_chainalign::align(
        reference,
        query,
        range,
        &performance_parameters,
        &|c| {
            AlphabetType::character_to_ascii(
                AlphabetType::ascii_to_character(c).unwrap().complement(),
            )
        },
        reference_name,
        query_name,
        &chaining_lower_bounds,
    )?;
    info!("Finished aligning");

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(
            std::fs::File::create(&output)
                .with_context(|| format!("Creating output file {output:?}"))?,
        );
        write!(output, "{}", toml::to_string(&alignment)?)?;
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
