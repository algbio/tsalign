use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use clap::{ArgMatches, parser::ValueSource};
use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::AStarError;
use generic_a_star::cost::{AStarCost, U32Cost};
use lib_ts_chainalign::{
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, DEFAULT_CONFIG_TSA},
};
use lib_tsalign::{arena_index::ArenaIndexType, config::TemplateSwitchConfig};
use log::info;

pub fn load_tsa_config<AlphabetType: Alphabet, Cost: AStarCost>(
    configuration_directory: impl AsRef<Path>,
) -> Result<TemplateSwitchConfig<AlphabetType, Cost>> {
    let mut config_path = configuration_directory.as_ref().to_path_buf();
    info!("Loading alignment config directory {config_path:?}");

    config_path.push("config.tsa");
    let config_file = std::io::BufReader::new(
        std::fs::File::open(&config_path)
            .with_context(|| format!("Error opening config file {config_path:?}"))?,
    );
    let costs = TemplateSwitchConfig::read_plain(config_file)
        .with_context(|| format!("Error parsing template switch config from {config_path:?}"))?;
    Ok(costs)
}

/// Load the config of tschainalign from the configuration directory,
/// or the built-in default config of tschainalign if no directory is given.
pub fn load_tschain_config<AlphabetType: Alphabet, Cost: AStarCost>(
    configuration_directory: Option<&Path>,
) -> Result<TemplateSwitchConfig<AlphabetType, Cost>> {
    if let Some(configuration_directory) = configuration_directory {
        load_tsa_config(configuration_directory)
    } else {
        info!("Using the built-in default config of tschainalign");
        TemplateSwitchConfig::read_plain(DEFAULT_CONFIG_TSA.as_bytes())
            .context("Error parsing the built-in default config of tschainalign")
    }
}

/// Returns true if the argument with the given id was given on the command line,
/// as opposed to being set by its default value or not at all.
pub fn is_explicit(matches: &ArgMatches, id: &str) -> bool {
    matches.value_source(id) == Some(ValueSource::CommandLine)
}

/// Reject arguments that only affect inexact anchors, if they are given explicitly while the anchors are exact.
///
/// `arguments` are pairs of argument id and the argument as written on the command line.
pub fn reject_inexact_anchor_arguments_for_exact_anchors(
    matches: &ArgMatches,
    max_anchor_mutations: u8,
    arguments: &[(&str, &str)],
) -> Result<()> {
    if max_anchor_mutations != 0 {
        return Ok(());
    }
    for (id, flag) in arguments {
        if is_explicit(matches, id) {
            bail!(
                "{flag} has no effect with exact anchors (--max-anchor-mutations 0); \
                 remove {flag} or allow mutations in anchors"
            );
        }
    }
    Ok(())
}

pub use lib_ts_chainalign::cache::infer_max_n as infer_tschain_max_n;

pub fn infer_tschain_k(
    reference_length: usize,
    query_length: usize,
    max_anchor_mutations: u8,
) -> Result<u32> {
    Ok(lib_ts_chainalign::cache::infer_k(
        reference_length,
        query_length,
        max_anchor_mutations,
    ))
}

/// Check that tschainalign supports the given `k` and `--max-anchor-mutations`,
/// before looking up cached lower bounds or computing them.
///
/// `k_given` tells if `k` was given with `-k`, or inferred from the sequence lengths.
pub fn validate_tschain_config(config: &ChainingLowerBoundConfig, k_given: bool) -> Result<()> {
    config.validate().map_err(|error| {
        let k_source = if k_given {
            "given with -k"
        } else {
            "inferred from the sequence lengths and --max-anchor-mutations; pass -k to choose it explicitly"
        };
        anyhow!("Unsupported tschainalign parameters: {error} (k was {k_source})")
    })
}

/// Write a cache file such that concurrent readers never observe a partially written file.
///
/// See [`lib_ts_chainalign::cache::write_cache_file_atomically`].
pub fn write_cache_file_atomically(
    cache_file: impl AsRef<Path>,
    write: impl FnOnce(&mut std::io::BufWriter<&std::fs::File>) -> std::io::Result<()>,
) -> Result<()> {
    Ok(lib_ts_chainalign::cache::write_cache_file_atomically(
        cache_file, write,
    )?)
}

pub fn tschain_preprocess_cache_file(
    alignment_costs: &AlignmentCosts<U32Cost>,
    cache_directory: impl AsRef<Path>,
    k: u32,
    max_anchor_mutations: u8,
    inexact_lower_bound: InexactLowerBoundKind,
    max_n: usize,
) -> Result<PathBuf> {
    Ok(lib_ts_chainalign::cache::cache_file_path(
        alignment_costs,
        cache_directory,
        k,
        max_anchor_mutations,
        &inexact_lower_bound,
        max_n,
    ))
}

/// Converts the error of a failed A* search into the error message of the CLI.
///
/// If the node arena is full and a larger `--arena-index` exists, then the message suggests it.
pub fn search_error(error: &AStarError, arena_index: ArenaIndexType) -> anyhow::Error {
    match (error, arena_index) {
        (AStarError::ArenaFull(_), ArenaIndexType::U32) => anyhow!(
            "{error}; rerun with `--arena-index u64` to allow more nodes (this needs a bit more memory per node)"
        ),
        _ => anyhow!("{error}"),
    }
}

#[cfg(test)]
mod tests {
    use generic_a_star::{AStarError, arena::ArenaFullError};
    use lib_tsalign::arena_index::ArenaIndexType;

    use super::search_error;

    #[test]
    fn full_arena_suggests_u64() {
        let error = AStarError::ArenaFull(ArenaFullError {
            capacity: 4_294_967_295,
            index_type: "u32",
        });
        assert_eq!(
            search_error(&error, ArenaIndexType::U32).to_string(),
            "the A* node arena is full: it can hold at most 4294967295 nodes with u32 node indices; \
             rerun with `--arena-index u64` to allow more nodes (this needs a bit more memory per node)"
        );

        let error = AStarError::ArenaFull(ArenaFullError {
            capacity: usize::MAX,
            index_type: "u64",
        });
        assert_eq!(
            search_error(&error, ArenaIndexType::U64).to_string(),
            format!(
                "the A* node arena is full: it can hold at most {} nodes with u64 node indices",
                usize::MAX
            )
        );
    }
}
