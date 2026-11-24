use std::{collections::BTreeSet, fmt::Debug, fs::create_dir_all, path::PathBuf};

use anyhow::{Context, Result, bail, ensure};
use clap::{ArgMatches, Parser};
use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::U32Cost;
use lib_ts_chainalign::{
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::AlignmentCosts,
    preprocess,
};
use log::{LevelFilter, info, warn};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

use crate::{
    align::InexactLowerBound,
    alphabet::{InputAlphabet, WithAlphabet},
    util::{
        infer_tschain_k, infer_tschain_max_n, load_tschain_config,
        reject_inexact_anchor_arguments_for_exact_anchors, tschain_preprocess_cache_file,
        validate_tschain_config, write_cache_file_atomically,
    },
};

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// The directory in which preprocessed data is stored.
    #[clap(long)]
    cache_directory: Option<PathBuf>,

    /// The alphabet present in the input files.
    ///
    /// This must also match the alphabet used in the config.
    #[clap(long, short = 'a', default_value = "dna")]
    alphabet: InputAlphabet,

    /// A directory containing the configuration files.
    ///
    /// See the README for its layout.
    /// Default: the built-in default config of tschainalign
    /// (equal to `test_files/config/chainalign` of the repository), independent of the working directory.
    #[clap(long, short = 'c')]
    configuration_directory: Option<PathBuf>,

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

    /// Maximum sequence length for which to preprocess.
    max_length: usize,
}

pub fn cli(cli: Cli, matches: &ArgMatches) -> Result<()> {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    reject_inexact_anchor_arguments_for_exact_anchors(
        matches,
        cli.max_anchor_mutations,
        &[(
            "inexact_chaining_lower_bound",
            "--inexact-chaining-lower-bound",
        )],
    )?;
    ensure!(cli.max_length > 0, "max_length must be at least 1");

    match cli.alphabet {
        #[cfg(feature = "alphabet_dna")]
        InputAlphabet::Dna => {}
        #[allow(
            unreachable_patterns,
            reason = "all alphabets other than DNA may be disabled"
        )]
        alphabet => {
            // Only A*-TS algo supports alphabets other than DNA.
            bail!("Unsupported alphabet type: {alphabet:?}");
        }
    }

    cli.alphabet.dispatch(Execute(cli))
}

struct Execute(Cli);

impl WithAlphabet for Execute {
    type Output = Result<()>;

    fn run<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(self) -> Result<()> {
        execute_with_alphabet::<AlphabetType>(self.0)
    }
}

fn execute_with_alphabet<AlphabetType: Alphabet>(cli: Cli) -> Result<()> {
    let costs =
        load_tschain_config::<AlphabetType, U32Cost>(cli.configuration_directory.as_deref())?;
    let costs = AlignmentCosts::try_from(costs)?;

    let cache_directory = cli.cache_directory.clone().unwrap_or_else(|| {
        warn!("No cache directory specified, dropping files into current working directory.");
        PathBuf::new()
    });
    create_dir_all(&cache_directory)
        .with_context(|| format!("Creating cache directory {cache_directory:?}"))?;

    let config = |k| ChainingLowerBoundConfig {
        anchor_k: k,
        max_anchor_mutations: cli.max_anchor_mutations,
        inexact_lower_bound: InexactLowerBoundKind::from(cli.inexact_chaining_lower_bound),
    };

    if let Some(k) = cli.k {
        info!("Using fixed k = {k} for all lengths");
        validate_tschain_config(&config(k), true)?;
        let max_n = infer_tschain_max_n(cli.max_length, cli.max_length);
        info!("Computing up to max_n = {max_n}");

        let mut current_max_n = max_n;
        // Also preprocess for max_n < k, since the aligner needs lower bounds for sequences shorter than k as well.
        while current_max_n >= 1 {
            info!("Preprocessing for max_n = {current_max_n}...");
            let lower_bounds = preprocess(current_max_n, config(k), costs.clone())?;
            let cache_file = tschain_preprocess_cache_file(
                &costs,
                &cache_directory,
                k,
                cli.max_anchor_mutations,
                cli.inexact_chaining_lower_bound.into(),
                current_max_n,
            )?;
            write_cache_file_atomically(&cache_file, |file| lower_bounds.write(file))?;

            current_max_n /= 2;
        }
    } else {
        info!("Using all possible inferred ks for each length");
        let max_n = infer_tschain_max_n(cli.max_length, cli.max_length);
        info!("Computing up to max_n = {max_n}");

        // Infer and validate all ks before preprocessing anything.
        let mut ks_per_max_n = Vec::new();
        let mut current_max_n = max_n;
        while current_max_n >= 1 {
            let next_max_n = current_max_n / 2;
            let mut ks = BTreeSet::new();
            for max_n in next_max_n + 1..=current_max_n {
                debug_assert_eq!(infer_tschain_max_n(max_n, max_n), current_max_n);
                ks.insert(infer_tschain_k(max_n, max_n, cli.max_anchor_mutations)?);
            }
            for &k in &ks {
                validate_tschain_config(&config(k), false)?;
            }
            ks_per_max_n.push((current_max_n, ks));

            current_max_n = next_max_n;
        }

        for (current_max_n, ks) in ks_per_max_n {
            for k in ks {
                info!("Preprocessing for max_n = {current_max_n} and k = {k}...");
                let lower_bounds = preprocess(current_max_n, config(k), costs.clone())?;
                let cache_file = tschain_preprocess_cache_file(
                    &costs,
                    &cache_directory,
                    k,
                    cli.max_anchor_mutations,
                    cli.inexact_chaining_lower_bound.into(),
                    current_max_n,
                )?;
                write_cache_file_atomically(&cache_file, |file| lower_bounds.write(file))?;
            }
        }
    }

    info!("Finished");
    Ok(())
}
