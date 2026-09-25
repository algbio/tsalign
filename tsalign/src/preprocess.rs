use std::{
    collections::BTreeSet,
    fs::{File, create_dir_all},
    path::PathBuf,
};

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};
use compact_genome::{
    implementation::alphabets::{
        dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
        dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet, rna_alphabet::RnaAlphabet,
        rna_alphabet_or_n::RnaAlphabetOrN,
        rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet,
    },
    interface::alphabet::Alphabet,
};
use generic_a_star::cost::U32Cost;
use lib_ts_chainalign::{costs::AlignmentCosts, preprocess};
use log::{LevelFilter, info, warn};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

use crate::util::{
    infer_tschain_k, infer_tschain_max_n, load_tsa_config, tschain_preprocess_cache_file,
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
    #[clap(long, short = 'c', default_value = "sample_tsa_config")]
    configuration_directory: PathBuf,

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

    /// Maximum sequence length for which to preprocess.
    max_length: usize,
}

#[derive(Debug, Clone, Eq, PartialEq, ValueEnum)]
enum InputAlphabet {
    Dna,
    DnaN,
    Rna,
    RnaN,
    DnaIupac,
    RnaIupac,
}

pub fn cli(cli: Cli) -> Result<()> {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    if cli.alphabet != InputAlphabet::Dna {
        // Only A*-TS algo supports alphabets other than DNA.
        bail!("Unsupported alphabet type: {:?}", cli.alphabet);
    }

    match cli.alphabet {
        InputAlphabet::Dna => execute_with_alphabet::<DnaAlphabet>(cli)?,
        InputAlphabet::DnaN => execute_with_alphabet::<DnaAlphabetOrN>(cli)?,
        InputAlphabet::Rna => execute_with_alphabet::<RnaAlphabet>(cli)?,
        InputAlphabet::RnaN => execute_with_alphabet::<RnaAlphabetOrN>(cli)?,
        InputAlphabet::DnaIupac => execute_with_alphabet::<DnaIupacNucleicAcidAlphabet>(cli)?,
        InputAlphabet::RnaIupac => execute_with_alphabet::<RnaIupacNucleicAcidAlphabet>(cli)?,
    }

    Ok(())
}

fn execute_with_alphabet<AlphabetType: Alphabet>(cli: Cli) -> Result<()> {
    let costs = load_tsa_config::<AlphabetType, U32Cost>(&cli.configuration_directory)?;
    let costs = AlignmentCosts::from(costs);

    let cache_directory = cli.cache_directory.clone().unwrap_or_else(|| {
        warn!("No cache directory specified, dropping files into current working directory.");
        PathBuf::new()
    });
    create_dir_all(&cache_directory)
        .with_context(|| format!("Creating cache directory {cache_directory:?}"))?;

    if let Some(k) = cli.k {
        info!("Using fixed k = {k} for all lengths");
        let max_n = infer_tschain_max_n(cli.max_length, cli.max_length);
        info!("Computing up to max_n = {max_n}");

        let mut current_max_n = max_n;
        while current_max_n >= 1 && current_max_n >= k.try_into()? {
            info!("Preprocessing for max_n = {current_max_n}...");
            let lower_bounds =
                preprocess(current_max_n, k, cli.max_anchor_mutations, costs.clone());
            let cache_file = tschain_preprocess_cache_file(
                &costs,
                &cache_directory,
                k,
                cli.max_anchor_mutations,
                current_max_n,
            )?;
            let mut file = File::create(&cache_file)
                .with_context(|| format!("Creating cache file {cache_file:?}"))?;
            lower_bounds.write(&mut file)?;

            current_max_n /= 2;
        }

        if current_max_n >= 1 {
            info!(
                "Not preprocessing for max_n = {current_max_n} or below, because it is smaller than k = {k}"
            );
        }
    } else {
        info!("Using all possible inferred ks for each length");
        let max_n = infer_tschain_max_n(cli.max_length, cli.max_length);
        info!("Computing up to max_n = {max_n}");

        let mut current_max_n = max_n;
        while current_max_n >= 1 {
            let next_max_n = current_max_n / 2;
            let mut ks = BTreeSet::new();
            for max_n in next_max_n + 1..=current_max_n {
                debug_assert_eq!(infer_tschain_max_n(max_n, max_n), current_max_n);
                ks.insert(infer_tschain_k(max_n, max_n, cli.max_anchor_mutations)?);
            }

            for k in ks {
                info!("Preprocessing for max_n = {current_max_n} and k = {k}...");
                let lower_bounds =
                    preprocess(current_max_n, k, cli.max_anchor_mutations, costs.clone());
                let cache_file = tschain_preprocess_cache_file(
                    &costs,
                    &cache_directory,
                    k,
                    cli.max_anchor_mutations,
                    current_max_n,
                )?;
                let mut file = File::create(&cache_file)
                    .with_context(|| format!("Creating cache file {cache_file:?}"))?;
                lower_bounds.write(&mut file)?;
            }

            current_max_n = next_max_n;
        }
    }

    info!("Finished");
    Ok(())
}
