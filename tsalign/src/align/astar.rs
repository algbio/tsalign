use std::fmt::Debug;

use anyhow::{Context, Result};
use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::{
    Aligner, aligner::configurable::PrimaryRangeStrategySelector,
    aligner::ts_14_out_of_range::AdditionalExplicitTSMStartsAndEnds,
    alignment_geometry::AlignmentRange,
};
use log::info;

use super::Cli;
use crate::util::load_tsa_config;

/// The node ordering of the classic aligner.
///
/// Only the anti-diagonal ordering is implemented, so the option has no effect.
/// It is kept as a hidden option so that existing command lines keep working.
#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchNodeOrdStrategySelector {
    AntiDiagonal,
}

/// Configure the classic aligner from the command line, except for its costs.
pub fn configure_aligner<AlphabetType: Alphabet>(cli: &Cli) -> Aligner<AlphabetType> {
    let mut aligner = Aligner::<AlphabetType>::default();
    aligner
        .set_min_length_strategy(cli.a_star.ts_min_length_strategy)
        .set_chaining_strategy(cli.a_star.ts_chaining_strategy)
        .set_total_length_strategy(cli.a_star.ts_total_length_strategy)
        .set_descendant_strategy(cli.a_star.ts_descendant_strategy)
        .set_no_ts(cli.a_star.no_ts)
        .set_ts_14_out_of_range_strategy(cli.a_star.allow_ts_14_out_of_range)
        .set_primary_range_strategy(PrimaryRangeStrategySelector::Prune)
        .set_force_label_correcting(cli.a_star.force_label_correcting);
    aligner
}

#[expect(clippy::too_many_arguments)]
pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    mut aligner: Aligner<AlphabetType>,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: AlignmentRange,
    additional_tsm_starts_and_ends: AdditionalExplicitTSMStartsAndEnds,
    reference_name: &str,
    query_name: &str,
) -> Result<()> {
    // Without `-c`, the aligner keeps its built-in default costs.
    if let Some(configuration_directory) = &cli.common.configuration_directory {
        aligner.set_costs(load_tsa_config(configuration_directory)?);
    } else {
        info!("Using the built-in default config");
    }

    let alignment = aligner.align_sequences(
        reference_name,
        reference,
        query_name,
        query,
        range,
        additional_tsm_starts_and_ends,
        cli.common.cost_limit,
        cli.a_star.memory_limit,
        !cli.a_star.dont_extend_beyond_range,
        cli.a_star.ts_uncertainty_range_extension_mode,
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
