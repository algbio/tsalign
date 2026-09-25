use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use lib_ts_chainalign::{
    chain_align::performance_parameters::AlignmentPerformanceParameters,
    chaining_lower_bounds::ChainingLowerBounds, costs::AlignmentCosts,
};
use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;
use log::{debug, info, warn};

use std::{fmt::Debug, fs::File, path::PathBuf};

use crate::{
    align::Cli,
    util::{infer_tschain_k, infer_tschain_max_n, load_tsa_config, tschain_preprocess_cache_file},
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
) {
    let alignment_costs = load_tsa_config::<AlphabetType, _>(&cli.configuration_directory).unwrap();
    let alignment_costs: AlignmentCosts<_> = alignment_costs.into();

    let cache_directory = cli.cache_directory.clone().unwrap_or_else(|| {
        warn!("No cache directory specified, dropping files into current working directory.");
        PathBuf::new()
    });
    let max_n = infer_tschain_max_n(reference.len(), query.len());
    let k = if let Some(k) = cli.k {
        k
    } else {
        infer_tschain_k(reference.len(), query.len(), cli.max_anchor_mutations).unwrap()
    };
    debug!("Using max_n = {max_n}");
    info!("Using k = {k}");
    let cache_file =
        tschain_preprocess_cache_file(&alignment_costs, &cache_directory, k, max_n).unwrap();

    let chaining_lower_bounds = if let Ok(mut file) = File::open(&cache_file) {
        info!("Loading preprocessed data from cache at {cache_file:?}");
        let chaining_lower_bounds = ChainingLowerBounds::read(&mut file).unwrap();
        assert_eq!(chaining_lower_bounds.alignment_costs(), &alignment_costs);
        assert_eq!(chaining_lower_bounds.anchor_k(), k);
        chaining_lower_bounds
    } else {
        assert!(
            !cli.force_no_preprocessing,
            "Forcing no preprocessing, but no cache file found."
        );

        info!("Preprocessing...");
        let chaining_lower_bounds =
            lib_ts_chainalign::preprocess(max_n, k, cli.max_anchor_mutations, alignment_costs);

        info!("Storing preprocessed data into cache at {cache_file:?}");
        let mut file = File::create(&cache_file).unwrap();
        chaining_lower_bounds.write(&mut file).unwrap();
        chaining_lower_bounds
    };

    let reference = reference.clone_as_vec();
    let query = query.clone_as_vec();
    let performance_parameters = AlignmentPerformanceParameters {
        max_successors: cli.max_chaining_successors,
        max_exact_cost_function_cost: cli.max_exact_cost_function_cost.into(),
        closed_list: cli.chaining_closed_list,
        open_list: cli.chaining_open_list,
    };

    let alignment = lib_ts_chainalign::align::<AlphabetType>(
        reference,
        query,
        range,
        cli.max_anchor_mutations,
        &performance_parameters,
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
}
