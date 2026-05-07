use std::{
    fmt::Debug,
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use anyhow::{Result, anyhow, ensure};
use clap::{Args, Parser, ValueEnum};
use compact_genome::{
    implementation::{
        alphabets::{
            dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
            dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet,
            rna_alphabet::RnaAlphabet, rna_alphabet_or_n::RnaAlphabetOrN,
            rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet,
        },
        vec_sequence::VectorGenome,
        vec_sequence_store::VectorSequenceStore,
    },
    interface::{
        alphabet::Alphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
        sequence_store::SequenceStore,
    },
};
use lib_ts_chainalign::chain_align::performance_parameters::{
    ChainingClosedList, ChainingOpenList,
};
use lib_tsalign::{
    a_star_aligner::{
        alignment_geometry::{AlignmentCoordinates, AlignmentRange},
        gap_affine_edit_distance, gap_affine_edit_distance_a_star_align,
    },
    alignment_configuration::AlignmentConfiguration,
    alignment_matrix::AlignmentMatrix,
    costs::U64Cost,
};
use log::{LevelFilter, debug, info};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use template_switch_distance_type_selectors::{
    TemplateSwitchChainingStrategySelector, TemplateSwitchMinLengthStrategySelector,
    TemplateSwitchNodeOrdStrategySelector, TemplateSwitchTotalLengthStrategySelector,
    align_a_star_template_switch_distance,
};

use crate::align::{
    a_star_chain_ts::align_a_star_chain_ts,
    fasta_parser::{parse_pair_fasta_file, parse_single_fasta_file},
};

mod a_star_chain_ts;
mod fasta_parser;
mod template_switch_distance_type_selectors;

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    #[command(flatten)]
    input: CliInput,

    /// The file to store the alignment statistics in toml format.
    #[clap(long, short = 'o')]
    output: Option<PathBuf>,

    /// The directory in which preprocessed data is stored.
    #[clap(long)]
    cache_directory: Option<PathBuf>,

    /// If set, then the aligner performs no preprocessing, but instead crashes if the cached file is missing.
    #[clap(long)]
    force_no_preprocessing: bool,

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
    #[clap(long, short = 'c', default_value = "sample_tsa_config")]
    configuration_directory: PathBuf,

    #[clap(long, default_value = "a-star-template-switch")]
    alignment_method: AlignmentMethod,

    /// k-mer size for tschainalign.
    ///
    /// If it is not specified, it is inferred from the sequence lengths.
    #[clap(short)]
    k: Option<u32>,

    #[clap(long, default_value = "anti-diagonal")]
    ts_node_ord_strategy: TemplateSwitchNodeOrdStrategySelector,

    #[clap(long, default_value = "lookahead")]
    ts_min_length_strategy: TemplateSwitchMinLengthStrategySelector,

    #[clap(long, default_value = "none")]
    ts_chaining_strategy: TemplateSwitchChainingStrategySelector,

    /// If set to maximise, the total length of template switches becomes a secondary search criterion.
    ///
    /// As a result, instead of reporting an arbitrary alignment with optimal cost,
    /// tsalign selects an alignment with maximal total TS length out of all the alignments with optimal cost.
    #[clap(long, default_value = "maximise")]
    ts_total_length_strategy: TemplateSwitchTotalLengthStrategySelector,

    /// The maximum amount of successors to generate while processing a node during chaining.
    /// Can be tuned to optimise performance.
    ///
    /// This applies only to tschainalign.
    #[clap(long, default_value = "1")]
    max_chaining_successors: usize,

    /// The maximum cost until which the cost function is initialised exactly.
    ///
    /// Setting this to a higher value will increase the time required to initialise the cost function, but decrease the amount of chains computed.
    ///
    /// This applies only to tschainalign.
    #[clap(long, default_value = "1")]
    max_exact_cost_function_cost: u32,

    /// The closed list type to use for chaining.
    ///
    /// This applies only to tschainalign.
    #[clap(long, default_value = "special")]
    chaining_closed_list: ChainingClosedList,

    /// The open list type to use for chaining.
    ///
    /// This applies only to tschainalign.
    #[clap(long, default_value = "linear-heap")]
    chaining_open_list: ChainingOpenList,

    /// If set, template switches are not allowed.
    ///
    /// Use this to compare a template switch alignment against an alignment with out template switches.
    #[clap(long)]
    no_ts: bool,

    /// A cost limit for the alignment.
    ///
    /// If there is no alignment with at most that cost, the aligner will abort without result.
    #[clap(long)]
    cost_limit: Option<U64Cost>,

    /// An approximate memory limit in bytes for the aligner.
    ///
    /// If it is exceeded, then the aligner will abort without result.
    #[clap(long)]
    memory_limit: Option<usize>,

    /// Force the search to run in label-correcting mode.
    ///
    /// This is for debug purposes only.
    #[clap(long)]
    force_label_correcting: bool,

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

    /// If set, the aligner will not attempt to extend the alignment beyond the specified range
    /// without increasing cost.
    #[clap(long)]
    dont_extend_beyond_range: bool,
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
    Matrix,
    AStarGapAffine,
    AStarTemplateSwitch,
    AStarChainTS,
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

    if cli.alignment_method != AlignmentMethod::AStarTemplateSwitch
        && cli.alignment_method != AlignmentMethod::AStarChainTS
        && cli.alphabet != InputAlphabet::Dna
    {
        // Only A*-TS algo supports alphabets other than DNA.
        panic!("Unsupported alphabet type: {:?}", cli.alphabet);
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

fn execute_with_alphabet<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(
    cli: Cli,
) -> Result<()> {
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

    // Remove skip characters.
    let skip_characters = cli.skip_characters.chars().collect::<Vec<_>>();
    ensure!(
        !cli.use_embedded_rq_ranges || !skip_characters.contains(&'|'),
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
    let range = if cli.use_embedded_rq_ranges {
        ensure!(
            cli.rq_ranges.is_none()
                && cli.reference_offset.is_none()
                && cli.reference_limit.is_none()
                && cli.query_offset.is_none()
                && cli.query_limit.is_none(),
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
        )
    };

    // Move sequences into sequence store.
    let mut sequence_store = VectorSequenceStore::<AlphabetType>::new();
    let reference_record =
        reference_record.try_transform_handle::<_, anyhow::Error>(|sequence| {
            Ok(sequence_store.add(
                &VectorGenome::from_slice_u8(sequence.as_bytes()).map_err(|error| {
                    anyhow!("Reference contains non-alphabet character: {error}")
                })?,
            ))
        })?;
    let query_record = query_record.try_transform_handle::<_, anyhow::Error>(|sequence| {
        Ok(sequence_store.add(
            &VectorGenome::from_slice_u8(sequence.as_bytes())
                .map_err(|error| anyhow!("Query contains non-alphabet character: {error}"))?,
        ))
    })?;
    let reference_sequence = sequence_store.get(&reference_record.sequence_handle);
    let query_sequence = sequence_store.get(&query_record.sequence_handle);

    debug!("Choosing alignment method...");
    match cli.alignment_method {
        AlignmentMethod::Matrix => align_matrix(cli, reference_sequence, query_sequence),
        AlignmentMethod::AStarGapAffine => {
            align_a_star_gap_affine_edit_distance(cli, reference_sequence, query_sequence)
        }
        AlignmentMethod::AStarTemplateSwitch => align_a_star_template_switch_distance(
            cli,
            reference_sequence,
            query_sequence,
            range,
            &format!("{} {}", reference_record.id, reference_record.comment),
            &format!("{} {}", query_record.id, query_record.comment),
        ),
        AlignmentMethod::AStarChainTS => align_a_star_chain_ts(
            cli,
            reference_sequence,
            query_sequence,
            range,
            &format!("{} {}", reference_record.id, reference_record.comment),
            &format!("{} {}", query_record.id, query_record.comment),
        ),
    }

    Ok(())
}

fn align_matrix<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    if cli.output.is_some() {
        panic!("Outputting statistics not supported by matrix alignment");
    }

    #[derive(serde::Deserialize)]
    struct MatrixConfig {
        match_cost: u64,
        substitution_cost: u64,
        indel_cost: u64,
    }

    let mut config_path = cli.configuration_directory.clone();
    config_path.push("matrix.toml");
    let mut config_file = BufReader::new(File::open(config_path).unwrap());
    let mut config = String::new();
    config_file.read_to_string(&mut config).unwrap();
    let matrix_config: MatrixConfig = toml::from_str(&config).unwrap();

    let configuration = AlignmentConfiguration {
        match_cost: matrix_config.match_cost.into(),
        substitution_cost: matrix_config.substitution_cost.into(),
        insertion_cost: matrix_config.indel_cost.into(),
        deletion_cost: matrix_config.indel_cost.into(),
    };

    let mut alignment_matrix =
        AlignmentMatrix::<U64Cost>::new(configuration, reference.len(), query.len());
    let cost = alignment_matrix.align(reference, query);
    println!("Cost: {cost}");
}

fn align_a_star_gap_affine_edit_distance<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    #[derive(serde::Deserialize)]
    struct GapAffineConfig {
        match_cost: u64,
        substitution_cost: u64,
        gap_open_cost: u64,
        gap_extend_cost: u64,
    }

    let mut config_path = cli.configuration_directory.clone();
    config_path.push("a_star_gap_affine.toml");
    let mut config_file = BufReader::new(File::open(config_path).unwrap());
    let mut config = String::new();
    config_file.read_to_string(&mut config).unwrap();
    let gap_affine_config: GapAffineConfig = toml::from_str(&config).unwrap();

    let alignment = gap_affine_edit_distance_a_star_align(
        reference,
        query,
        gap_affine_edit_distance::ScoringTable::<U64Cost> {
            match_cost: gap_affine_config.match_cost.into(),
            substitution_cost: gap_affine_config.substitution_cost.into(),
            gap_open_cost: gap_affine_config.gap_open_cost.into(),
            gap_extend_cost: gap_affine_config.gap_extend_cost.into(),
        },
    );

    if let Some(output) = cli.output {
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{alignment}");
}

fn parse_range(cli: &Cli, reference_length: usize, query_length: usize) -> AlignmentRange {
    let complete_reference_range = 0..reference_length;
    let complete_query_range = 0..query_length;

    let (reference_range, query_range) = if let Some(rq_ranges) = cli.rq_ranges.as_ref() {
        let mut rq_ranges = rq_ranges.chars().peekable();

        let mut reference_range = None;
        let mut query_range = None;

        while rq_ranges.peek().is_some() {
            let rq = rq_ranges.next().unwrap();

            while let Some(c) = rq_ranges.peek() {
                if c.is_whitespace() {
                    rq_ranges.next().unwrap();
                } else {
                    break;
                }
            }

            let mut offset = String::new();
            while let Some(c) = rq_ranges.peek() {
                if c.is_numeric() {
                    offset.push(rq_ranges.next().unwrap());
                } else {
                    break;
                }
            }

            // Parse ..
            assert_eq!(rq_ranges.next(), Some('.'));
            assert_eq!(rq_ranges.next(), Some('.'));

            let mut limit = String::new();
            while let Some(c) = rq_ranges.peek() {
                if c.is_numeric() {
                    limit.push(rq_ranges.next().unwrap());
                } else {
                    break;
                }
            }

            let offset = offset.parse().unwrap();
            let limit = limit.parse().unwrap();

            match rq {
                'R' => {
                    assert!(reference_range.is_none());
                    reference_range = Some(offset..limit)
                }
                'Q' => {
                    assert!(query_range.is_none());
                    query_range = Some(offset..limit)
                }
                _ => panic!(),
            }
        }

        assert!(
            reference_range.is_none()
                || (cli.reference_offset.is_none() && cli.reference_limit.is_none())
        );
        assert!(query_range.is_none() || (cli.query_offset.is_none() && cli.query_limit.is_none()));

        (
            reference_range.unwrap_or(complete_reference_range),
            query_range.unwrap_or(complete_query_range),
        )
    } else {
        (complete_reference_range, complete_query_range)
    };

    AlignmentRange::new_offset_limit(
        AlignmentCoordinates::new(
            cli.reference_offset.unwrap_or(reference_range.start),
            cli.query_offset.unwrap_or(query_range.start),
        ),
        AlignmentCoordinates::new(
            cli.reference_limit.unwrap_or(reference_range.end),
            cli.query_limit.unwrap_or(query_range.end),
        ),
    )
}
