use std::{
    fs::File,
    io::{Read, Write, stdout},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::Parser;
use lib_tsalign::{alignment_result::AlignmentResult, costs::U64Cost};
use lib_tsshow::{
    plain_text::show_template_switches,
    svg::{SvgConfig, UncertaintyRangeMode, create_error_svg, create_ts_svg},
    svg_to_png,
};
use log::{LevelFilter, error, info, warn};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// Path to a toml output file of `tsalign align`.
    ///
    /// This is the input from which template switches are visualised.
    #[clap(long, short = 'i')]
    input: PathBuf,

    /// Path to a toml output file of `tsalign align` created with the `--no-ts` flag.
    ///
    /// This input can be used to provide an alternative alignment without template switches.
    /// It is expected to contain no template switches, otherwise it will be rejected.
    #[clap(long, short = 'n')]
    no_ts_input: Option<PathBuf>,

    /// Create an SVG image file showing the complete alignment with all template switches.
    #[clap(long, short = 's')]
    svg: Option<PathBuf>,

    /// If set together with --svg, then the SVG image is rendered to a PNG image as well.
    #[clap(long, short = 'p')]
    png: bool,

    /// Output the template switches in plain text format to stdout.
    #[clap(long, short = 't')]
    plain_text: bool,

    /// Always render the SVG (and optionally PNG), even if an error occurs.
    /// In the case of an error, the SVG simply contains the error message.
    #[clap(long, short = 'r')]
    render_always: bool,

    /// Draw the SVG image with arrows connecting the switchpoints.
    /// Ignored if --svg is not set.
    #[clap(long, short = 'a')]
    svg_arrows: bool,

    /// Draw the SVG image with more complement characters than just the bare minimum needed to visualise the template switch.
    #[clap(long, short = 'c')]
    more_svg_complement: bool,

    /// Draw only the given number of characters around the TSMs.
    /// Ignored if --svg is not set.
    #[clap(long, short = 'z')]
    restrict_context_to: Option<usize>,

    /// Add a visualisation of the TSM uncertainty ranges to the SVG output.
    #[clap(long, short = 'e', default_value = "none", default_missing_value = "full", num_args = 0..=1)]
    uncertainty_range_mode: UncertaintyRangeMode,
}

pub fn cli(cli: Cli) -> Result<()> {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        // Log to stderr, such that stdout contains only the plain text output.
        TerminalMode::Stderr,
        ColorChoice::Auto,
    )
    .unwrap();

    if cli.svg.is_none() && !cli.plain_text {
        error!("Neither --svg nor --plain-text is set. Nothing to do.");
        bail!("Neither --svg nor --plain-text is set. Nothing to do.");
    }

    info!("Reading tsalign output toml file {:?}", cli.input);
    let result = read_alignment(&cli.input)?;
    let no_ts_result = cli
        .no_ts_input
        .as_ref()
        .map(|no_ts_input| {
            info!("Reading tsalign no-ts output toml file {no_ts_input:?}");
            read_alignment(no_ts_input)
        })
        .transpose()?;

    if cli.plain_text {
        show_template_switches(stdout(), &result, &no_ts_result)?;
    }

    if let Some(svg_out_path) = cli.svg.as_ref() {
        let mut svg = Vec::new();
        if let Err(error) = create_ts_svg(
            &mut svg,
            &result,
            &no_ts_result,
            &SvgConfig {
                render_arrows: cli.svg_arrows,
                render_more_complement: cli.more_svg_complement,
                restrict_context: cli.restrict_context_to,
                uncertainty_range_mode: cli.uncertainty_range_mode,
            },
        ) {
            if cli.render_always {
                create_error_svg(&mut svg, error).with_context(|| "Error creating error SVG")?;
            } else {
                return Err(error).with_context(|| "Error creating SVG.");
            }
        }
        let svg = svg;

        info!("Writing svg to {svg_out_path:?}");
        write_file(svg_out_path, &svg)?;

        if cli.png {
            let png = svg_to_png(&svg, 20.0);

            let png_out_path = svg_out_path.with_extension("png");
            if &png_out_path == svg_out_path {
                warn!(
                    "SVG and PNG filenames are the same ({png_out_path:?}). Do not use '.png' as the extension for your SVG file. Overwriting SVG file..."
                );
            }

            info!("Writing png to {png_out_path:?}");
            write_file(&png_out_path, &png)?;
        }
    }

    Ok(())
}

/// Reads a toml output file of `tsalign align`.
fn read_alignment(path: &Path) -> Result<AlignmentResult<U64Cost>> {
    let mut buffer = String::new();
    File::open(path)
        .with_context(|| format!("Error opening input file {path:?}"))?
        .read_to_string(&mut buffer)
        .with_context(|| format!("Error reading input file {path:?}"))?;
    toml::from_str(&buffer).with_context(|| format!("Error parsing input file {path:?}"))
}

fn write_file(path: &Path, content: &[u8]) -> Result<()> {
    File::create(path)
        .with_context(|| format!("Error creating output file {path:?}"))?
        .write_all(content)
        .with_context(|| format!("Error writing output file {path:?}"))
}
