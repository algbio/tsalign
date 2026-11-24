//! End-to-end tests of `tsalign show` with all its options and error cases (issue #279).
//!
//! The inputs are the alignment fixtures of the lib_tsshow tests (`lib_tsshow/tests/fixtures`)
//! and alignments in `test_files`. The outputs are written to `target/tsalign-show-e2e/<test>/`.

use std::{
    fs,
    path::{Path, PathBuf},
    process::Output,
};

use anyhow::Result;
use util::run_in_repo_root_output;

mod util;

const FIXTURES: &str = "lib_tsshow/tests/fixtures";
const TS: &str = "lib_tsshow/tests/fixtures/classic_both_descendants.toml";
const NO_TS: &str = "lib_tsshow/tests/fixtures/classic_both_descendants.no_ts.toml";

/// The output directory of a test, relative to the repository root (where tsalign runs).
fn output_directory(test: &str) -> PathBuf {
    let relative = PathBuf::from("target").join("tsalign-show-e2e").join(test);
    let absolute = repository_root().join(&relative);
    if absolute.exists() {
        fs::remove_dir_all(&absolute).unwrap();
    }
    fs::create_dir_all(&absolute).unwrap();
    relative
}

fn repository_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_owned()
}

fn show(args: &str) -> Result<Output> {
    run_in_repo_root_output(&format!("show {args}"))
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn assert_success(args: &str) -> Result<Output> {
    let output = show(args)?;
    assert!(
        output.status.success(),
        "`tsalign show {args}` failed\nstderr:\n{}",
        stderr(&output)
    );
    Ok(output)
}

/// Asserts that tsalign rejects the arguments with exit code 1 and an error message containing `expected`,
/// without panicking.
fn assert_error(args: &str, expected: &str) -> Result<()> {
    let output = show(args)?;
    let stderr = stderr(&output);
    assert!(
        !stderr.contains("panicked"),
        "`tsalign show {args}` panicked instead of reporting an error\nstderr:\n{stderr}"
    );
    assert_eq!(
        output.status.code(),
        Some(1),
        "`tsalign show {args}`: unexpected exit status\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains(expected),
        "`tsalign show {args}`: the error does not contain {expected:?}\nstderr:\n{stderr}"
    );
    Ok(())
}

fn read(relative_path: &Path) -> Vec<u8> {
    fs::read(repository_root().join(relative_path))
        .unwrap_or_else(|error| panic!("reading {relative_path:?}: {error}"))
}

fn assert_svg(relative_path: &Path) -> String {
    let svg = String::from_utf8(read(relative_path)).unwrap();
    assert!(svg.starts_with("<svg"), "{relative_path:?} is not an SVG");
    assert!(
        svg.trim_end().ends_with("</svg>"),
        "{relative_path:?} is truncated"
    );
    assert!(
        svg.contains("@font-face"),
        "{relative_path:?} embeds no fonts"
    );
    svg
}

#[test]
fn plain_text() -> Result<()> {
    let output = assert_success(&format!("-i {TS} -t"))?;
    let stdout = stdout(&output);
    assert_eq!(stdout.matches("Switch process:").count(), 2, "{stdout}");
    assert!(!stdout.contains("No-ts alignment:"), "{stdout}");
    Ok(())
}

/// The log messages go to stderr, so that stdout contains only the plain text output.
#[test]
fn plain_text_contains_no_log_messages() -> Result<()> {
    let output = assert_success(&format!("-i {TS} -t --log-level trace"))?;
    let stdout = stdout(&output);
    assert!(stdout.starts_with("Child: "), "{stdout}");
    for level in ["INFO", "DEBUG", "TRACE"] {
        assert!(!stdout.contains(&format!("[{level}]")), "{stdout}");
    }
    Ok(())
}

#[test]
fn plain_text_with_no_ts_alignment() -> Result<()> {
    let output = assert_success(&format!("-i {TS} -n {NO_TS} --plain-text"))?;
    let stdout = stdout(&output);
    assert_eq!(stdout.matches("Switch process:").count(), 2, "{stdout}");
    assert_eq!(stdout.matches("No-ts alignment:").count(), 2, "{stdout}");
    Ok(())
}

#[test]
fn svg() -> Result<()> {
    let directory = output_directory("svg");
    let svg_path = directory.join("out.svg");
    let output = assert_success(&format!("-i {TS} -s {}", svg_path.display()))?;
    // Without -t, nothing is written to stdout.
    assert!(!stdout(&output).contains("Switch process:"));
    let svg = assert_svg(&svg_path);
    assert!(svg.contains("TS-A inner:") && svg.contains("TS-B inner:"));
    assert!(!svg.contains("arrow_head_red)"), "arrows without -a");
    assert!(
        !repository_root().join(directory.join("out.png")).exists(),
        "PNG without -p"
    );
    Ok(())
}

#[test]
fn svg_with_all_options() -> Result<()> {
    let directory = output_directory("svg_with_all_options");
    let svg_path = directory.join("out.svg");
    assert_success(&format!(
        "-i {TS} -n {NO_TS} -s {} -a -c -z 3 -e full -t",
        svg_path.display()
    ))?;
    let svg = assert_svg(&svg_path);
    assert_eq!(svg.matches("url(#arrow_head_red)").count(), 4);
    assert!(svg.contains("Uncertainty range of the TSM"));
    Ok(())
}

/// `-e` without a value means `full`, and `inner-only` and `none` are accepted.
#[test]
fn uncertainty_range_modes() -> Result<()> {
    let directory = output_directory("uncertainty_range_modes");
    // A fixture with non-empty uncertainty ranges.
    let input = format!("{FIXTURES}/classic_large_uncertainty.toml");
    let mut svgs = Vec::new();
    for (name, option) in [
        ("default", ""),
        ("none", "-e none"),
        ("inner_only", "-e inner-only"),
        ("full", "-e full"),
        ("missing_value", "-e"),
    ] {
        let svg_path = directory.join(format!("{name}.svg"));
        // `-e` without a value must be the last argument.
        assert_success(&format!("-i {input} -s {} {option}", svg_path.display()))?;
        svgs.push(assert_svg(&svg_path));
    }
    assert_eq!(svgs[0], svgs[1], "the default is none");
    assert_eq!(svgs[3], svgs[4], "-e without a value is full");
    assert_ne!(svgs[1], svgs[3]);
    assert_ne!(svgs[2], svgs[3]);
    Ok(())
}

#[test]
fn png() -> Result<()> {
    let directory = output_directory("png");
    let svg_path = directory.join("out.svg");
    assert_success(&format!(
        "-i {TS} -n {NO_TS} -s {} -p -a",
        svg_path.display()
    ))?;
    assert_svg(&svg_path);
    let png = read(&directory.join("out.png"));
    assert!(
        png.starts_with(b"\x89PNG\r\n\x1a\n"),
        "the PNG file has no PNG header"
    );
    assert!(png.len() > 1000, "the PNG file is too small");
    Ok(())
}

/// `-p` without `-s` does nothing.
#[test]
fn png_without_svg_is_ignored() -> Result<()> {
    assert_success(&format!("-i {TS} -t -p"))?;
    Ok(())
}

/// The output does not depend on the process (e.g. hash map ordering or randomised hashing).
#[test]
fn output_is_deterministic() -> Result<()> {
    let directory = output_directory("output_is_deterministic");
    let mut outputs = Vec::new();
    for run in 0..2 {
        let svg_path = directory.join(format!("{run}.svg"));
        let output = assert_success(&format!(
            "-i {TS} -n {NO_TS} -s {} -a -c -e full -t",
            svg_path.display()
        ))?;
        outputs.push((read(&svg_path), output.stdout));
    }
    assert!(outputs[0] == outputs[1], "two runs gave different outputs");
    Ok(())
}

/// Every fixture can be shown with every output.
#[test]
fn all_fixtures() -> Result<()> {
    let directory = output_directory("all_fixtures");
    let mut names: Vec<String> = fs::read_dir(repository_root().join(FIXTURES))?
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter_map(|name| name.strip_suffix(".toml").map(str::to_owned))
        .filter(|name| !name.ends_with(".no_ts"))
        .collect();
    names.sort();
    assert!(names.len() >= 20);

    for name in names {
        let svg_path = directory.join(format!("{name}.svg"));
        assert_success(&format!(
            "-i {FIXTURES}/{name}.toml -n {FIXTURES}/{name}.no_ts.toml -s {} -t -a -e -l warn",
            svg_path.display()
        ))?;
        assert_svg(&svg_path);
    }
    Ok(())
}

/// The commands of the old show-bug scripts (`test_files/ari_200_show_bug.sh`, `test_files/ari_868_show_bug.sh`,
/// PRs #219 and #220).
#[test]
fn show_bug_scripts() -> Result<()> {
    let directory = output_directory("show_bug_scripts");
    for (name, no_ts) in [("ari_200", "ari_200_nots"), ("ari_868", "ari_868_nots")] {
        let svg_path = directory.join(format!("{name}.svg"));
        assert_success(&format!(
            "-i test_files/{name}.toml -n test_files/{no_ts}.toml -s {} -pare -z 5 --log-level warn",
            svg_path.display()
        ))?;
        assert_svg(&svg_path);
        assert!(read(&directory.join(format!("{name}.png"))).starts_with(b"\x89PNG"));
    }
    Ok(())
}

#[test]
fn nothing_to_do() -> Result<()> {
    assert_error(&format!("-i {TS}"), "Neither --svg nor --plain-text is set")
}

#[test]
fn missing_input_file() -> Result<()> {
    assert_error(
        "-i test_files/does_not_exist.toml -t",
        "does_not_exist.toml",
    )
}

#[test]
fn missing_no_ts_input_file() -> Result<()> {
    assert_error(
        &format!("-i {TS} -n test_files/does_not_exist.toml -t"),
        "does_not_exist.toml",
    )
}

/// A file that is not a toml output of `tsalign align`.
#[test]
fn invalid_input_file() -> Result<()> {
    assert_error("-i test_files/twin_a.fa -t", "twin_a.fa")
}

#[test]
fn missing_input_argument() -> Result<()> {
    let output = show("-t")?;
    // Clap exits with code 2 on invalid arguments.
    assert_eq!(output.status.code(), Some(2), "{}", stderr(&output));
    assert!(stderr(&output).contains("--input"));
    Ok(())
}

/// A "no-TS" alignment that contains template switches is rejected.
#[test]
fn no_ts_alignment_with_template_switches() -> Result<()> {
    let directory = output_directory("no_ts_alignment_with_template_switches");
    assert_error(&format!("-i {TS} -n {TS} -t"), "template switches")?;
    assert_error(
        &format!("-i {TS} -n {TS} -s {}", directory.join("out.svg").display()),
        "template switches",
    )
}

/// A no-TS alignment of other sequences is rejected.
#[test]
fn mismatched_sequences() -> Result<()> {
    let directory = output_directory("mismatched_sequences");
    let other = format!("{FIXTURES}/classic_qr_reverse.no_ts.toml");
    assert_error(&format!("-i {TS} -n {other} -t"), "other sequences")?;
    assert_error(
        &format!(
            "-i {TS} -n {other} -s {}",
            directory.join("out.svg").display()
        ),
        "other sequences",
    )
}

/// An alignment without template switches has nothing to show, but is not an error.
#[test]
fn alignment_without_template_switches() -> Result<()> {
    let directory = output_directory("alignment_without_template_switches");
    let svg_path = directory.join("out.svg");
    let output = assert_success(&format!("-i {NO_TS} -t -s {}", svg_path.display()))?;
    assert!(!stdout(&output).contains("Switch process:"));
    let svg = assert_svg(&svg_path);
    assert!(!svg.contains("inner:"));
    Ok(())
}

/// An alignment without target (the aligner stopped at its cost limit) cannot be rendered as SVG,
/// unless `--render-always` (`-r`) renders the error into the SVG.
#[test]
fn alignment_without_target() -> Result<()> {
    let directory = output_directory("alignment_without_target");
    let input = "test_files/twin_ari_chrX_146823507_146823598_no_ts.toml";
    let svg_path = directory.join("out.svg");

    // The plain text output shows nothing.
    assert_success(&format!("-i {input} -t"))?;
    assert_error(
        &format!("-i {input} -s {}", svg_path.display()),
        "Error creating SVG",
    )?;

    let output = assert_success(&format!("-i {input} -s {} -r -p", svg_path.display()))?;
    let svg = assert_svg(&svg_path);
    assert!(
        svg.contains("Alignment is incomplete"),
        "the SVG does not show the error\n{}",
        stderr(&output)
    );
    assert!(read(&directory.join("out.png")).starts_with(b"\x89PNG"));
    Ok(())
}

/// A no-TS alignment without target makes the SVG fail, but `-r` renders the error.
#[test]
fn no_ts_alignment_without_target() -> Result<()> {
    let directory = output_directory("no_ts_alignment_without_target");
    let input = "test_files/twin_ari_chrX_146823507_146823598.toml";
    let no_ts = "test_files/twin_ari_chrX_146823507_146823598_no_ts.toml";
    let svg_path = directory.join("out.svg");

    assert_error(
        &format!("-i {input} -n {no_ts} -s {}", svg_path.display()),
        "Error creating SVG",
    )?;
    assert_success(&format!(
        "-i {input} -n {no_ts} -s {} -r",
        svg_path.display()
    ))?;
    assert!(assert_svg(&svg_path).contains("No-TS alignment is incomplete"));

    // The plain text output skips the no-TS alignment.
    let output = assert_success(&format!("-i {input} -n {no_ts} -t"))?;
    assert_eq!(stdout(&output).matches("Switch process:").count(), 3);
    Ok(())
}

/// The SVG cannot be written into a directory that does not exist.
#[test]
fn unwritable_svg_path() -> Result<()> {
    assert_error(
        &format!("-i {TS} -s target/tsalign-show-e2e/does/not/exist/out.svg"),
        "out.svg",
    )
}
