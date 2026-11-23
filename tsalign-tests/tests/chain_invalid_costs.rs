//! `tsalign align --alignment-method a-star-chain-ts` rejects cost configs that tschainalign cannot use
//! with an error message, like `tsalign preprocess` does, instead of panicking.

use std::{fs, path::Path};

use chain_common::{CONFIG, repo_root, run_tsalign_raw, scratch_path};

mod chain_common;

/// Run `tsalign align` with tschainalign and the given config directory,
/// and assert that it exits with an error mentioning `expected_message`, without a panic and without output.
fn assert_align_rejects_config(description: &str, config_directory: &Path, expected_message: &str) {
    let input = scratch_path("invalid-costs.fa");
    fs::write(
        &input,
        ">reference\nACGTACGTTGCAACGTTGCAAGGCTT\n>query\nACGTACGATGCAACGTTCCAAGGCTT\n",
    )
    .unwrap();
    let cache_directory = scratch_path("invalid-costs-cache");
    let output = scratch_path("invalid-costs.toml");

    let args: Vec<String> = [
        "align",
        "--alignment-method",
        "a-star-chain-ts",
        "--alphabet",
        "dna",
        "-c",
        &config_directory.display().to_string(),
        "-l",
        "warn",
        "-k",
        "6",
        "--max-anchor-mutations",
        "0",
        "--cache-directory",
        &cache_directory.display().to_string(),
        "-p",
        &input.display().to_string(),
        "-o",
        &output.display().to_string(),
    ]
    .map(String::from)
    .into();

    let (status, stderr) = run_tsalign_raw(&args);
    assert!(
        !status.success(),
        "{description}: tsalign accepted the config\nstderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked"),
        "{description}: tsalign panicked instead of reporting an error ({status})\nstderr:\n{stderr}"
    );
    assert_eq!(
        status.code(),
        Some(1),
        "{description}: unexpected exit status\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains(expected_message),
        "{description}: the error does not contain {expected_message:?}\nstderr:\n{stderr}"
    );
    assert!(
        !output.exists(),
        "{description}: align wrote an alignment to {output:?}"
    );
}

/// A `Length` cost that is finite outside of its zero range cannot be converted into tschainalign costs.
#[test]
fn align_rejects_length_cost_without_single_zero_range() {
    let config = fs::read_to_string(repo_root().join(CONFIG).join("config.tsa")).unwrap();
    let valid_length = "Length\n   0 10 200\n inf  0 inf\n";
    assert!(
        config.contains(valid_length),
        "unexpected {CONFIG}/config.tsa"
    );
    let config = config.replace(valid_length, "Length\n   0 10 200\n inf  0   1\n");

    let config_directory = scratch_path("invalid-length-config");
    fs::create_dir_all(&config_directory).unwrap();
    fs::write(config_directory.join("config.tsa"), config).unwrap();

    assert_align_rejects_config(
        "invalid Length costs",
        &config_directory,
        "the Length costs must be zero within a single range and infinite outside of it",
    );
}

/// A missing config file is reported as an error.
#[test]
fn align_rejects_missing_config() {
    let config_directory = scratch_path("missing-config");
    assert_align_rejects_config(
        "missing config",
        &config_directory,
        "Error opening config file",
    );
}
