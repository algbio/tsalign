//! Without `-c`, `tsalign align --alignment-method a-star-chain-ts` and `tsalign preprocess` use the built-in
//! default config of tschainalign, which equals `test_files/config/chainalign/config.tsa`,
//! instead of `sample_tsa_config`, which tschainalign does not support.

use std::fs;

use chain_common::{CONFIG, run_tsalign_raw, scratch_path};

mod chain_common;

const INPUT: &str = "test_files/twin_10_ts_only.fa";

/// Run `tsalign align` with tschainalign and the given extra arguments, and return the statistics output.
fn align(cache_directory: &str, extra_args: &[&str]) -> String {
    let output = scratch_path("default-config.toml");
    let mut args: Vec<String> = [
        "align",
        "--alignment-method",
        "a-star-chain-ts",
        "--alphabet",
        "dna",
        "-l",
        "warn",
        "-k",
        "4",
        "--cache-directory",
        cache_directory,
        "-p",
        INPUT,
        "-o",
        &output.display().to_string(),
    ]
    .map(String::from)
    .into();
    args.extend(extra_args.iter().map(|arg| arg.to_string()));

    let (status, stderr) = run_tsalign_raw(&args);
    assert!(
        status.success(),
        "tsalign align {extra_args:?} failed ({status})\nstderr:\n{stderr}"
    );
    fs::read_to_string(&output).unwrap()
}

fn cost(statistics: &str) -> &str {
    statistics
        .lines()
        .find_map(|line| line.strip_prefix("cost = "))
        .unwrap_or_else(|| panic!("no cost in the output\n{statistics}"))
}

#[test]
fn align_without_config_uses_chainalign_config() {
    let cache_directory = scratch_path("default-config-cache");
    fs::create_dir_all(&cache_directory).unwrap();
    let cache_directory = cache_directory.display().to_string();

    // The test runs in the repository root, where `sample_tsa_config` exists, so this fails if align uses it.
    let default = align(&cache_directory, &[]);
    // The cache file name depends on the costs, so this fails if the built-in costs differ from CONFIG.
    let explicit = align(
        &cache_directory,
        &["-c", CONFIG, "--force-no-preprocessing"],
    );
    assert_eq!(cost(&default), cost(&explicit));
}

#[test]
fn preprocess_without_config_uses_chainalign_config() {
    let cache_directory = scratch_path("default-config-preprocess-cache");
    fs::create_dir_all(&cache_directory).unwrap();
    let cache_directory = cache_directory.display().to_string();

    let args: Vec<String> = [
        "preprocess",
        "-l",
        "warn",
        "--alphabet",
        "dna",
        "--cache-directory",
        &cache_directory,
        "-k",
        "4",
        "--max-anchor-mutations",
        "0",
        "16",
    ]
    .map(String::from)
    .into();
    let (status, stderr) = run_tsalign_raw(&args);
    assert!(
        status.success(),
        "tsalign preprocess failed ({status})\nstderr:\n{stderr}"
    );

    // Uses the cache written with the built-in costs, which fails if they differ from CONFIG.
    align(
        &cache_directory,
        &["-c", CONFIG, "--force-no-preprocessing"],
    );
}
