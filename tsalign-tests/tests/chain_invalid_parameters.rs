//! tschainalign rejects unsupported combinations of `k` and `max_anchor_mutations` with an error message,
//! both in `tsalign align --alignment-method a-star-chain-ts` and in `tsalign preprocess`.
//!
//! Without validation, such parameters either panic on an arithmetic overflow, or, without overflow checks,
//! silently compute lower bounds from wrapped-around values (e.g. `k - 1 - max_anchor_mutations`),
//! write them into the cache and report wrong alignments.

use std::fs;

use chain_common::{CONFIG, run_tsalign_raw, scratch_path};

mod chain_common;

/// Unsupported `(k, max_anchor_mutations)` combinations. `None` means that `k` is inferred from the sequence lengths.
const INVALID_PARAMETERS: &[(Option<u32>, u8)] = &[
    // max_anchor_mutations must be less than k.
    (Some(6), 6),
    (Some(6), 7),
    (Some(1), 1),
    // k must be at least one.
    (Some(0), 0),
    // Exact anchors (max_anchor_mutations = 0) support only k <= 64.
    (Some(65), 0),
    // Inexact anchors support only k <= 255.
    (Some(256), 1),
    // The inferred k grows with the number of anchor mutations, and exceeds 255 here.
    (None, 200),
];

/// Assert that tsalign rejected the parameters with an error instead of a panic or a silent success,
/// and that it did not write any lower bounds into the cache.
fn assert_rejected(description: &str, args: &[String], cache_directory: &std::path::Path) {
    let (status, stderr) = run_tsalign_raw(args);
    let cache_files: Vec<_> = fs::read_dir(cache_directory)
        .map(|entries| entries.map(|entry| entry.unwrap().file_name()).collect())
        .unwrap_or_default();

    assert!(
        !status.success(),
        "{description}: tsalign accepted the parameters\nstderr:\n{stderr}"
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
        stderr.contains("max_anchor_mutations") || stderr.contains("k = "),
        "{description}: the error does not explain the problem with k or max_anchor_mutations\nstderr:\n{stderr}"
    );
    assert!(
        cache_files.is_empty(),
        "{description}: tsalign wrote cache files {cache_files:?}"
    );
}

fn k_args(k: Option<u32>, max_anchor_mutations: u8) -> Vec<String> {
    let mut args = Vec::new();
    if let Some(k) = k {
        args.extend(["-k".to_string(), k.to_string()]);
    }
    args.extend([
        "--max-anchor-mutations".to_string(),
        max_anchor_mutations.to_string(),
    ]);
    args
}

#[test]
fn align_rejects_invalid_k_and_max_anchor_mutations() {
    let input = scratch_path("invalid-parameters.fa");
    fs::write(
        &input,
        ">reference\nACGTACGTTGCAACGTTGCAAGGCTT\n>query\nACGTACGATGCAACGTTCCAAGGCTT\n",
    )
    .unwrap();

    for &(k, max_anchor_mutations) in INVALID_PARAMETERS {
        let cache_directory = scratch_path("invalid-parameters-cache");
        fs::create_dir_all(&cache_directory).unwrap();
        let output = scratch_path("invalid-parameters.toml");
        let mut args: Vec<String> = [
            "align",
            "--alignment-method",
            "a-star-chain-ts",
            "--alphabet",
            "dna",
            "-c",
            CONFIG,
            "-l",
            "warn",
        ]
        .map(String::from)
        .into();
        args.extend(k_args(k, max_anchor_mutations));
        args.extend([
            "--cache-directory".to_string(),
            cache_directory.display().to_string(),
            "-p".to_string(),
            input.display().to_string(),
            "-o".to_string(),
            output.display().to_string(),
        ]);

        assert_rejected(
            &format!("align with k = {k:?}, max_anchor_mutations = {max_anchor_mutations}"),
            &args,
            &cache_directory,
        );
        assert!(!output.exists(), "align wrote an alignment to {output:?}");
    }
}

#[test]
fn preprocess_rejects_invalid_k_and_max_anchor_mutations() {
    for &(k, max_anchor_mutations) in INVALID_PARAMETERS {
        let cache_directory = scratch_path("invalid-parameters-cache");
        let mut args: Vec<String> = [
            "preprocess",
            "--alphabet",
            "dna",
            "-c",
            CONFIG,
            "-l",
            "warn",
        ]
        .map(String::from)
        .into();
        args.extend(k_args(k, max_anchor_mutations));
        args.extend([
            "--cache-directory".to_string(),
            cache_directory.display().to_string(),
            "16".to_string(),
        ]);

        assert_rejected(
            &format!("preprocess with k = {k:?}, max_anchor_mutations = {max_anchor_mutations}"),
            &args,
            &cache_directory,
        );
    }
}

/// The largest supported parameters are still accepted.
#[test]
fn preprocess_accepts_boundary_k_and_max_anchor_mutations() {
    for (k, max_anchor_mutations) in [(1, 0), (64, 0), (2, 1), (255, 254)] {
        let cache_directory = scratch_path("boundary-parameters-cache");
        let mut args: Vec<String> = [
            "preprocess",
            "--alphabet",
            "dna",
            "-c",
            CONFIG,
            "-l",
            "warn",
        ]
        .map(String::from)
        .into();
        args.extend(k_args(Some(k), max_anchor_mutations));
        args.extend([
            "--cache-directory".to_string(),
            cache_directory.display().to_string(),
            "16".to_string(),
        ]);

        let (status, stderr) = run_tsalign_raw(&args);
        assert!(
            status.success(),
            "preprocess with k = {k}, max_anchor_mutations = {max_anchor_mutations} failed ({status})\nstderr:\n{stderr}"
        );
        assert!(
            fs::read_dir(&cache_directory).unwrap().next().is_some(),
            "preprocess with k = {k}, max_anchor_mutations = {max_anchor_mutations} wrote no cache files"
        );
    }
}
