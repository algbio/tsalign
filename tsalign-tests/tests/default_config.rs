//! Without `-c`, `tsalign` uses the default configs built into the binary (#280),
//! independent of the working directory: the classic aligner the one equal to `sample_tsa_config`,
//! tschainalign the one equal to `test_files/config/chainalign`.

use std::{fs, path::Path};

use chain_common::{CONFIG, repo_root, run_tsalign_output_in, scratch_path};

mod chain_common;

/// An input with a template switch, so that the template switch costs matter.
const INPUT: &str = "test_files/twin_10_ts_only.fa";

/// An empty working directory outside the repository root, so no config can be found relative to it.
fn empty_working_directory() -> std::path::PathBuf {
    let directory = scratch_path("empty-working-directory");
    fs::create_dir_all(&directory).unwrap();
    directory
}

/// Run `tsalign align` in `working_directory` and return the `alignment` and `cost` lines of the statistics output.
fn align(working_directory: &Path, extra_args: &[&str]) -> String {
    let output = scratch_path("default-config.toml");
    let input = repo_root().join(INPUT);
    let mut args: Vec<String> = vec![
        "align".into(),
        "-l".into(),
        "warn".into(),
        "-p".into(),
        input.display().to_string(),
        "-o".into(),
        output.display().to_string(),
    ];
    args.extend(extra_args.iter().map(|arg| arg.to_string()));

    let result = run_tsalign_output_in(&args, working_directory);
    assert!(
        result.status.success(),
        "tsalign align {extra_args:?} failed ({})\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );
    let statistics = fs::read_to_string(&output).unwrap();
    let lines: Vec<_> = statistics
        .lines()
        .filter(|line| line.starts_with("alignment = ") || line.starts_with("cost = "))
        .collect();
    assert_eq!(lines.len(), 3, "unexpected output\n{statistics}");
    lines.join("\n")
}

#[test]
fn classic_align_without_config_outside_repository() {
    let working_directory = empty_working_directory();
    let sample_config = repo_root().join("sample_tsa_config").display().to_string();

    let default = align(&working_directory, &[]);
    let explicit = align(&working_directory, &["-c", &sample_config]);
    assert_eq!(default, explicit);
}

#[test]
fn classic_align_without_config_supports_dna_alphabet() {
    // The built-in default config contains `N`, which the DNA alphabet lacks.
    // The input contains only A, C, G and T, so the result is the same as with the default alphabet `dna-n`.
    let working_directory = empty_working_directory();
    let dna = align(&working_directory, &["--alphabet", "dna"]);
    let dna_n = align(&working_directory, &["--alphabet", "dna-n"]);
    assert_eq!(dna, dna_n);
}

#[test]
fn chain_align_without_config_outside_repository() {
    let working_directory = empty_working_directory();
    let cache_directory = scratch_path("default-config-cache");
    fs::create_dir_all(&cache_directory).unwrap();
    let cache_directory = cache_directory.display().to_string();
    let chain_config = repo_root().join(CONFIG).display().to_string();
    let chain_args = [
        "--alignment-method",
        "a-star-chain-ts",
        "--alphabet",
        "dna",
        "-k",
        "4",
        "--cache-directory",
        &cache_directory,
    ];

    let default = align(&working_directory, &chain_args);
    // The cache file name depends on the costs, so this fails if the built-in costs differ from CONFIG.
    let mut explicit_args = chain_args.to_vec();
    explicit_args.extend(["-c", &chain_config, "--force-no-preprocessing"]);
    let explicit = align(&working_directory, &explicit_args);
    assert_eq!(default, explicit);
}
