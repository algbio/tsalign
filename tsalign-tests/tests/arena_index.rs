//! `--arena-index` chooses the integer type of the node indices of the A* search.
//! It limits how many nodes the search can store at the same time, but does not change the alignment.
//!
//! A full arena itself cannot be tested at the CLI level: it needs more than 4294967295 nodes at the same time
//! (hundreds of GB of memory), and the binary has no smaller index type.
//! The libraries test it with 16-bit indices instead
//! (`lib_tsalign` and `lib_ts_chainalign`: `full_arena_is_an_error`),
//! and the CLI's error message is tested in `tsalign/src/util.rs`.

use std::fs;

use chain_common::{CONFIG, run_tsalign_output, scratch_path};

mod chain_common;

/// Run `tsalign align` with the given arguments and return the lines of the statistics output that describe the
/// alignment and the search, without the times.
fn align(args: &[&str]) -> String {
    let output = scratch_path("arena-index.toml");
    let mut all_args: Vec<String> = vec![
        "align".into(),
        "-l".into(),
        "warn".into(),
        "-o".into(),
        output.display().to_string(),
    ];
    all_args.extend(args.iter().map(|arg| arg.to_string()));

    let result = run_tsalign_output(&all_args);
    assert!(
        result.status.success(),
        "tsalign {all_args:?} failed ({})\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );
    let statistics = fs::read_to_string(&output).unwrap();
    let lines: Vec<_> = statistics
        .lines()
        .filter(|line| {
            [
                "alignment = ",
                "cost = ",
                "opened_nodes = ",
                "closed_nodes = ",
                "chaining_rounds = ",
                "chaining_opened_nodes = ",
                "chaining_closed_nodes = ",
            ]
            .iter()
            .any(|prefix| line.starts_with(prefix))
        })
        .collect();
    assert!(
        lines.iter().any(|line| line.starts_with("alignment = ")),
        "no alignment in the output\n{statistics}"
    );
    lines.join("\n")
}

/// Check that `--arena-index u32` and `--arena-index u64` give the same output as the default.
fn check_same_output(args: &[&str]) {
    let default = align(args);
    for arena_index in ["u32", "u64"] {
        let mut args = args.to_vec();
        args.extend(["--arena-index", arena_index]);
        assert_eq!(align(&args), default, "--arena-index {arena_index}");
    }
}

#[test]
fn classic_aligner_gives_same_alignment_with_u64() {
    for input in [
        "test_files/twin_100_0.01.fa",
        "test_files/twin_10_ts_only.fa",
        "test_files/twin_100_gap10.fa",
    ] {
        check_same_output(&["-p", input]);
    }
}

#[test]
fn chain_aligner_gives_same_alignment_with_u64() {
    let cache_directory = scratch_path("arena-index-cache");
    fs::create_dir_all(&cache_directory).unwrap();
    let cache_directory = cache_directory.display().to_string();
    for (input, max_anchor_mutations) in [
        ("test_files/twin_100_0.01.fa", "0"),
        ("test_files/twin_100_0.01.fa", "1"),
        ("test_files/twin_10_ts_only.fa", "0"),
    ] {
        check_same_output(&[
            "--alignment-method",
            "a-star-chain-ts",
            "-c",
            CONFIG,
            "--alphabet",
            "dna",
            "--cache-directory",
            &cache_directory,
            "--max-anchor-mutations",
            max_anchor_mutations,
            "-p",
            input,
        ]);
    }
}

#[test]
fn help_documents_arena_index() {
    let output = run_tsalign_output(&["align".into(), "--help".into()]);
    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("--arena-index <ARENA_INDEX>"), "{help}");
    assert!(help.contains("[default: u32]"), "{help}");
    assert!(help.contains("- u32: "), "{help}");
    assert!(help.contains("- u64: "), "{help}");
}

#[test]
fn invalid_arena_index_is_rejected() {
    let output = run_tsalign_output(
        &[
            "align",
            "-p",
            "test_files/twin_10_eq.fa",
            "--arena-index",
            "u16",
        ]
        .map(String::from),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panicked"), "{stderr}");
    // Clap exits with code 2 on invalid arguments.
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("invalid value 'u16'"), "{stderr}");
}
