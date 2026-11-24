//! `tsalign align` and `tsalign preprocess` reject invalid arguments, and arguments that would be ignored,
//! with an error message and exit code 1, instead of panicking or silently ignoring them (issue #249).

use anyhow::Result;
use util::{run_in_repo_root, run_in_repo_root_output};

mod util;

const A_STAR: &str = "align -p test_files/twin_a.fa --alignment-method a-star-template-switch";
/// tschainalign without a cache directory, which the tests add themselves.
const CHAIN_WITHOUT_CACHE: &str = "align -p test_files/twin_a.fa --alignment-method a-star-chain-ts \
                                   -c test_files/config/chainalign --alphabet dna";
const CHAIN: &str = "align -p test_files/twin_a.fa --alignment-method a-star-chain-ts \
                     -c test_files/config/chainalign --alphabet dna --cache-directory tmp/cli_invalid_arguments";

/// Run tsalign and assert that it rejected the arguments with an error that contains all `expected` strings.
fn assert_rejected(args: &str, expected: &[&str]) -> Result<()> {
    let output = run_in_repo_root_output(args)?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("panicked"),
        "`{args}`: tsalign panicked instead of reporting an error\nstderr:\n{stderr}"
    );
    assert_eq!(
        output.status.code(),
        Some(1),
        "`{args}`: unexpected exit status\nstderr:\n{stderr}"
    );
    for expected in expected {
        assert!(
            stderr.contains(expected),
            "`{args}`: the error does not contain {expected:?}\nstderr:\n{stderr}"
        );
    }
    Ok(())
}

/// Run tsalign and assert that it succeeded. The output is printed by [`run_in_repo_root`].
fn assert_accepted(args: &str) -> Result<()> {
    run_in_repo_root(args)
}

#[test]
fn a_star_rejects_chain_options() -> Result<()> {
    for option in [
        "-k 9",
        "--max-anchor-mutations 0",
        "--cache-directory tmp/cli_invalid_arguments",
        "--force-no-preprocessing",
        "--inexact-chaining-lower-bound max-match-run",
        "--max-chaining-successors 5",
        "--max-exact-cost-function-cost 0",
        "--chaining-closed-list special",
        "--chaining-open-list linear-heap",
        "--no-anchor-domination-pruning",
    ] {
        let name = option.split_whitespace().next().unwrap();
        assert_rejected(&format!("{A_STAR} {option}"), &[name, "a-star-chain-ts"])?;
    }
    Ok(())
}

#[test]
fn chain_rejects_a_star_options() -> Result<()> {
    for option in [
        "--ts-min-length-strategy lookahead",
        "--ts-chaining-strategy none",
        "--ts-total-length-strategy maximise",
        "--ts-descendant-strategy allow-any",
        "--allow-ts-14-out-of-range",
        "--ts-uncertainty-range-extension-mode none",
        "--no-ts",
        "--force-label-correcting",
        "--dont-extend-beyond-range",
    ] {
        let name = option.split_whitespace().next().unwrap();
        assert_rejected(
            &format!("{CHAIN} {option}"),
            &[name, "a-star-template-switch"],
        )?;
    }
    Ok(())
}

#[test]
fn a_star_rejects_ts_options_with_no_ts() -> Result<()> {
    for option in [
        "--ts-min-length-strategy none",
        "--ts-chaining-strategy none",
        "--ts-total-length-strategy none",
        "--ts-descendant-strategy allow-only-all-equal",
        "--allow-ts-14-out-of-range --rq-ranges R2..20Q2..20",
        "--ts-uncertainty-range-extension-mode equal-cost",
    ] {
        let name = option.split_whitespace().next().unwrap();
        assert_rejected(&format!("{A_STAR} --no-ts {option}"), &[name, "--no-ts"])?;
    }
    Ok(())
}

#[test]
fn a_star_rejects_dont_extend_beyond_range_without_range() -> Result<()> {
    assert_rejected(
        &format!("{A_STAR} --dont-extend-beyond-range"),
        &["--dont-extend-beyond-range", "range"],
    )?;
    assert_accepted(&format!(
        "{A_STAR} --dont-extend-beyond-range --rq-ranges R2..20Q2..20"
    ))?;
    assert_accepted(&format!(
        "{A_STAR} --dont-extend-beyond-range --reference-offset 2"
    ))
}

#[test]
fn chain_rejects_inexact_options_with_exact_anchors() -> Result<()> {
    for option in [
        "--inexact-chaining-lower-bound max-match-run",
        "--no-anchor-domination-pruning",
    ] {
        let name = option.split_whitespace().next().unwrap();
        assert_rejected(
            &format!("{CHAIN} {option}"),
            &[name, "--max-anchor-mutations"],
        )?;
        assert_rejected(
            &format!("{CHAIN} --max-anchor-mutations 0 {option}"),
            &[name, "--max-anchor-mutations"],
        )?;
    }
    Ok(())
}

#[test]
fn preprocess_rejects_inexact_options_with_exact_anchors() -> Result<()> {
    assert_rejected(
        "preprocess -c test_files/config/chainalign --cache-directory tmp/cli_invalid_arguments \
         -k 4 --inexact-chaining-lower-bound max-match-run 8",
        &["--inexact-chaining-lower-bound", "--max-anchor-mutations"],
    )
}

#[test]
fn preprocess_rejects_zero_max_length() -> Result<()> {
    assert_rejected(
        "preprocess -c test_files/config/chainalign --cache-directory tmp/cli_invalid_arguments -k 4 0",
        &["max_length"],
    )
}

#[test]
fn rejects_malformed_rq_ranges() -> Result<()> {
    for method in [A_STAR, CHAIN] {
        for rq_ranges in ["X0..5", "R0-5", "R0..", "R0..5R1..6", "Q..5"] {
            assert_rejected(
                &format!("{method} --rq-ranges {rq_ranges}"),
                &["--rq-ranges"],
            )?;
        }
    }
    Ok(())
}

#[test]
fn rejects_redundant_rq_ranges() -> Result<()> {
    for method in [A_STAR, CHAIN] {
        assert_rejected(
            &format!("{method} --rq-ranges R0..5 --reference-offset 1"),
            &["--rq-ranges", "--reference-offset"],
        )?;
        assert_rejected(
            &format!("{method} --rq-ranges Q0..5 --query-limit 4"),
            &["--rq-ranges", "--query-limit"],
        )?;
    }
    Ok(())
}

#[test]
fn rejects_range_outside_of_sequences() -> Result<()> {
    // The reference has length 30 and the query has length 29.
    for method in [A_STAR, CHAIN] {
        for range in [
            "--reference-limit 31",
            "--query-offset 30",
            "--rq-ranges R0..30Q0..31",
            "--reference-offset 10 --reference-limit 5",
        ] {
            assert_rejected(&format!("{method} {range}"), &["range"])?;
        }
    }
    Ok(())
}

/// The classic aligner does not support TSMs that start or end outside of the alignment range yet.
/// It used to panic in this case.
#[test]
fn a_star_rejects_allow_ts_14_out_of_range() -> Result<()> {
    for range in [
        "",
        "--rq-ranges R2..20Q2..20",
        "--reference-offset 40 --reference-limit 41",
    ] {
        assert_rejected(
            &format!("{A_STAR} --allow-ts-14-out-of-range {range}"),
            &["--allow-ts-14-out-of-range", "not supported"],
        )?;
    }
    Ok(())
}

#[test]
fn chain_rejects_non_dna_characters() -> Result<()> {
    assert_rejected(
        "align -p test_files/twin_n.fa --alignment-method a-star-chain-ts -c test_files/config/chainalign \
         --alphabet dna-n --cache-directory tmp/cli_invalid_arguments",
        &["N"],
    )
}

#[test]
fn chain_rejects_missing_cache_with_force_no_preprocessing() -> Result<()> {
    assert_rejected(
        &format!(
            "{CHAIN_WITHOUT_CACHE} --cache-directory tmp/cli_invalid_arguments/empty --force-no-preprocessing"
        ),
        &["--force-no-preprocessing"],
    )
}
