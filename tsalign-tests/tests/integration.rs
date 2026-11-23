use anyhow::Result;
use util::{run_in_repo_root, run_in_repo_root_output};

mod util;

#[test]
fn align_default_cfg_twin() -> Result<()> {
    run_in_repo_root("align -p test_files/twin_a.fa")
}

#[test]
fn align_default_cfg_qr() -> Result<()> {
    run_in_repo_root("align -q test_files/query_a.fa -r test_files/reference_a.fa")
}

#[test]
fn align_with_cost_limit() -> Result<()> {
    run_in_repo_root("align -p test_files/twin_100_0.01.fa --cost-limit 0")
}

#[test]
fn align_with_memory_limit() -> Result<()> {
    run_in_repo_root("align -p test_files/twin_100_0.01.fa --memory-limit 1000")
}

#[test]
fn align_with_embedded_rq_ranges() -> Result<()> {
    run_in_repo_root("align -p test_files/twin_embedded.fa --use-embedded-rq-ranges")
}

/// Issue #248: the lower bound chaining strategy is disabled, so the CLI rejects it with a usage error.
#[test]
fn test_align_lower_bound_chaining_is_rejected() -> Result<()> {
    let output = run_in_repo_root_output(
        "align -p test_files/twin_10_eq.fa --ts-chaining-strategy lower-bound",
    )?;
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        !stderr.contains("panicked"),
        "tsalign panicked instead of rejecting the chaining strategy"
    );
    // Clap exits with code 2 on invalid arguments.
    assert_eq!(output.status.code(), Some(2));
    assert!(
        stderr.contains("invalid value 'lower-bound'"),
        "the error does not name the rejected value"
    );
    Ok(())
}
