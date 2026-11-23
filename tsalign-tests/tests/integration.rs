use anyhow::Result;
use util::run_in_repo_root;

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
