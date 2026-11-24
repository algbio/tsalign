//! The cost limit and the iteration limit of tschainalign (`tsalign align --alignment-method a-star-chain-ts`),
//! see [issue #259](https://github.com/algbio/tsalign/issues/259).
//!
//! * `--cost-limit` is inclusive and applies to the chain cost: an alignment whose chain cost equals the limit
//!   is returned, otherwise the result reports that the cost limit was exceeded.
//! * `--iteration-limit` limits the number of chaining rounds.
//! * `--memory-limit` is not supported by tschainalign and is rejected with an error.
//!
//! Reaching a limit is a regular result (exit status 0) without an alignment, like for the classic aligner.

use std::fs;

use chain_common::{CONFIG, preprocessed_cache, run_tsalign_output, scratch_path};

mod chain_common;

/// The optimal alignment contains a template switch.
const TS_PAIR: (&str, &str) = (
    "TCTAACTCCAGTCAACTTCCGAGAGCACCA",
    "TCTAACTCCAGGAAGTAAACGAGAGATCCA",
);

/// Finding the optimal chain takes several chaining rounds.
const ROUNDS_PAIR: (&str, &str) = (
    "GGACATTCATCAAGGCCCAGTACCTAATTC",
    "GAACATTCATCTGGGCCTTGTGCCTAATTC",
);

struct Run {
    status: std::process::ExitStatus,
    stdout: String,
    stderr: String,
    /// The toml output, if any.
    output: Option<toml::Table>,
}

fn run(pair: (&str, &str), alignment_method: &str, extra_args: &[&str]) -> Run {
    let input = scratch_path("limits.fa");
    let output = input.with_extension("toml");
    fs::write(
        &input,
        format!(">reference\n{}\n>query\n{}\n", pair.0, pair.1),
    )
    .unwrap();

    let mut args: Vec<String> = [
        "align",
        "--alignment-method",
        alignment_method,
        "--alphabet",
        "dna",
        "-c",
        CONFIG,
        "-l",
        "warn",
    ]
    .map(String::from)
    .into();
    if alignment_method == "a-star-chain-ts" {
        let (k, max_anchor_mutations) = (4, 0);
        let cache = preprocessed_cache(k, max_anchor_mutations).unwrap();
        args.extend([
            "-k".to_string(),
            k.to_string(),
            "--max-anchor-mutations".to_string(),
            max_anchor_mutations.to_string(),
            "--cache-directory".to_string(),
            cache.display().to_string(),
            "--force-no-preprocessing".to_string(),
        ]);
    }
    args.extend(extra_args.iter().map(|arg| arg.to_string()));
    args.extend([
        "-p".to_string(),
        input.display().to_string(),
        "-o".to_string(),
        output.display().to_string(),
    ]);

    let result = run_tsalign_output(&args);
    let toml_output = fs::read_to_string(&output)
        .ok()
        .map(|content| toml::from_str(&content).unwrap());
    let _ = fs::remove_file(&input);
    let _ = fs::remove_file(&output);
    Run {
        status: result.status,
        stdout: String::from_utf8_lossy(&result.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&result.stderr).into_owned(),
        output: toml_output,
    }
}

fn run_chain(pair: (&str, &str), extra_args: &[&str]) -> Run {
    run(pair, "a-star-chain-ts", extra_args)
}

impl Run {
    fn assert_success(&self, description: &str) -> &toml::Table {
        assert!(
            self.status.success(),
            "{description}: tsalign failed ({})\nstderr:\n{}",
            self.status,
            self.stderr
        );
        self.output
            .as_ref()
            .unwrap_or_else(|| panic!("{description}: no toml output"))
    }

    fn result_type(&self, description: &str) -> String {
        self.assert_success(description)["type"]
            .as_str()
            .unwrap()
            .to_string()
    }

    fn astar_result_type(&self, description: &str) -> String {
        self.assert_success(description)["result"]["astar_result_type"]
            .as_str()
            .unwrap()
            .to_string()
    }

    fn chain_cost(&self, description: &str) -> u32 {
        self.assert_success(description)["chain_cost"]
            .as_float()
            .unwrap() as u32
    }

    fn assert_rejected(&self, description: &str, expected_message: &str) {
        assert_eq!(
            self.status.code(),
            Some(1),
            "{description}: unexpected exit status\nstderr:\n{}",
            self.stderr
        );
        assert!(
            !self.stderr.contains("panicked"),
            "{description}: tsalign panicked\nstderr:\n{}",
            self.stderr
        );
        assert!(
            self.stderr.contains(expected_message),
            "{description}: the error does not contain {expected_message:?}\nstderr:\n{}",
            self.stderr
        );
        assert!(
            self.output.is_none(),
            "{description}: tsalign wrote an output"
        );
    }
}

#[test]
fn cost_limit_is_inclusive() {
    for pair in [TS_PAIR, ROUNDS_PAIR] {
        let unlimited = run_chain(pair, &[]);
        assert_eq!(unlimited.result_type("no limit"), "WithTarget");
        let optimum = unlimited.chain_cost("no limit");
        assert!(optimum > 0);

        for limit in [optimum, optimum + 1] {
            let description = format!("{pair:?} with cost limit {limit}");
            let limited = run_chain(pair, &["--cost-limit", &limit.to_string()]);
            assert_eq!(limited.result_type(&description), "WithTarget");
            assert_eq!(limited.chain_cost(&description), optimum, "{description}");
            assert_eq!(
                limited.output.as_ref().unwrap()["cost"],
                unlimited.output.as_ref().unwrap()["cost"],
                "{description}"
            );
        }

        for limit in [0, optimum - 1] {
            let description = format!("{pair:?} with cost limit {limit}");
            let limited = run_chain(pair, &["--cost-limit", &limit.to_string()]);
            assert_eq!(limited.result_type(&description), "WithoutTarget");
            assert_eq!(limited.astar_result_type(&description), "ExceededCostLimit");
            assert_eq!(
                limited.output.as_ref().unwrap()["result"]["cost_limit"].as_integer(),
                Some(i64::from(limit)),
                "{description}"
            );
            assert!(
                limited.stdout.contains("No alignment found")
                    && limited
                        .stdout
                        .contains(&format!("Exceeded cost limit of {limit}")),
                "{description}: stdout:\n{}",
                limited.stdout
            );
        }
    }
}

/// A cost limit larger than the largest chain cost (`u32::MAX`) is no limit.
#[test]
fn huge_cost_limit() {
    let result = run_chain(TS_PAIR, &["--cost-limit", &u64::MAX.to_string()]);
    assert_eq!(result.result_type("huge cost limit"), "WithTarget");
}

#[test]
fn iteration_limit() {
    let unlimited = run_chain(ROUNDS_PAIR, &[]);
    let optimum = unlimited.chain_cost("no limit");

    let description = "iteration limit 0";
    let limited = run_chain(ROUNDS_PAIR, &["--iteration-limit", "0"]);
    assert_eq!(limited.result_type(description), "WithoutTarget");
    assert_eq!(limited.astar_result_type(description), "NoTarget");
    let output = limited.output.as_ref().unwrap();
    assert_eq!(output["exceeded_iteration_limit"].as_integer(), Some(0));
    assert_eq!(output["cost"].as_float(), Some(0.0));
    assert!(
        limited.stdout.contains("No alignment found")
            && limited
                .stdout
                .contains("Exceeded iteration limit of 0, the chain cost is at least 0"),
        "stdout:\n{}",
        limited.stdout
    );

    // The number of rounds needed for the optimal chain.
    let mut rounds = 1;
    loop {
        let description = format!("iteration limit {rounds}");
        let limited = run_chain(ROUNDS_PAIR, &["--iteration-limit", &rounds.to_string()]);
        if limited.result_type(&description) == "WithTarget" {
            assert_eq!(limited.chain_cost(&description), optimum);
            assert!(
                !limited
                    .output
                    .as_ref()
                    .unwrap()
                    .contains_key("exceeded_iteration_limit"),
                "{description}"
            );
            break;
        }
        let output = limited.output.as_ref().unwrap();
        assert_eq!(
            output["exceeded_iteration_limit"].as_integer(),
            Some(rounds),
            "{description}"
        );
        assert!(
            output["cost"].as_float().unwrap() <= f64::from(optimum),
            "{description}: the lower bound exceeds the optimum"
        );
        rounds += 1;
        assert!(rounds < 1000, "the chaining does not terminate");
    }
    assert!(rounds > 1, "the input needs only one chaining round");
}

#[test]
fn memory_limit_is_rejected() {
    run_chain(TS_PAIR, &["--memory-limit", "1000000000"]).assert_rejected(
        "tschainalign with memory limit",
        "--memory-limit is only used by --alignment-method a-star-template-switch",
    );
}

#[test]
fn iteration_limit_is_rejected_by_the_classic_aligner() {
    run(
        TS_PAIR,
        "a-star-template-switch",
        &["--iteration-limit", "10"],
    )
    .assert_rejected(
        "classic aligner with iteration limit",
        "--iteration-limit is only used by --alignment-method a-star-chain-ts",
    );
}
