//! End-to-end tests of tschainalign on the pair fasta files in `test_files`.
//!
//! Only files with fewer than 200 characters per sequence (after skipping gap characters) and only the characters `ACGT` are used.
//!
//! Some files are too slow for some parameter combinations (issue #223).
//! Files that do not test a specific edge case (rendering demos, benchmarks) are tested without the slow combinations.
//! Files on which even the remaining combinations are too slow, or which test a specific edge case,
//! are replaced by shorter files built the same way, which are stored in `test_files/e2e_short` (see the README there).
//! See [`chain_common`] for the checks.

use chain_common::{Case, Params, read_pair_fasta, run_case};

mod chain_common;

fn test_file(path: &str) {
    test_file_skipping(path, b"");
}

/// Test a file whose sequences contain alignment gap characters `-`, which are skipped.
fn test_file_skipping(path: &str, skip_characters: &[u8]) {
    let (reference, query) = read_pair_fasta(path, skip_characters);
    run_case(&Case::new(path, reference, query));
}

/// Test a file with the default grid without the given parameter combinations.
///
/// Used for files on which these combinations take minutes (issue #223).
/// These files were made as rendering demos or as a benchmark for the classic aligner, not for chaining,
/// so little chaining coverage is lost.
fn case_without(path: &str, slow: &[Params]) -> Case {
    let (reference, query) = read_pair_fasta(path, b"");
    let case = Case::new(path, reference, query);
    let grid: Vec<Params> = case
        .grid
        .iter()
        .copied()
        .filter(|params| !slow.contains(params))
        .collect();
    case.with_grid(grid)
}

#[test]
fn ari_41_short() {
    test_file_skipping("test_files/ari-41-short.fa", b"-");
}

#[test]
fn twin_minus() {
    test_file_skipping("test_files/twin_minus.fa", b"-");
}

#[test]
fn homo_sapiens_12_74743744_74973891_sub() {
    test_file_skipping("test_files/homo_sapiens.12.74743744.74973891.sub.fa", b"-");
}

/// Rendering demo for an optional template switch.
/// Runs with inexact anchors take between five and more than ten minutes each on its long homopolymer runs,
/// and the classic aligner with template switches takes about two minutes (issue #223).
#[test]
fn twin_show_ts_optional() {
    run_case(
        &case_without(
            "test_files/twin_show_ts_optional.fa",
            &[(6, 1), (9, 1), (12, 1)],
        )
        .without_ts_lower_bound(),
    );
}

/// Rendering demo for RR and QQ template switches. k = 6 with one anchor mutation takes about two minutes (issue #223).
#[test]
fn twin_show_ts_small_rr_qq() {
    run_case(&case_without(
        "test_files/twin_show_ts_small_rr_qq.fa",
        &[(6, 1)],
    ));
}

/// Rendering demo for template switch inners. k = 6 with one anchor mutation takes about 20 seconds (issue #223).
#[test]
fn twin_show_ts_small3() {
    run_case(&case_without(
        "test_files/twin_show_ts_small3.fa",
        &[(6, 1)],
    ));
}

/// Rendering demo for indels next to a template switch. k = 6 with one anchor mutation takes more than a minute (issue #223).
#[test]
fn twin_show_ts_indel2() {
    run_case(&case_without(
        "test_files/twin_show_ts_indel2.fa",
        &[(6, 1)],
    ));
}

/// Demo for visualising the uncertainty of template switches of variable length.
/// Most runs with inexact anchors take 20 to 50 seconds each on its homopolymer runs (issue #223).
#[test]
fn twin_variable_length_ts() {
    run_case(&case_without(
        "test_files/twin_variable_length_ts.fa",
        &[(6, 1), (12, 1), (10, 2), (13, 2)],
    ));
}

/// Performance example (named `perf`) for the classic aligner. k = 6 with one anchor mutation takes more than a minute (issue #223).
#[test]
fn twin_perf() {
    run_case(&case_without("test_files/twin_perf.fa", &[(6, 1)]));
}

macro_rules! test_files {
    ($($name:ident => $path:literal,)*) => {
        $(
            #[test]
            fn $name() {
                test_file($path);
            }
        )*
    };
}

test_files! {
    error_pair => "test_files/error.pair.fa",
    error2_pair => "test_files/error2.pair.fa",
    twin_100_0_01 => "test_files/twin_100_0.01.fa",
    twin_100_0_02 => "test_files/twin_100_0.02.fa",
    twin_100_0_03 => "test_files/twin_100_0.03.fa",
    twin_100_eq => "test_files/twin_100_eq.fa",
    twin_100_gap10 => "test_files/twin_100_gap10.fa",
    twin_10_0_1 => "test_files/twin_10_0.1.fa",
    twin_10_0_2 => "test_files/twin_10_0.2.fa",
    twin_10_0_3 => "test_files/twin_10_0.3.fa",
    twin_10_eq => "test_files/twin_10_eq.fa",
    twin_10_no_anchors => "test_files/twin_10_no_anchors.fa",
    twin_10_repetitive => "test_files/twin_10_repetitive.fa",
    twin_10_ts => "test_files/twin_10_ts.fa",
    twin_10_ts_only => "test_files/twin_10_ts_only.fa",
    twin_20_badends => "test_files/twin_20_badends.fa",
    twin_a => "test_files/twin_a.fa",
    twin_both_descendants => "test_files/e2e_short/twin_both_descendants.fa",
    twin_double_gap => "test_files/twin_double_gap.fa",
    twin_equal_cost_ignore_geometry => "test_files/e2e_short/twin_equal_cost_ignore_geometry.fa",
    twin_in_place_inversion => "test_files/twin_in_place_inversion.fa",
    twin_in_place_inversion_1s => "test_files/twin_in_place_inversion_1S.fa",
    twin_in_place_inversion_2s => "test_files/twin_in_place_inversion_2S.fa",
    twin_in_place_inversion_large => "test_files/twin_in_place_inversion_large.fa",
    twin_intra_forward_jump => "test_files/twin_intra_forward_jump.fa",
    twin_min_length_test1 => "test_files/twin_min_length_test1.fa",
    twin_min_length_test2 => "test_files/twin_min_length_test2.fa",
    twin_range => "test_files/twin_range.fa",
    twin_repeat_rq => "test_files/twin_repeat_rq.fa",
    twin_repeat_rr => "test_files/twin_repeat_rr.fa",
    twin_repeat_rr2 => "test_files/twin_repeat_rr2.fa",
    twin_show_repeated_chars => "test_files/twin_show_repeated_chars.fa",
    twin_show_ts_copy => "test_files/twin_show_ts_copy.fa",
    twin_show_ts_indel1 => "test_files/twin_show_ts_indel1.fa",
    twin_show_ts_indel3 => "test_files/twin_show_ts_indel3.fa",
    twin_show_ts_small => "test_files/twin_show_ts_small.fa",
    twin_show_ts_small2 => "test_files/twin_show_ts_small2.fa",
    twin_show_ts_small_indel1 => "test_files/twin_show_ts_small_indel1.fa",
    twin_show_ts_small_indel2 => "test_files/twin_show_ts_small_indel2.fa",
    twin_show_ts_small_rr2 => "test_files/twin_show_ts_small_rr2.fa",
    twin_small_range => "test_files/twin_small_range.fa",
    twin_two_step_modification => "test_files/twin_two_step_modification.fa",
}
