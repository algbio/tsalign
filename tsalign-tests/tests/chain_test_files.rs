//! End-to-end tests of tschainalign on the pair fasta files in `test_files`.
//!
//! Only files with fewer than 200 characters per sequence (after skipping gap characters) and only the characters `ACGT` are used.
//! See [`chain_common`] for the checks.

use chain_common::{Case, read_pair_fasta, run_case};

mod chain_common;

fn test_file(path: &str) {
    test_file_skipping(path, b"");
}

/// Test a file whose sequences contain alignment gap characters `-`, which are skipped.
fn test_file_skipping(path: &str, skip_characters: &[u8]) {
    let (reference, query) = read_pair_fasta(path, skip_characters);
    run_case(&Case::new(path, reference, query));
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
    twin_both_descendants => "test_files/twin_both_descendants.fa",
    twin_double_gap => "test_files/twin_double_gap.fa",
    twin_equal_cost_ignore_geometry => "test_files/twin_equal_cost_ignore_geometry.fa",
    twin_in_place_inversion => "test_files/twin_in_place_inversion.fa",
    twin_in_place_inversion_1s => "test_files/twin_in_place_inversion_1S.fa",
    twin_in_place_inversion_2s => "test_files/twin_in_place_inversion_2S.fa",
    twin_in_place_inversion_large => "test_files/twin_in_place_inversion_large.fa",
    twin_intra_forward_jump => "test_files/twin_intra_forward_jump.fa",
    twin_min_length_test1 => "test_files/twin_min_length_test1.fa",
    twin_min_length_test2 => "test_files/twin_min_length_test2.fa",
    twin_perf => "test_files/twin_perf.fa",
    twin_range => "test_files/twin_range.fa",
    twin_repeat_rq => "test_files/twin_repeat_rq.fa",
    twin_repeat_rr => "test_files/twin_repeat_rr.fa",
    twin_repeat_rr2 => "test_files/twin_repeat_rr2.fa",
    twin_show_repeated_chars => "test_files/twin_show_repeated_chars.fa",
    twin_show_ts_copy => "test_files/twin_show_ts_copy.fa",
    twin_show_ts_indel1 => "test_files/twin_show_ts_indel1.fa",
    twin_show_ts_indel2 => "test_files/twin_show_ts_indel2.fa",
    twin_show_ts_indel3 => "test_files/twin_show_ts_indel3.fa",
    twin_show_ts_optional => "test_files/twin_show_ts_optional.fa",
    twin_show_ts_small => "test_files/twin_show_ts_small.fa",
    twin_show_ts_small2 => "test_files/twin_show_ts_small2.fa",
    twin_show_ts_small3 => "test_files/twin_show_ts_small3.fa",
    twin_show_ts_small_indel1 => "test_files/twin_show_ts_small_indel1.fa",
    twin_show_ts_small_indel2 => "test_files/twin_show_ts_small_indel2.fa",
    twin_show_ts_small_rr2 => "test_files/twin_show_ts_small_rr2.fa",
    twin_show_ts_small_rr_qq => "test_files/twin_show_ts_small_rr_qq.fa",
    twin_small_range => "test_files/twin_small_range.fa",
    twin_two_step_modification => "test_files/twin_two_step_modification.fa",
    twin_variable_length_ts => "test_files/twin_variable_length_ts.fa",
}
