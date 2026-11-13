//! Hand-crafted end-to-end edge cases for tschainalign with inexact anchors.
//!
//! See [`chain_common`] for the checks that are applied to every case.
//! Some tests additionally assert concrete expectations, such as a cost of zero for identical sequences
//! or that a template switch is found if the anchors support it.

use chain_common::{Case, Params, RunOutcome, run_case, run_chain, run_classic};

mod chain_common;

/// The reverse complement of a DNA sequence.
fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => panic!("unexpected character {}", *other as char),
        })
        .collect()
}

/// A deterministic pseudo-random DNA sequence.
fn random_dna(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(0x2545_F491_4F6C_DD1D);
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 32) as usize % 4]
        })
        .collect()
}

/// Substitute the character at `position` with a different one.
fn substitute(sequence: &mut [u8], position: usize) {
    sequence[position] = match sequence[position] {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        other => panic!("unexpected character {}", other as char),
    };
}

/// Replace `sequence[offset..offset + length]` with the reverse complement of itself (an in-place inversion).
fn invert(sequence: &[u8], offset: usize, length: usize) -> Vec<u8> {
    let mut result = sequence.to_vec();
    result[offset..offset + length]
        .copy_from_slice(&reverse_complement(&sequence[offset..offset + length]));
    result
}

fn concat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

fn chain_result(case: &Case, (k, max_anchor_mutations): Params) -> chain_common::ChainResult {
    match run_chain(&case.reference, &case.query, k, max_anchor_mutations) {
        RunOutcome::Success(result) => result,
        other => panic!(
            "{}: k = {k}, max_anchor_mutations = {max_anchor_mutations}: {other}",
            case.name
        ),
    }
}

fn classic_cost(case: &Case, template_switches: bool) -> u32 {
    match run_classic(&case.reference, &case.query, template_switches) {
        RunOutcome::Success(result) => result.cost,
        other => panic!("{}: classic aligner: {other}", case.name),
    }
}

/// Assert that the alignment has the given cost.
fn assert_cost(case: &Case, params: Params, expected: u32) {
    let result = chain_result(case, params);
    assert_eq!(
        result.cost, expected,
        "{}: k = {}, max_anchor_mutations = {}: expected cost {expected}, got {} with alignment {}",
        case.name, params.0, params.1, result.cost, result.alignment
    );
}

/// Assert that the anchors for `params` support a template switch that makes the alignment cheaper
/// than the optimal alignment without template switches.
fn assert_finds_template_switch(case: &Case, params: Params) {
    let no_ts_cost = classic_cost(case, false);
    let result = chain_result(case, params);
    assert!(
        result.template_switch_amount >= 1 && result.cost < no_ts_cost,
        "{}: k = {}, max_anchor_mutations = {}: expected a template switch that is cheaper than the optimal alignment \
         without template switches (cost {no_ts_cost}), but got cost {} with {} template switches and alignment {}",
        case.name,
        params.0,
        params.1,
        result.cost,
        result.template_switch_amount,
        result.alignment
    );
}

/// Identical sequences align without any cost, for every `k` and number of anchor mutations.
#[test]
fn identical_sequences() {
    let sequence = random_dna(1, 40);
    let case = Case::new("identical_sequences", &sequence, &sequence);
    run_case(&case);
    assert_cost(&case, (6, 0), 0);
    assert_cost(&case, (9, 1), 0);
}

/// A single equal character in both sequences, much shorter than any k.
#[test]
fn single_character_equal() {
    let case = Case::new("single_character_equal", b"A", b"A");
    run_case(&case);
    assert_cost(&case, (4, 0), 0);
}

/// A single differing character in both sequences: one substitution.
#[test]
fn single_character_different() {
    run_case(&Case::new("single_character_different", b"A", b"C"));
}

/// One sequence of length 1 against a sequence of length 20 that contains the character.
#[test]
fn length_one_versus_twenty() {
    let mut long = random_dna(2, 20);
    long[10] = b'G';
    run_case(&Case::new("length_one_versus_twenty", b"G", &long));
}

/// Both sequences are much shorter than k, so there are no anchors at all.
#[test]
fn much_shorter_than_k() {
    run_case(
        &Case::new("much_shorter_than_k", b"ACG", b"ACTG").with_grid([
            (9, 0),
            (12, 0),
            (9, 1),
            (12, 1),
            (13, 2),
        ]),
    );
}

/// The sequences have length exactly k and differ by one substitution in the middle,
/// so there is exactly one inexact anchor candidate and no exact one.
#[test]
fn length_equals_k() {
    let reference = random_dna(3, 9);
    let mut query = reference.clone();
    substitute(&mut query, 4);
    run_case(&Case::new("length_equals_k", &reference, &query).with_grid([(9, 0), (9, 1), (9, 2)]));
}

/// The sequences have length k - 1, so no k-mer fits into them.
#[test]
fn k_is_length_plus_one() {
    let reference = random_dna(4, 8);
    let mut query = reference.clone();
    substitute(&mut query, 3);
    run_case(
        &Case::new("k_is_length_plus_one", &reference, &query).with_grid([(9, 0), (9, 1), (9, 2)]),
    );
}

/// The sequences share no k-mer at all, neither forward nor reverse complement.
#[test]
fn no_shared_kmer() {
    run_case(&Case::new("no_shared_kmer", [b'A'; 20], [b'C'; 20]));
}

/// Two unrelated sequences that share just a single 6-mer in the middle.
#[test]
fn nearly_no_shared_kmer() {
    let shared = b"GATTCA";
    let reference = concat(&[&random_dna(5, 15), shared, &random_dna(6, 15)]);
    let query = concat(&[&random_dna(7, 15), shared, &random_dna(8, 15)]);
    run_case(&Case::new("nearly_no_shared_kmer", &reference, &query));
}

/// Homopolymers of different length: the optimal alignment is a single gap.
#[test]
fn homopolymer_length_difference() {
    run_case(&Case::new(
        "homopolymer_length_difference",
        [b'A'; 12],
        [b'A'; 20],
    ));
}

/// A dinucleotide tandem repeat with a copy number difference of three, between unique flanks.
#[test]
fn dinucleotide_repeat_copy_number() {
    let left = random_dna(9, 10);
    let right = random_dna(10, 10);
    let reference = concat(&[&left, &b"CA".repeat(8), &right]);
    let query = concat(&[&left, &b"CA".repeat(11), &right]);
    run_case(&Case::new(
        "dinucleotide_repeat_copy_number",
        &reference,
        &query,
    ));
}

/// A trinucleotide tandem repeat with a copy number difference of three, between unique flanks.
#[test]
fn trinucleotide_repeat_copy_number() {
    let left = random_dna(11, 10);
    let right = random_dna(12, 10);
    let reference = concat(&[&left, &b"CAG".repeat(6), &right]);
    let query = concat(&[&left, &b"CAG".repeat(9), &right]);
    run_case(&Case::new(
        "trinucleotide_repeat_copy_number",
        &reference,
        &query,
    ));
}

/// The query is the reverse complement of the whole reference,
/// so a single template switch spans both sequences completely.
#[test]
fn full_reverse_complement() {
    let reference = random_dna(13, 30);
    let query = reverse_complement(&reference);
    let case = Case::new("full_reverse_complement", &reference, &query);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

fn inversion_case(name: &str, seed: u64, offset: usize, length: usize) -> Case {
    let reference = random_dna(seed, 40);
    let query = invert(&reference, offset, length);
    Case::new(name, &reference, &query)
}

/// An exact inversion of the minimum template switch length 10 at the start of the sequences.
#[test]
fn inversion_10_start() {
    let case = inversion_case("inversion_10_start", 14, 0, 10);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

/// An exact inversion of the minimum template switch length 10 in the middle of the sequences.
#[test]
fn inversion_10_middle() {
    let case = inversion_case("inversion_10_middle", 15, 15, 10);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

/// An exact inversion of the minimum template switch length 10 at the end of the sequences.
#[test]
fn inversion_10_end() {
    let case = inversion_case("inversion_10_end", 16, 30, 10);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

/// An exact inversion of length 9, one below the minimum template switch length, at the start of the sequences.
#[test]
fn inversion_9_start() {
    run_case(&inversion_case("inversion_9_start", 17, 0, 9));
}

/// An exact inversion of length 9, one below the minimum template switch length, in the middle of the sequences.
#[test]
fn inversion_9_middle() {
    run_case(&inversion_case("inversion_9_middle", 18, 15, 9));
}

/// An exact inversion of length 9, one below the minimum template switch length, at the end of the sequences.
#[test]
fn inversion_9_end() {
    run_case(&inversion_case("inversion_9_end", 19, 31, 9));
}

/// An inversion of length 12 whose inner part has one mismatch in the middle (position 6).
/// The longest exact run inside is 6, so for k in 7..=12 the template switch needs an inexact anchor.
/// It is supported by k = 12 and k = 10 with one mutation (single anchor),
/// and by k = 6 with one mutation (two adjacent anchors).
#[test]
fn inversion_with_one_mismatch() {
    let reference = random_dna(20, 40);
    let mut query = invert(&reference, 14, 12);
    substitute(&mut query, 14 + 6);
    let case = Case::new("inversion_with_one_mismatch", &reference, &query).with_grid([
        (12, 0),
        (9, 0),
        (6, 0),
        (12, 1),
        (10, 1),
        (6, 1),
        (13, 2),
    ]);
    run_case(&case);
    assert_finds_template_switch(&case, (12, 1));
    assert_finds_template_switch(&case, (10, 1));
    assert_finds_template_switch(&case, (6, 1));
}

/// An inversion of length 10 (the minimum) with one mismatch in the middle.
/// It is supported by k = 10 with one mutation, but (unless the flanks happen to match) neither by k = 10
/// without mutations nor by k = 6 with one mutation (two non-overlapping 6-mers do not fit into 10 characters).
/// So here a smaller k may legitimately lead to a higher cost.
#[test]
fn inversion_10_with_middle_mismatch() {
    let reference = random_dna(21, 40);
    let mut query = invert(&reference, 15, 10);
    substitute(&mut query, 15 + 5);
    let case = Case::new("inversion_10_with_middle_mismatch", &reference, &query).with_grid([
        (10, 1),
        (6, 1),
        (10, 0),
    ]);
    run_case(&case);
    assert_finds_template_switch(&case, (10, 1));
}

/// An inversion of length 12 whose mismatches are at its first and last character,
/// i.e. right at the boundaries of the anchors that cover the inner part.
#[test]
fn inversion_with_boundary_mismatches() {
    let reference = random_dna(22, 40);
    let mut query = invert(&reference, 14, 12);
    substitute(&mut query, 14);
    substitute(&mut query, 14 + 11);
    run_case(&Case::new(
        "inversion_with_boundary_mismatches",
        &reference,
        &query,
    ));
}

/// Two separate exact inversions of length 12 in an 80 bp sequence.
#[test]
fn two_template_switches() {
    let reference = random_dna(23, 80);
    let query = invert(&invert(&reference, 15, 12), 50, 12);
    let case = Case::new("two_template_switches", &reference, &query);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
    let result = chain_result(&case, (6, 0));
    assert_eq!(
        result.template_switch_amount, 2,
        "expected two template switches, got alignment {}",
        result.alignment
    );
}

/// An exact inversion of length 12 with a 1 bp deletion directly before and a 2 bp insertion directly after it.
#[test]
fn template_switch_with_flanking_indels() {
    let reference = random_dna(24, 44);
    let inverted = invert(&reference, 16, 12);
    let query = concat(&[&inverted[..15], &inverted[16..28], b"TT", &inverted[28..]]);
    let case = Case::new("template_switch_with_flanking_indels", &reference, &query);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

/// Substitutions exactly at multiples of 6 and 9, i.e. at the boundaries of tiling k-mers.
#[test]
fn substitutions_at_kmer_boundaries() {
    let reference = random_dna(25, 45);
    let mut query = reference.clone();
    for position in [6, 9, 18, 27, 36] {
        substitute(&mut query, position);
    }
    run_case(&Case::new(
        "substitutions_at_kmer_boundaries",
        &reference,
        &query,
    ));
}

/// Very permissive anchors with max_anchor_mutations = k - 1 on a tiny case with an inversion of length 10
/// that contains one mismatch.
#[test]
fn max_anchor_mutations_k_minus_one() {
    let reference = random_dna(26, 16);
    let mut query = invert(&reference, 3, 10);
    substitute(&mut query, 3 + 4);
    run_case(
        &Case::new("max_anchor_mutations_k_minus_one", &reference, &query)
            .with_grid([(3, 2), (4, 3)]),
    );
}

/// A 30 bp deletion in the middle of a 70 bp sequence.
#[test]
fn large_gap() {
    let reference = random_dna(27, 70);
    let query = concat(&[&reference[..20], &reference[50..]]);
    run_case(&Case::new("large_gap", &reference, &query));
}

/// A query segment of length 12 replaced by the reverse complement of a different, nearby segment of the reference
/// (an inter-sequence template switch that is not an in-place inversion).
#[test]
fn inter_sequence_template_switch() {
    let reference = random_dna(28, 50);
    let mut query = reference.clone();
    query[28..40].copy_from_slice(&reverse_complement(&reference[8..20]));
    let case = Case::new("inter_sequence_template_switch", &reference, &query);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}

/// A query segment of length 12 replaced by the reverse complement of a different, nearby segment of the query
/// (an intra-sequence template switch).
/// The source segment carries two substitutions relative to the reference,
/// such that copying from the query is cheaper than copying from the reference.
#[test]
fn intra_sequence_template_switch() {
    let reference = random_dna(29, 50);
    let mut query = reference.clone();
    substitute(&mut query, 11);
    substitute(&mut query, 16);
    let source = reverse_complement(&query[8..20]);
    query[28..40].copy_from_slice(&source);
    let case = Case::new("intra_sequence_template_switch", &reference, &query);
    run_case(&case);
    assert_finds_template_switch(&case, (6, 0));
}
