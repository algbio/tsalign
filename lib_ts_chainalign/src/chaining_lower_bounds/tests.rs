//! Tests for the lower bounds of gaps that touch the start or end of the alignment.

use generic_a_star::cost::U32Cost;
use num_traits::bounds::UpperBounded;

use crate::{
    alignment::{coordinates::PrimaryAlignmentCoordinates, sequences::AlignmentSequences},
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    max_match_run_chaining::gap_affine::GapAffineAligner,
    panic_on_extend::PanicOnExtend,
};

const MAX_N: usize = 16;

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => panic!("Unsupported character: {c}"),
    }
}

fn costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    )
}

fn exact_lower_bounds(k: u32) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(
        MAX_N,
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        costs(),
    )
}

/// A small deterministic random generator (SplitMix64).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn below(&mut self, limit: usize) -> usize {
        (self.next() % limit as u64) as usize
    }

    /// A random sequence over a small alphabet, such that it contains many matches.
    fn sequence(&mut self, length: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..length)
            .map(|_| alphabet[self.below(alphabet.len())])
            .collect()
    }

    /// A copy of `sequence` with a few random substitutions, insertions and deletions.
    fn mutate(&mut self, sequence: &[u8], alphabet: &[u8]) -> Vec<u8> {
        let mut result = Vec::new();
        for &c in sequence {
            match self.below(8) {
                0 => result.push(alphabet[self.below(alphabet.len())]),
                1 => {}
                2 => {
                    result.push(c);
                    result.push(alphabet[self.below(alphabet.len())]);
                }
                _ => result.push(c),
            }
        }
        result
    }
}

/// Pairs of short sequences with many matches at their starts and ends.
fn test_sequence_pairs() -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = vec![
        (b"CCCTAGTT".to_vec(), b"CCCAAGTT".to_vec()),
        (
            b"A".to_vec(),
            b"TTTTTTTTTATTTTTTTTTT".to_vec()[..12].to_vec(),
        ),
        (b"ACGTACGTAC".to_vec(), b"ACGTACGTAC".to_vec()),
        (b"ACGTTACGTA".to_vec(), b"ACGTACGTA".to_vec()),
        (b"AAAAAAAAAAAA".to_vec(), b"AAAAAAAAA".to_vec()),
        (Vec::new(), b"ACGT".to_vec()),
    ];

    let mut rng = Rng(0x5eed);
    for _ in 0..40 {
        let alphabet: &[u8] = if rng.below(2) == 0 { b"AC" } else { b"ACGT" };
        let length = rng.below(MAX_N - 3) + 1;
        let sequence_a = rng.sequence(length, alphabet);
        let sequence_b = rng.mutate(&sequence_a, alphabet);
        if sequence_b.len() <= MAX_N {
            pairs.push((sequence_a, sequence_b));
        }
    }
    pairs
}

/// For exact anchors, the lower bounds for gaps from the start, to the end and from start to end
/// must not be larger than the costs computed by the aligner that fills these gaps during chaining.
///
/// The lower bound tables for gaps between anchors assume that a gap starts and ends with a non-match,
/// which does not hold at the start and end of the alignment.
#[test]
fn start_and_end_lower_bounds_are_admissible() {
    let costs = costs();
    let mut table_bound_violated = false;

    for k in 2..=10 {
        let lower_bounds = exact_lower_bounds(k);

        for (sequence_a, sequence_b) in test_sequence_pairs() {
            let sequences =
                AlignmentSequences::new_complete(sequence_a.clone(), sequence_b.clone());
            let start = sequences.primary_start();
            let end = sequences.primary_end();
            let mut aligner =
                GapAffineAligner::new(&sequences, &costs.primary_costs, &rc_fn, k - 1, 0);
            let context = || {
                format!(
                    "k = {k}, a = {}, b = {}",
                    String::from_utf8_lossy(&sequence_a),
                    String::from_utf8_lossy(&sequence_b)
                )
            };

            let (cost, _) = aligner.align(start, end, &mut Vec::new(), &mut PanicOnExtend);
            let bound = lower_bounds.primary_lower_bound_start_to_end(end.a(), end.b());
            assert!(
                bound <= cost,
                "Start to end: bound {bound} > cost {cost} for {}",
                context()
            );
            table_bound_violated |= lower_bounds.primary_lower_bound(end.a(), end.b()) > cost;

            for a in 0..=end.a() {
                for b in 0..=end.b() {
                    let point = PrimaryAlignmentCoordinates::new(a, b);

                    if point != end {
                        let (cost, _) =
                            aligner.align(start, point, &mut Vec::new(), &mut PanicOnExtend);
                        if cost < U32Cost::max_value() {
                            let bound = lower_bounds.primary_lower_bound_from_start(a, b);
                            assert!(
                                bound <= cost,
                                "Start to ({a}, {b}): bound {bound} > cost {cost} for {}",
                                context()
                            );
                            table_bound_violated |= lower_bounds.primary_lower_bound(a, b) > cost;
                        }
                    }

                    if point != start {
                        let (cost, _) =
                            aligner.align(point, end, &mut Vec::new(), &mut PanicOnExtend);
                        if cost < U32Cost::max_value() {
                            let (gap1, gap2) = (end.a() - a, end.b() - b);
                            let bound = lower_bounds.primary_lower_bound_to_end(gap1, gap2);
                            assert!(
                                bound <= cost,
                                "({a}, {b}) to end: bound {bound} > cost {cost} for {}",
                                context()
                            );
                            table_bound_violated |=
                                lower_bounds.primary_lower_bound(gap1, gap2) > cost;
                        }
                    }
                }
            }
        }
    }

    // Make sure that the test cases actually exercise the difference to the bounds for gaps between anchors.
    assert!(table_bound_violated);
}

/// The bounds for gaps that touch the start or end are never larger than the bounds for gaps between anchors,
/// and for inexact anchors they are the same.
#[test]
fn start_and_end_lower_bounds_are_at_most_the_inner_bounds() {
    for (k, max_anchor_mutations) in [(4, 0), (8, 0), (4, 1), (8, 2)] {
        let lower_bounds = ChainingLowerBounds::new(
            MAX_N,
            ChainingLowerBoundConfig {
                anchor_k: k,
                max_anchor_mutations,
                inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
            },
            costs(),
        );

        for gap1 in 0..=MAX_N {
            for gap2 in 0..=MAX_N {
                let inner = lower_bounds.primary_lower_bound(gap1, gap2);
                let boundary = [
                    lower_bounds.primary_lower_bound_from_start(gap1, gap2),
                    lower_bounds.primary_lower_bound_to_end(gap1, gap2),
                    lower_bounds.primary_lower_bound_start_to_end(gap1, gap2),
                ];
                for bound in boundary {
                    if max_anchor_mutations == 0 {
                        assert!(bound <= inner);
                    } else {
                        assert_eq!(bound, inner);
                    }
                }
            }

            let gap = gap1;
            let jump_12 = lower_bounds.jump_12_lower_bound_from_start(gap);
            let jump_34 = lower_bounds.jump_34_lower_bound_to_end(gap);
            if max_anchor_mutations == 0 {
                assert!(jump_12 <= lower_bounds.jump_12_lower_bound(gap));
                assert!(jump_34 <= lower_bounds.jump_34_lower_bound(gap));
            } else {
                assert_eq!(jump_12, lower_bounds.jump_12_lower_bound(gap));
                assert_eq!(jump_34, lower_bounds.jump_34_lower_bound(gap));
            }
        }
    }
}
