//! Utilities for testing algorithms that use alignment histories.

use std::collections::BTreeSet;

use crate::{
    alignment::{Alignment, AlignmentType, GapType},
    alignment_history::{
        history_alignment_operations::AlignmentHistoryOperation, history_vec::AlignmentHistory,
    },
    costs::GapAffineCosts,
};

pub const OPERATIONS: [AlignmentHistoryOperation; 4] = [
    AlignmentHistoryOperation::Match,
    AlignmentHistoryOperation::Substitution,
    AlignmentHistoryOperation::GapInA,
    AlignmentHistoryOperation::GapInB,
];

/// A small deterministic pseudo-random number generator (xorshift64*).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut rng = Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        rng.next_u64();
        rng
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Returns a uniformly random number in `0..n`.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % u64::try_from(n).unwrap()) as usize
    }

    /// Returns a random DNA sequence of the given length.
    pub fn dna(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[self.below(4)]).collect()
    }

    /// Returns a copy of `sequence` where each character is substituted, deleted or followed by an insertion with probability `mutation_percent / 100`.
    pub fn mutate(&mut self, sequence: &[u8], mutation_percent: usize) -> Vec<u8> {
        let mut result = Vec::new();
        for &c in sequence {
            if self.below(100) >= mutation_percent {
                result.push(c);
                continue;
            }

            match self.below(3) {
                0 => result.push(b"ACGT"[self.below(4)]),
                1 => {}
                _ => {
                    result.push(c);
                    result.push(b"ACGT"[self.below(4)]);
                }
            }
        }
        result
    }
}

/// Returns true if the alignment (given from oldest to newest operation) contains a subalignment that could be an anchor.
///
/// Such a subalignment contains at most `max_anchor_mutations` mutations and covers exactly `anchor_k` characters of at least one of the sequences.
/// This is a brute force check over all subalignments.
pub fn contains_anchor(
    alignment: &[AlignmentHistoryOperation],
    anchor_k: usize,
    max_anchor_mutations: usize,
) -> bool {
    (0..alignment.len()).any(|start| {
        let mut len_a = 0;
        let mut len_b = 0;
        let mut mutations = 0;
        alignment[start..].iter().any(|operation| {
            len_a += operation.len_a();
            len_b += operation.len_b();
            mutations += operation.mutations();
            mutations <= max_anchor_mutations && (len_a == anchor_k || len_b == anchor_k)
        })
    })
}

/// Returns true if the alignment (given from oldest to newest operation) has a suffix that could be an anchor.
///
/// If the alignment without its last operation contains no anchor, then this is equivalent to [`contains_anchor`].
pub fn ends_with_anchor(
    alignment: &[AlignmentHistoryOperation],
    anchor_k: usize,
    max_anchor_mutations: usize,
) -> bool {
    let mut len_a = 0;
    let mut len_b = 0;
    let mut mutations = 0;
    alignment.iter().rev().any(|operation| {
        len_a += operation.len_a();
        len_b += operation.len_b();
        mutations += operation.mutations();
        mutations <= max_anchor_mutations && (len_a == anchor_k || len_b == anchor_k)
    })
}

/// Returns the longest suffix of the alignment (given from oldest to newest operation) with at most `max_anchor_mutations` mutations.
///
/// The suffix is returned from newest to oldest operation.
pub fn longest_suffix_newest_first(
    alignment: &[AlignmentHistoryOperation],
    max_anchor_mutations: usize,
) -> Vec<AlignmentHistoryOperation> {
    let mut mutations = 0;
    alignment
        .iter()
        .rev()
        .copied()
        .take_while(|operation| {
            mutations += operation.mutations();
            mutations <= max_anchor_mutations
        })
        .collect()
}

/// Returns all alignment histories that are reachable from the empty history via [`AlignmentHistory::try_extend`].
pub fn reachable_histories<AlignmentHistoryVec: AlignmentHistory>(
    anchor_k: u8,
    max_anchor_mutations: u8,
) -> Vec<AlignmentHistoryVec> {
    let mut reached = BTreeSet::new();
    let mut stack = vec![AlignmentHistoryVec::default()];
    reached.insert(AlignmentHistoryVec::default());

    while let Some(history) = stack.pop() {
        for operation in OPERATIONS {
            if let Some(extension) = history.try_extend(operation, anchor_k, max_anchor_mutations) {
                if reached.insert(extension) {
                    stack.push(extension);
                }
            }
        }
    }

    reached.into_iter().collect()
}

/// The gap-affine cost of the alignment (given from oldest to newest operation), assuming that it does not start within a gap.
pub fn alignment_cost(alignment: &[AlignmentHistoryOperation], costs: &GapAffineCosts<u32>) -> u32 {
    let mut gap_type = GapType::None;
    alignment
        .iter()
        .map(|operation| {
            let (cost, next_gap_type) = match operation {
                AlignmentHistoryOperation::Match => (0, GapType::None),
                AlignmentHistoryOperation::Substitution => (costs.substitution, GapType::None),
                AlignmentHistoryOperation::GapInA => (
                    if gap_type == GapType::InA {
                        costs.gap_extend
                    } else {
                        costs.gap_open
                    },
                    GapType::InA,
                ),
                AlignmentHistoryOperation::GapInB => (
                    if gap_type == GapType::InB {
                        costs.gap_extend
                    } else {
                        costs.gap_open
                    },
                    GapType::InB,
                ),
            };
            gap_type = next_gap_type;
            cost
        })
        .sum()
}

/// Splits the alignment at template switch boundaries and converts each segment into alignment history operations.
pub fn alignment_segments(alignment: &Alignment) -> Vec<Vec<AlignmentHistoryOperation>> {
    let mut segments = vec![Vec::new()];
    for alignment_type in alignment.iter_unpacked() {
        match alignment_type {
            AlignmentType::Match => segments
                .last_mut()
                .unwrap()
                .push(AlignmentHistoryOperation::Match),
            AlignmentType::Substitution => segments
                .last_mut()
                .unwrap()
                .push(AlignmentHistoryOperation::Substitution),
            AlignmentType::GapA => segments
                .last_mut()
                .unwrap()
                .push(AlignmentHistoryOperation::GapInA),
            AlignmentType::GapB => segments
                .last_mut()
                .unwrap()
                .push(AlignmentHistoryOperation::GapInB),
            AlignmentType::TsStart { .. } | AlignmentType::TsEnd { .. } => {
                segments.push(Vec::new())
            }
        }
    }
    segments
}
