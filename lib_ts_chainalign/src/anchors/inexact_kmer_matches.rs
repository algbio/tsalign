use std::{iter, mem, ops::RangeBounds};

use extend_map::ExtendMap;
use generic_a_star::cost::AStarCost;
use integer_composition::WeakKCompositionIterator;
use itertools::{Itertools, iproduct};
use num_traits::bounds::LowerBounded;

use crate::{
    anchors::{
        exact_kmer_matches::compute_exact_kmers,
        kmers::{Kmer, KmerStore},
    },
    costs::GapAffineCosts,
};

#[cfg(test)]
mod tests;

struct Cluster<Store, Cost> {
    indexes_a: Vec<(usize, Cost)>,
    indexes_b: Vec<usize>,
    kmer: Kmer<Store>,
}

impl<Store, Cost> Cluster<Store, Cost> {
    /// Creates a new cluster with the given k-mer.
    fn new(kmer: Kmer<Store>) -> Self {
        Self {
            indexes_a: Vec::new(),
            indexes_b: Vec::new(),
            kmer,
        }
    }

    /// Resets the cluster with a new k-mer and returns an iterator over the previous matches.
    fn reset(&mut self, kmer: Kmer<Store>) -> impl Iterator<Item = (usize, usize, Cost)>
    where
        Cost: Clone,
    {
        self.kmer = kmer;
        let indexes_a = mem::take(&mut self.indexes_a);
        let indexes_b = mem::take(&mut self.indexes_b);

        iproduct!(indexes_a, indexes_b).map(|((a, c), b)| (a, b, c))
    }

    /// Returns the current k-mer of the cluster.
    fn kmer(&self) -> &Kmer<Store> {
        &self.kmer
    }

    /// Pushes a k-mer incidence from sequence A into the cluster.
    fn push_a(&mut self, index: usize, cost: Cost) {
        self.indexes_a.push((index, cost));
    }

    /// Pushes a k-mer incidence from sequence B into the cluster.
    fn push_b(&mut self, index: usize) {
        self.indexes_b.push(index);
    }
}

/// Compute the k-mers of the given `sequence`, and all l-mers that can be created by mutating these k-mers up to `max_mutations` times.
///
/// * `applied_sequence_offset` is the offset to apply to the k-mer positions.
///   If the sequence is for example a subsequence starting at position `i`, then `applied_sequence_offset` should be set to `i`.
pub fn compute_inexact_kmers<Store: KmerStore, Cost: AStarCost>(
    sequence: &[u8],
    applied_sequence_offset: usize,
    k: usize,
    max_mutations: usize,
    costs: &GapAffineCosts<Cost>,
) -> Vec<Vec<(Kmer<Store>, usize, Cost)>> {
    let mut kmers = vec![Vec::new(); 2 * max_mutations + 1];
    let mut deleted_buffer = Vec::new();
    let mut del_sub_buffer = Vec::new();

    for offset in 0..sequence.len().saturating_sub(k) + 1 {
        let kmer = Kmer::<Store>::from(&sequence[offset..offset + k]);

        for mutations in 0..=max_mutations {
            for (i, target_length) in
                ((max_mutations - mutations)..).zip(k.saturating_sub(mutations)..=k + mutations)
            {
                // Apply insertions or deletions to reach the target length.
                let insertion_count = target_length.saturating_sub(k);
                let deletion_count = k.saturating_sub(target_length);

                // The remainder of the mutations can be filled with pairs of insertions and deletions, or substitutions.
                let remaining_mutations = mutations - insertion_count - deletion_count;

                for substitution_count in (0..=remaining_mutations).rev().step_by(2) {
                    let remaining_mutations = remaining_mutations - substitution_count;
                    debug_assert_eq!(remaining_mutations % 2, 0);
                    let insertion_count = insertion_count + remaining_mutations / 2;
                    let deletion_count = deletion_count + remaining_mutations / 2;

                    debug_assert_eq!(
                        insertion_count + deletion_count + substitution_count,
                        mutations
                    );

                    generate_kmer_deletions(kmer, k, deletion_count, costs, &mut deleted_buffer);
                    for (kmer, deletion_cost) in deleted_buffer.drain(..) {
                        generate_kmer_substitutions(
                            kmer,
                            k - deletion_count,
                            substitution_count,
                            costs,
                            &mut ExtendMap::new(
                                &mut del_sub_buffer,
                                |(kmer, substitution_cost)| {
                                    (kmer, deletion_cost + substitution_cost)
                                },
                            ),
                        );
                    }
                    for (kmer, del_sub_cost) in del_sub_buffer.drain(..) {
                        generate_kmer_insertions(
                            kmer,
                            k - deletion_count,
                            insertion_count,
                            costs,
                            &mut ExtendMap::new(&mut kmers[i], |(kmer, insertion_cost)| {
                                (
                                    kmer,
                                    applied_sequence_offset + offset,
                                    del_sub_cost + insertion_cost,
                                )
                            }),
                        );
                    }
                }
            }
        }
    }

    // Sort and deduplicate k-mers, and filter suboptimal k-mers.
    for kmers in &mut kmers {
        kmers.sort_unstable();
        let mut previous = (Kmer::<Store>::default(), usize::MAX, Cost::max_value());
        kmers.retain(|(kmer, offset, cost)| {
            if *kmer == previous.0 && *offset == previous.1 {
                debug_assert!(previous.2 <= *cost);
                false
            } else {
                previous = (*kmer, *offset, *cost);
                true
            }
        });
    }

    kmers
}

fn generate_kmer_insertions<Store: KmerStore, Cost: AStarCost>(
    kmer: Kmer<Store>,
    k: usize,
    insertion_count: usize,
    costs: &GapAffineCosts<Cost>,
    output: &mut impl Extend<(Kmer<Store>, Cost)>,
) {
    if insertion_count == 0 {
        output.extend(iter::once((kmer, Cost::zero())));
        return;
    }

    let kmer = kmer.to_vec(k);

    // Iterate over all possible combinations of insertion locations and lengths.
    let mut insertion_pattern_iterator =
        WeakKCompositionIterator::new(k + 1, insertion_count).unwrap();
    while let Some(insertion_pattern) = insertion_pattern_iterator.next_borrowing() {
        let cost: Cost = insertion_pattern
            .iter()
            .map(|insertion_size| costs.gap_cost(*insertion_size))
            .sum();

        // Iterate over all possible combinations of insertion sequences.
        for insertion_sequence in (0..insertion_count)
            .map(|_| b"ACGT".iter().copied())
            .multi_cartesian_product()
        {
            let mut insertion_sequence = insertion_sequence.into_iter();
            let mut result_kmer = Kmer::<Store>::default();
            for (i, c) in kmer.iter().copied().enumerate() {
                for _ in 0..insertion_pattern[i] {
                    result_kmer.push(insertion_sequence.next().unwrap());
                }
                result_kmer.push(c);
            }
            for _ in 0..insertion_pattern[k] {
                result_kmer.push(insertion_sequence.next().unwrap());
            }

            output.extend(iter::once((result_kmer, cost)));
        }
    }
}

fn generate_kmer_deletions<Store: KmerStore, Cost: AStarCost>(
    kmer: Kmer<Store>,
    k: usize,
    deletion_count: usize,
    costs: &GapAffineCosts<Cost>,
    output: &mut impl Extend<(Kmer<Store>, Cost)>,
) {
    if deletion_count == 0 {
        output.extend(iter::once((kmer, Cost::zero())));
        return;
    }

    let kmer = kmer.to_vec(k);

    for mut deletion_pattern in (0..k).rev().combinations(deletion_count) {
        let cost = deletion_pattern
            .iter()
            .rev()
            .fold(
                (usize::MAX - 1, Cost::zero()),
                |(last_deletion, total_cost), &deletion_index| {
                    if last_deletion + 1 == deletion_index {
                        (deletion_index, total_cost + costs.gap_extend)
                    } else {
                        (deletion_index, total_cost + costs.gap_open)
                    }
                },
            )
            .1;
        let result_kmer = kmer
            .iter()
            .enumerate()
            .filter_map(|(i, c)| {
                if Some(i) == deletion_pattern.last().copied() {
                    deletion_pattern.pop();
                    None
                } else {
                    Some(*c)
                }
            })
            .collect();
        output.extend(iter::once((result_kmer, cost)));
    }
}

/// Rotates a byte DNA character by n positions in alphabetic order.
fn rotate_byte(b: u8, n: u8) -> u8 {
    debug_assert!(n < 4);
    match n {
        0 => b,
        1 => match b {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            _ => unreachable!(),
        },
        2 => match b {
            b'A' => b'G',
            b'C' => b'T',
            b'G' => b'A',
            b'T' => b'C',
            _ => unreachable!(),
        },
        3 => match b {
            b'A' => b'T',
            b'C' => b'A',
            b'G' => b'C',
            b'T' => b'G',
            _ => unreachable!(),
        },
        _ => unreachable!(),
    }
}

fn generate_kmer_substitutions<Store: KmerStore, Cost: AStarCost>(
    kmer: Kmer<Store>,
    k: usize,
    substitution_count: usize,
    costs: &GapAffineCosts<Cost>,
    output: &mut impl Extend<(Kmer<Store>, Cost)>,
) {
    if substitution_count == 0 {
        output.extend(iter::once((kmer, Cost::zero())));
        return;
    }

    let kmer = kmer.to_vec(k);
    let cost = costs.substitution * Cost::from_usize(substitution_count);

    for substitution_sequence in (0..substitution_count)
        .map(|_| 1..4)
        .multi_cartesian_product()
    {
        for mut substitution_pattern in (0..k).rev().combinations(substitution_count) {
            let result_kmer = kmer
                .iter()
                .copied()
                .enumerate()
                .map(|(i, c)| {
                    if Some(i) == substitution_pattern.last().copied() {
                        substitution_pattern.pop();
                        rotate_byte(c, substitution_sequence[substitution_pattern.len()])
                    } else {
                        c
                    }
                })
                .collect();
            output.extend(iter::once((result_kmer, cost)));
        }
    }
}

/// For each k in the given range, computes all k-mers of a sequence and returns them in alphabetical order.
/// The result is a [`Vec`] of `Vec`s, with one inner `Vec` for each k.
///
/// * `applied_sequence_offset` is the offset to apply to the k-mer positions.
///   If the sequence is for example a subsequence starting at position `i`, then `applied_sequence_offset` should be set to `i`.
pub fn compute_exact_kmers_range<Store: KmerStore>(
    sequence: &[u8],
    applied_sequence_offset: usize,
    k_range: &impl RangeBounds<usize>,
) -> Vec<Vec<(Kmer<Store>, usize)>> {
    let start_inclusive = match k_range.start_bound() {
        std::ops::Bound::Included(start) => *start,
        std::ops::Bound::Excluded(start) => start + 1,
        std::ops::Bound::Unbounded => 0,
    };
    let end_exclusive = match k_range.end_bound() {
        std::ops::Bound::Included(end) => end + 1,
        std::ops::Bound::Excluded(end) => *end,
        std::ops::Bound::Unbounded => sequence.len() + 1,
    };

    (start_inclusive..end_exclusive)
        .map(|k| compute_exact_kmers::<Store>(sequence, applied_sequence_offset, k))
        .collect()
}

/// Computes the inexact k-mer matches between the two sequences.
///
/// The inexact k-mer matches are exact matches between the inexact k-mers of the first sequence and the exact k-mers of the second sequence.
///
/// # Return value
/// (`first_sequence_offset`, `second_sequence_offset`, `cost`)
pub fn find_inexact_kmer_matches<Store: KmerStore, Cost: AStarCost>(
    mut a: &[(Kmer<Store>, usize, Cost)],
    mut b: &[(Kmer<Store>, usize)],
) -> Vec<(usize, usize, Cost)> {
    debug_assert!(a.is_sorted());
    debug_assert!(b.is_sorted());

    let mut result = Vec::new();
    let mut cluster = Cluster::new(Kmer::min_value());

    while let (Some((kmer_a, index_a, cost)), Some((kmer_b, index_b))) = (a.first(), b.first()) {
        if kmer_a < kmer_b {
            a = &a[1..];

            if kmer_a != cluster.kmer() {
                result.extend(cluster.reset(*kmer_a));
            }
            cluster.push_a(*index_a, *cost);
        } else {
            b = &b[1..];

            if kmer_b != cluster.kmer() {
                result.extend(cluster.reset(*kmer_b));
            }
            cluster.push_b(*index_b);
        }
    }

    for (kmer, index, cost) in a {
        if kmer != cluster.kmer() {
            result.extend(cluster.reset(*kmer));
        }
        cluster.push_a(*index, *cost);
    }

    for (kmer, index) in b {
        if kmer != cluster.kmer() {
            result.extend(cluster.reset(*kmer));
        }
        cluster.push_b(*index);
    }

    result.extend(cluster.reset(*cluster.kmer()));
    result
}
