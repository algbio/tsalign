use std::iter;

use extend_map::ExtendMap;
use generic_a_star::cost::AStarCost;
use integer_composition::WeakKCompositionIterator;
use itertools::Itertools;

use crate::{
    anchors::kmers::{Kmer, KmerStore},
    costs::GapAffineCosts,
};

#[cfg(test)]
mod tests;

/// Compute the k-mers of the given `sequence`, and all l-mers that can be created by mutating these k-mers up to `max_mutations` times.
pub fn compute_inexact_kmers<Store: KmerStore, Cost: AStarCost>(
    sequence: &[u8],
    k: usize,
    max_mutations: usize,
    costs: &GapAffineCosts<Cost>,
) -> Vec<Vec<(Kmer<Store>, usize, Cost)>> {
    let mut kmers = vec![Vec::new(); 2 * max_mutations + 1];
    let mut kmer_buffer = Vec::new();

    for (kmer, offset) in (0..sequence.len().saturating_sub(k) + 1)
        .map(|offset| (Kmer::<Store>::from(&sequence[offset..offset + k]), offset))
    {
        for mutations in 0..=max_mutations {
            for (i, target_length) in
                ((max_mutations - mutations) * 2..).zip(k.saturating_sub(mutations)..=k + mutations)
            {
                let insertion_count = target_length.saturating_sub(k);
                let deletion_count = k.saturating_sub(target_length);
                let substitution_count = mutations.saturating_sub(insertion_count + deletion_count);
                debug_assert_eq!(
                    insertion_count + deletion_count + substitution_count,
                    mutations
                );
                debug_assert!(insertion_count == 0 || deletion_count == 0);

                if deletion_count > 0 {
                    generate_kmer_deletions(kmer, k, deletion_count, costs, &mut kmer_buffer);
                    for (kmer, deletion_cost) in kmer_buffer.drain(..) {
                        generate_kmer_substitutions(
                            kmer,
                            k,
                            substitution_count,
                            costs,
                            &mut ExtendMap::new(&mut kmers[i], |(kmer, substitution_cost)| {
                                (kmer, offset, deletion_cost + substitution_cost)
                            }),
                        );
                    }
                } else if insertion_count > 0 {
                    generate_kmer_substitutions(
                        kmer,
                        k,
                        substitution_count,
                        costs,
                        &mut kmer_buffer,
                    );
                    for (kmer, substitution_cost) in kmer_buffer.drain(..) {
                        generate_kmer_insertions(
                            kmer,
                            k,
                            insertion_count,
                            costs,
                            &mut ExtendMap::new(&mut kmers[i], |(kmer, insertion_cost)| {
                                (kmer, offset, substitution_cost + insertion_cost)
                            }),
                        );
                    }
                } else {
                    generate_kmer_substitutions(
                        kmer,
                        k,
                        substitution_count,
                        costs,
                        &mut ExtendMap::new(&mut kmers[i], |(kmer, substitution_cost)| {
                            (kmer, offset, substitution_cost)
                        }),
                    );
                }

                kmer_buffer.clear();
            }
        }
    }

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
