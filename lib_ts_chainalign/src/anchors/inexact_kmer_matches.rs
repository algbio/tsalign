use generic_a_star::cost::AStarCost;

use crate::{
    anchors::kmers::{Kmer, KmerStore},
    costs::GapAffineCosts,
};

/// Compute the k-mers of the given `sequence`, and all l-mers that can be created by mutating these k-mers up to `max_mutations` times.
pub fn compute_inexact_kmers<Store: KmerStore, Cost: AStarCost>(
    sequence: &[u8],
    k: usize,
    max_mutations: usize,
    costs: &GapAffineCosts<Cost>,
) -> Vec<Vec<(Kmer<Store>, usize, Cost)>> {
    let mut kmers = vec![Vec::new(); 2 * max_mutations + 1];

    for (kmer, offset) in (0..sequence.len().saturating_sub(k) + 1)
        .map(|offset| (Kmer::<Store>::from(&sequence[offset..offset + k]), offset))
    {
        for target_length in k.saturating_sub(max_mutations)..=k + max_mutations {
            let insertion_count = target_length.saturating_sub(k);
            let deletion_count = k.saturating_sub(target_length);
            let substitution_count = max_mutations.saturating_sub(insertion_count + deletion_count);
            debug_assert_eq!(
                insertion_count + deletion_count + substitution_count,
                max_mutations
            );
            debug_assert!(insertion_count == 0 || deletion_count == 0);
        }

        todo!()
    }

    todo!("Sort and remove non-minimal cost kmers");

    kmers
}

fn generate_kmer_insertions<Store: KmerStore, Cost: AStarCost>(
    kmer: Kmer<Store>,
    k: usize,
    insertion_count: usize,
    output: &mut impl Extend<(Kmer<Store>, Cost)>,
) {
    // TODO Use integer-composition crate to generate all possible combinations of insertion positions and sizes.
    //      Then generate all the corresponding possible insertions.
    todo!()
}
