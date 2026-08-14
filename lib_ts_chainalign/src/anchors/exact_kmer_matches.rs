use std::mem;

use itertools::iproduct;
use num_traits::bounds::LowerBounded;

use crate::anchors::kmers::{Kmer, KmerStore};

struct Cluster<Store> {
    indexes_a: Vec<usize>,
    indexes_b: Vec<usize>,
    kmer: Kmer<Store>,
}

impl<Store> Cluster<Store> {
    fn new(kmer: Kmer<Store>) -> Self {
        Self {
            indexes_a: Vec::new(),
            indexes_b: Vec::new(),
            kmer,
        }
    }

    fn reset(&mut self, kmer: Kmer<Store>) -> impl Iterator<Item = (usize, usize)> {
        self.kmer = kmer;
        let indexes_a = mem::take(&mut self.indexes_a);
        let indexes_b = mem::take(&mut self.indexes_b);

        iproduct!(indexes_a, indexes_b)
    }
}

pub fn find_exact_kmer_matches<Store: KmerStore>(
    mut a: &[(Kmer<Store>, usize)],
    mut b: &[(Kmer<Store>, usize)],
) -> Vec<(usize, usize)> {
    debug_assert!(a.is_sorted());
    debug_assert!(b.is_sorted());

    let mut result = Vec::new();
    let mut cluster = Cluster::new(Kmer::min_value());

    while let (Some((kmer_a, index_a)), Some((kmer_b, index_b))) = (a.first(), b.first()) {
        if kmer_a < kmer_b {
            a = &a[1..];

            if kmer_a != &cluster.kmer {
                result.extend(cluster.reset(*kmer_a));
            }
            cluster.indexes_a.push(*index_a);
        } else {
            b = &b[1..];

            if kmer_b != &cluster.kmer {
                result.extend(cluster.reset(*kmer_b));
            }
            cluster.indexes_b.push(*index_b);
        }
    }

    for (kmer, index) in a {
        if kmer != &cluster.kmer {
            result.extend(cluster.reset(*kmer));
        }
        cluster.indexes_a.push(*index);
    }

    for (kmer, index) in b {
        if kmer != &cluster.kmer {
            result.extend(cluster.reset(*kmer));
        }
        cluster.indexes_b.push(*index);
    }

    result.extend(cluster.reset(cluster.kmer));
    result
}
