use std::{fmt::Display, time::Instant};

use generic_a_star::cost::AStarCost;
use log::{info, trace};
use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::{
            AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            range::{
                AnySecondaryAlignmentRange, PrimaryAlignmentRange, SpecificSecondaryAlignmentRange,
            },
        },
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    anchors::{
        exact_kmer_matches::{compute_exact_kmers, find_exact_kmer_matches},
        index::AnchorIndex,
        inexact_kmer_matches::{
            compute_exact_kmers_range, compute_inexact_kmers, find_inexact_kmer_matches,
        },
        kmers::KmerStore,
        primary::PrimaryAnchor,
        secondary::SecondaryAnchor,
    },
    costs::GapAffineCosts,
};

pub mod exact_kmer_matches;
pub mod index;
pub mod inexact_kmer_matches;
pub mod kmers;
pub mod primary;
pub mod reverse_lookup;
pub mod secondary;
#[cfg(test)]
mod tests;

#[derive(Debug, PartialEq, Eq)]
pub struct Anchors<Cost> {
    primary: Vec<PrimaryAnchor<Cost>>,
    secondaries: [Vec<SecondaryAnchor<Cost>>; 4],
}

impl<Cost> Anchors<Cost> {
    pub fn new_exact(sequences: &AlignmentSequences, k: u32, rc_fn: &dyn Fn(u8) -> u8) -> Self
    where
        Cost: Zero + Ord,
    {
        if k <= 8 {
            Self::new_exact_with_kmer_store::<u16>(sequences, k, rc_fn)
        } else if k <= 16 {
            Self::new_exact_with_kmer_store::<u32>(sequences, k, rc_fn)
        } else if k <= 32 {
            Self::new_exact_with_kmer_store::<u64>(sequences, k, rc_fn)
        } else if k <= 64 {
            Self::new_exact_with_kmer_store::<u128>(sequences, k, rc_fn)
        } else {
            panic!("Only k-mer sizes up to 64 are supported, but got {k}");
        }
    }

    fn new_exact_with_kmer_store<Store: KmerStore>(
        sequences: &AlignmentSequences,
        k: u32,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self
    where
        Cost: Zero + Ord,
    {
        let start_time = Instant::now();

        let k = usize::try_from(k).unwrap();
        let s1 = sequences.seq1();
        let s2 = sequences.seq2();
        let s1_rc: Vec<_> = s1.iter().copied().rev().map(rc_fn).collect();
        let s2_rc: Vec<_> = s2.iter().copied().rev().map(rc_fn).collect();

        // Compute k-mers.
        let s1_kmers = compute_exact_kmers::<Store>(
            &s1[sequences.primary_start().a()..sequences.primary_end().a()],
            sequences.primary_start().a(),
            k,
        );
        let s2_kmers = compute_exact_kmers::<Store>(
            &s2[sequences.primary_start().b()..sequences.primary_end().b()],
            sequences.primary_start().a(),
            k,
        );
        let s1_rc_kmers = compute_exact_kmers::<Store>(&s1_rc, 0, k);
        let s2_rc_kmers = compute_exact_kmers::<Store>(&s2_rc, 0, k);

        trace!("s1_kmers: {s1_kmers:?}");
        trace!("s2_kmers: {s2_kmers:?}");

        // Compute anchors.
        let mut primary: Vec<_> = find_exact_kmer_matches(&s1_kmers, &s2_kmers)
            .into_iter()
            .map(|(seq1, seq2)| {
                PrimaryAnchor::new_exact(PrimaryAlignmentCoordinates::new(seq1, seq2), k)
            })
            .collect();
        let secondary_11: Vec<_> = find_exact_kmer_matches(&s1_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| {
                SecondaryAnchor::new_exact(
                    AnySecondaryAlignmentCoordinates::new(s1.len() - ancestor, descendant),
                    k,
                )
            })
            .collect();
        let secondary_12: Vec<_> = find_exact_kmer_matches(&s1_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| {
                SecondaryAnchor::new_exact(
                    AnySecondaryAlignmentCoordinates::new(s1.len() - ancestor, descendant),
                    k,
                )
            })
            .collect();
        let secondary_21: Vec<_> = find_exact_kmer_matches(&s2_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| {
                SecondaryAnchor::new_exact(
                    AnySecondaryAlignmentCoordinates::new(s2.len() - ancestor, descendant),
                    k,
                )
            })
            .collect();
        let secondary_22: Vec<_> = find_exact_kmer_matches(&s2_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| {
                SecondaryAnchor::new_exact(
                    AnySecondaryAlignmentCoordinates::new(s2.len() - ancestor, descendant),
                    k,
                )
            })
            .collect();
        let mut secondaries = [secondary_11, secondary_12, secondary_21, secondary_22];

        // Sort anchors.
        primary.sort_unstable();
        for secondary in &mut secondaries {
            secondary.sort_unstable();
        }

        let duration = start_time.elapsed();

        info!(
            "Found {} anchors ({} + {} + {} + {} + {}) in {:.0}ms",
            primary.len() + secondaries.iter().map(Vec::len).sum::<usize>(),
            primary.len(),
            secondaries[0].len(),
            secondaries[1].len(),
            secondaries[2].len(),
            secondaries[3].len(),
            duration.as_secs_f64() * 1e3,
        );

        Self {
            primary,
            secondaries,
        }
    }

    pub fn new_inexact(
        sequences: &AlignmentSequences,
        k: u32,
        max_mutations: u8,
        costs: &GapAffineCosts<Cost>,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self
    where
        Cost: AStarCost,
    {
        if k + u32::from(max_mutations) <= 8 {
            Self::new_inexact_with_kmer_store::<u16>(sequences, k, max_mutations, costs, rc_fn)
        } else if k + u32::from(max_mutations) <= 16 {
            Self::new_inexact_with_kmer_store::<u32>(sequences, k, max_mutations, costs, rc_fn)
        } else if k + u32::from(max_mutations) <= 32 {
            Self::new_inexact_with_kmer_store::<u64>(sequences, k, max_mutations, costs, rc_fn)
        } else if k + u32::from(max_mutations) <= 64 {
            Self::new_inexact_with_kmer_store::<u128>(sequences, k, max_mutations, costs, rc_fn)
        } else {
            panic!(
                "Only k-mer sizes up to 64 are supported, but got k = {k} and max_mutations = {max_mutations}, resulting in a maximum k-mer size of k + max_mutations = {}",
                k + u32::from(max_mutations),
            );
        }
    }

    fn new_inexact_with_kmer_store<Store: KmerStore>(
        sequences: &AlignmentSequences,
        k: u32,
        max_mutations: u8,
        costs: &GapAffineCosts<Cost>,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self
    where
        Cost: AStarCost,
    {
        let start_time = Instant::now();

        let k = usize::try_from(k).unwrap();
        let max_mutations = usize::from(max_mutations);
        let k_range = k.saturating_sub(max_mutations)..=k + max_mutations;
        let s1 = sequences.seq1();
        let s2 = sequences.seq2();
        let s1_rc: Vec<_> = s1.iter().copied().rev().map(rc_fn).collect();
        let s2_rc: Vec<_> = s2.iter().copied().rev().map(rc_fn).collect();
        let s1_rc = s1_rc.as_slice();
        let s2_rc = s2_rc.as_slice();

        // Compute k-mers.
        let s1_inexact_kmers = compute_inexact_kmers::<Store, _>(
            &s1[sequences.primary_start().a()..sequences.primary_end().a()],
            sequences.primary_start().a(),
            k,
            max_mutations,
            costs,
        );
        let s2_inexact_kmers = compute_inexact_kmers::<Store, _>(
            &s2[sequences.primary_start().b()..sequences.primary_end().b()],
            sequences.primary_start().b(),
            k,
            max_mutations,
            costs,
        );
        let s2_exact_kmers: Vec<_> = compute_exact_kmers_range(
            &s2[sequences.primary_start().b()..sequences.primary_end().b()],
            sequences.primary_start().b(),
            &k_range,
        );

        let s1_rc_exact_kmers = compute_exact_kmers_range::<Store>(s1_rc, 0, &k_range);
        let s2_rc_exact_kmers = compute_exact_kmers_range::<Store>(s2_rc, 0, &k_range);

        // Compute anchors.
        let mut primary: Vec<_> = k_range
            .clone()
            .zip(s1_inexact_kmers.iter().zip(s2_exact_kmers.iter()))
            .flat_map(|(local_k, (s1_kmers, s2_kmers))| {
                find_inexact_kmer_matches(s1_kmers, s2_kmers)
                    .into_iter()
                    .map(move |(offset1, offset2, cost)| {
                        PrimaryAnchor::new(
                            PrimaryAlignmentRange::new(
                                PrimaryAlignmentCoordinates::new(offset1, offset2),
                                // Inexact l-mers from sequence 1 are produced by transforming a k-mer and matched against an exact l-mer from sequence 2.
                                // The anchor is with respect to the original sequences, and hence for sequence 1 we take k characters to get the k-mer that was transformed, and for sequence 2 we take local_k characters to get the l-mer that was matched.
                                PrimaryAlignmentCoordinates::new(offset1 + k, offset2 + local_k),
                            ),
                            cost,
                        )
                        .trim(s1, s2)
                    })
            })
            .collect();
        let secondary_11: Vec<_> = k_range
            .clone()
            .zip(s1_inexact_kmers.iter().zip(s1_rc_exact_kmers.iter()))
            .flat_map(|(local_k, (s1_inexact_kmers, s1_rc_exact_kmers))| {
                find_inexact_kmer_matches(s1_inexact_kmers, s1_rc_exact_kmers)
                    .into_iter()
                    .map(move |(offset_descendant, offset_ancestor, cost)| {
                        SecondaryAnchor::new(
                            AnySecondaryAlignmentRange::new(
                                AnySecondaryAlignmentCoordinates::new(
                                    s1.len() - offset_ancestor,
                                    offset_descendant,
                                ),
                                AnySecondaryAlignmentCoordinates::new(
                                    s1.len() - offset_ancestor - k,
                                    offset_descendant + local_k,
                                ),
                            ),
                            cost,
                        )
                        .trim(s1_rc, s1)
                    })
            })
            .collect();
        let secondary_12: Vec<_> = k_range
            .clone()
            .zip(s2_inexact_kmers.iter().zip(s1_rc_exact_kmers.iter()))
            .flat_map(|(local_k, (s2_inexact_kmers, s1_rc_exact_kmers))| {
                find_inexact_kmer_matches(s2_inexact_kmers, s1_rc_exact_kmers)
                    .into_iter()
                    .map(move |(offset_descendant, offset_ancestor, cost)| {
                        SecondaryAnchor::new(
                            AnySecondaryAlignmentRange::new(
                                AnySecondaryAlignmentCoordinates::new(
                                    s1.len() - offset_ancestor,
                                    offset_descendant,
                                ),
                                AnySecondaryAlignmentCoordinates::new(
                                    s1.len() - offset_ancestor - k,
                                    offset_descendant + local_k,
                                ),
                            ),
                            cost,
                        )
                        .trim(s1_rc, s2)
                    })
            })
            .collect();
        let secondary_21: Vec<_> = k_range
            .clone()
            .zip(s1_inexact_kmers.iter().zip(s2_rc_exact_kmers.iter()))
            .flat_map(|(local_k, (s1_inexact_kmers, s2_rc_exact_kmers))| {
                find_inexact_kmer_matches(s1_inexact_kmers, s2_rc_exact_kmers)
                    .into_iter()
                    .map(move |(offset_descendant, offset_ancestor, cost)| {
                        SecondaryAnchor::new(
                            AnySecondaryAlignmentRange::new(
                                AnySecondaryAlignmentCoordinates::new(
                                    s2.len() - offset_ancestor,
                                    offset_descendant,
                                ),
                                AnySecondaryAlignmentCoordinates::new(
                                    s2.len() - offset_ancestor - k,
                                    offset_descendant + local_k,
                                ),
                            ),
                            cost,
                        )
                        .trim(s2_rc, s1)
                    })
            })
            .collect();
        let secondary_22: Vec<_> = k_range
            .clone()
            .zip(s2_inexact_kmers.iter().zip(s2_rc_exact_kmers.iter()))
            .flat_map(|(local_k, (s2_inexact_kmers, s2_rc_exact_kmers))| {
                find_inexact_kmer_matches(s2_inexact_kmers, s2_rc_exact_kmers)
                    .into_iter()
                    .map(move |(offset_descendant, offset_ancestor, cost)| {
                        SecondaryAnchor::new(
                            AnySecondaryAlignmentRange::new(
                                AnySecondaryAlignmentCoordinates::new(
                                    s2.len() - offset_ancestor,
                                    offset_descendant,
                                ),
                                AnySecondaryAlignmentCoordinates::new(
                                    s2.len() - offset_ancestor - k,
                                    offset_descendant + local_k,
                                ),
                            ),
                            cost,
                        )
                        .trim(s2_rc, s2)
                    })
            })
            .collect();

        let mut secondaries = [secondary_11, secondary_12, secondary_21, secondary_22];

        // Sort anchors.
        primary.sort_unstable();
        for secondary in &mut secondaries {
            secondary.sort_unstable();
        }

        // Deduplicate anchors and remove anchors with suboptimal cost.
        let mut previous_primary =
            PrimaryAnchor::new_exact(PrimaryAlignmentCoordinates::new(usize::MAX, usize::MAX), 0);
        primary.retain(|primary| {
            if primary.is_at(previous_primary.range()) {
                debug_assert!(primary.cost() <= previous_primary.cost());
                false
            } else {
                previous_primary = *primary;
                true
            }
        });
        for secondary in &mut secondaries {
            let mut previous_secondary = SecondaryAnchor::new_exact(
                AnySecondaryAlignmentCoordinates::new(usize::MAX, usize::MAX),
                0,
            );
            secondary.retain(|secondary| {
                if secondary.is_at(previous_secondary.range()) {
                    debug_assert!(secondary.cost() <= previous_secondary.cost());
                    false
                } else {
                    previous_secondary = *secondary;
                    true
                }
            });
        }

        let duration = start_time.elapsed();

        info!(
            "Found {} anchors ({} + {} + {} + {} + {}) in {:.0}ms",
            primary.len() + secondaries.iter().map(Vec::len).sum::<usize>(),
            primary.len(),
            secondaries[0].len(),
            secondaries[1].len(),
            secondaries[2].len(),
            secondaries[3].len(),
            duration.as_secs_f64() * 1e3,
        );

        Self {
            primary,
            secondaries,
        }
    }

    pub fn primary(&self, index: AnchorIndex) -> &PrimaryAnchor<Cost> {
        &self.primary[index.as_usize()]
    }

    pub fn primary_len(&self) -> AnchorIndex {
        self.primary.len().into()
    }

    pub fn enumerate_primaries(&self) -> impl Iterator<Item = (AnchorIndex, PrimaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.primary
            .iter()
            .copied()
            .enumerate()
            .map(|(index, anchor)| (index.into(), anchor))
    }

    pub fn primary_index_from_range(&self, range: PrimaryAlignmentRange) -> Option<AnchorIndex>
    where
        Cost: Copy,
    {
        self.enumerate_primaries()
            .filter_map(|(index, anchor)| (anchor.is_at(range)).then_some(index))
            .next()
    }

    fn secondary_anchor_vec(&self, ts_kind: TsKind) -> &Vec<SecondaryAnchor<Cost>> {
        &self.secondaries[ts_kind.index()]
    }

    pub fn secondary(&self, index: AnchorIndex, ts_kind: TsKind) -> &SecondaryAnchor<Cost> {
        &self.secondary_anchor_vec(ts_kind)[index.as_usize()]
    }

    pub fn secondary_len(&self, ts_kind: TsKind) -> AnchorIndex {
        self.secondary_anchor_vec(ts_kind).len().into()
    }

    pub fn enumerate_secondaries(
        &self,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (AnchorIndex, SecondaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.secondary_anchor_vec(ts_kind)
            .iter()
            .copied()
            .enumerate()
            .map(|(index, anchor)| (index.into(), anchor))
    }

    pub fn secondary_index_from_range(
        &self,
        start: SpecificSecondaryAlignmentRange,
    ) -> Option<AnchorIndex>
    where
        Cost: Copy,
    {
        self.enumerate_secondaries(start.ts_kind())
            .filter_map(|(index, anchor)| (anchor.is_at(start.into())).then_some(index))
            .next()
    }
}

impl<Cost: Display> Display for Anchors<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "P: [")?;
        let mut once = true;
        for primary_anchor in &self.primary {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{primary_anchor}")?;
        }
        writeln!(f, "]")?;

        let mut ts_kind_once = true;
        for ts_kind in TsKind::iter() {
            if ts_kind_once {
                ts_kind_once = false;
            } else {
                writeln!(f)?;
            }

            write!(f, "S{}: [", ts_kind.digits())?;
            let mut once = true;
            for secondary_anchor in self.secondary_anchor_vec(ts_kind) {
                if once {
                    once = false;
                } else {
                    write!(f, ", ")?;
                }
                write!(f, "{secondary_anchor}")?;
            }
            write!(f, "]")?;
        }

        Ok(())
    }
}
