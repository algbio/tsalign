use std::{fmt::Display, time::Instant};

use generic_a_star::cost::AStarCost;
use inexact_alignment_anchors::generate_inexact_alignment_anchors_of_subsequences;
use log::{info, trace};
use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::{AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates},
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    anchors::{
        exact_kmer_matches::{compute_exact_kmers, find_exact_kmer_matches},
        index::AnchorIndex,
        kmers::KmerStore,
        primary::PrimaryAnchor,
        secondary::SecondaryAnchor,
    },
    costs::AlignmentCosts,
};

pub mod domination;
pub mod exact_kmer_matches;
pub mod index;
/// Superseded by the crate `inexact_alignment_anchors`, only kept for its tests.
#[cfg(test)]
mod inexact_kmer_matches;
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
            sequences.primary_start().b(),
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
        costs: &AlignmentCosts<Cost>,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self
    where
        Cost: AStarCost,
    {
        let start_time = Instant::now();

        let k = usize::try_from(k).unwrap();
        let max_mutations = usize::from(max_mutations);
        let s1 = sequences.seq1();
        let s2 = sequences.seq2();
        let s1_rc: Vec<_> = s1.iter().copied().rev().map(rc_fn).collect();
        let s2_rc: Vec<_> = s2.iter().copied().rev().map(rc_fn).collect();
        let s1_rc = s1_rc.as_slice();
        let s2_rc = s2_rc.as_slice();

        // Compute anchors.
        let mut primary: Vec<_> = generate_inexact_alignment_anchors_of_subsequences(
            s1,
            s2,
            sequences.range1(),
            sequences.range2(),
            &costs.primary_costs,
            k,
            max_mutations,
        )
        .unwrap()
        .map(PrimaryAnchor::from)
        .collect();
        let secondary_11: Vec<_> = generate_inexact_alignment_anchors_of_subsequences(
            s1_rc,
            s1,
            0..s1_rc.len(),
            sequences.range1(),
            &costs.secondary_costs,
            k,
            max_mutations,
        )
        .unwrap()
        .map(|anchor| SecondaryAnchor::new_from_inexact_alignment_anchor(anchor, s1.len()))
        .collect();
        let secondary_12: Vec<_> = generate_inexact_alignment_anchors_of_subsequences(
            s1_rc,
            s2,
            0..s1_rc.len(),
            sequences.range2(),
            &costs.secondary_costs,
            k,
            max_mutations,
        )
        .unwrap()
        .map(|anchor| SecondaryAnchor::new_from_inexact_alignment_anchor(anchor, s1.len()))
        .collect();
        let secondary_21: Vec<_> = generate_inexact_alignment_anchors_of_subsequences(
            s2_rc,
            s1,
            0..s2_rc.len(),
            sequences.range1(),
            &costs.secondary_costs,
            k,
            max_mutations,
        )
        .unwrap()
        .map(|anchor| SecondaryAnchor::new_from_inexact_alignment_anchor(anchor, s2.len()))
        .collect();
        let secondary_22: Vec<_> = generate_inexact_alignment_anchors_of_subsequences(
            s2_rc,
            s2,
            0..s2_rc.len(),
            sequences.range2(),
            &costs.secondary_costs,
            k,
            max_mutations,
        )
        .unwrap()
        .map(|anchor| SecondaryAnchor::new_from_inexact_alignment_anchor(anchor, s2.len()))
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

    /// Removes the dominated primary anchors, if that is known to not change the cost of the optimal chain.
    ///
    /// See the [`domination`] module for the definition of dominated anchors and the conditions for pruning.
    /// The proof of soundness holds only for the max-match-run chaining,
    /// so this must not be called if the anchors are chained with a different method.
    ///
    /// * `k` and `max_mutations` are the parameters that the anchors were created with by [`Self::new_inexact`].
    ///
    /// Returns the number of removed anchors.
    pub fn prune_dominated_primary_anchors(
        &mut self,
        sequences: &AlignmentSequences,
        k: u32,
        max_mutations: u8,
        costs: &AlignmentCosts<Cost>,
    ) -> usize
    where
        Cost: AStarCost,
    {
        if !domination::is_pruning_sound(usize::from(max_mutations), &costs.primary_costs) {
            return 0;
        }

        let pruned = domination::prune_dominated_primary_anchors(
            &mut self.primary,
            sequences.seq1(),
            sequences.seq2(),
            usize::try_from(k).unwrap(),
            &costs.primary_costs,
        );
        info!(
            "Pruned {pruned} dominated primary anchors, {} remain",
            self.primary.len()
        );
        pruned
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
        self.enumerate_primaries_from(0)
    }

    /// Enumerates the primary anchors starting from the given index.
    fn enumerate_primaries_from(
        &self,
        first_index: usize,
    ) -> impl Iterator<Item = (AnchorIndex, PrimaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.primary[first_index..]
            .iter()
            .copied()
            .enumerate()
            .map(move |(index, anchor)| ((first_index + index).into(), anchor))
    }

    /// Returns the index of the first primary anchor whose start is at least the given coordinates.
    ///
    /// The primary anchors are sorted by their start.
    fn first_primary_index_starting_at_or_after(
        &self,
        coordinates: PrimaryAlignmentCoordinates,
    ) -> usize {
        self.primary
            .partition_point(|anchor| anchor.start() < coordinates)
    }

    /// Enumerates the primary anchors whose start is at least the given coordinates
    /// in the order of [`PrimaryAlignmentCoordinates`], which compares the first ordinate first.
    ///
    /// For coordinates `(a, 0)`, these are a superset of the anchors that start at or after `a` in the first sequence.
    pub fn enumerate_primaries_starting_at_or_after(
        &self,
        coordinates: PrimaryAlignmentCoordinates,
    ) -> impl Iterator<Item = (AnchorIndex, PrimaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.enumerate_primaries_from(self.first_primary_index_starting_at_or_after(coordinates))
    }

    /// Enumerates the primary anchors that start exactly at the given coordinates.
    pub fn enumerate_primaries_starting_at(
        &self,
        coordinates: PrimaryAlignmentCoordinates,
    ) -> impl Iterator<Item = (AnchorIndex, PrimaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.enumerate_primaries_starting_at_or_after(coordinates)
            .take_while(move |(_, anchor)| anchor.start() == coordinates)
    }

    #[cfg(test)]
    pub fn primary_index_from_range(
        &self,
        range: crate::alignment::coordinates::range::PrimaryAlignmentRange,
    ) -> Option<AnchorIndex>
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
        self.enumerate_secondaries_from(0, ts_kind)
    }

    /// Enumerates the secondary anchors of the given kind starting from the given index.
    fn enumerate_secondaries_from(
        &self,
        first_index: usize,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (AnchorIndex, SecondaryAnchor<Cost>)>
    where
        Cost: Copy,
    {
        self.secondary_anchor_vec(ts_kind)[first_index..]
            .iter()
            .copied()
            .enumerate()
            .map(move |(index, anchor)| ((first_index + index).into(), anchor))
    }

    /// Returns the index of the first secondary anchor of the given kind whose start is at least the given coordinates.
    ///
    /// The secondary anchors are sorted by their start.
    fn first_secondary_index_starting_at_or_after(
        &self,
        coordinates: AnySecondaryAlignmentCoordinates,
        ts_kind: TsKind,
    ) -> usize {
        self.secondary_anchor_vec(ts_kind)
            .partition_point(|anchor| anchor.start() < coordinates)
    }

    #[cfg(test)]
    pub fn secondary_index_from_range(
        &self,
        start: crate::alignment::coordinates::range::SpecificSecondaryAlignmentRange,
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
