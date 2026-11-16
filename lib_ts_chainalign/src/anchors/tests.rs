use std::{
    collections::HashSet,
    fmt::{Debug, Display},
};

use generic_a_star::cost::U16Cost;
use itertools::Itertools;
use num_traits::Zero;

use crate::{
    alignment::{
        coordinates::{
            AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates,
            range::{AnySecondaryAlignmentRange, PrimaryAlignmentRange},
        },
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    anchors::{Anchors, PrimaryAnchor, SecondaryAnchor, index::AnchorIndex},
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
};

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => unimplemented!("Unsupported character {c}"),
    }
}

#[test]
fn coordinates() {
    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;

    let anchors = Anchors::new_exact(&sequences, k, &rc_fn);
    assert_eq!(
        anchors.primary,
        [
            (PrimaryAlignmentRange::new_from_ranges(0..2, 0..2), 0),
            (PrimaryAlignmentRange::new_from_ranges(2..4, 0..2), 0)
        ]
        .map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS11).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS22),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
}

#[test]
fn coordinates_rev() {
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;

    let anchors = Anchors::new_exact(&sequences, k, &rc_fn);
    assert_eq!(
        anchors.primary,
        [
            (PrimaryAlignmentRange::new_from_ranges(0..2, 0..2), 0),
            (PrimaryAlignmentRange::new_from_ranges(0..2, 2..4), 0)
        ]
        .map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS22).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS11),
        &[
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0)
        ]
        .map(SecondaryAnchor::from)
    );
}

#[test]
fn coordinates_inexact_0() {
    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors = Anchors::new_inexact(&sequences, k, 0, &costs, &rc_fn);
    assert_eq!(
        anchors.primary,
        [
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(2..4, 0..2),
                U16Cost::from(0u16)
            )
        ]
        .map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS11).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS22),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
}

#[test]
fn coordinates_rev_inexact_0() {
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors = Anchors::new_inexact(&sequences, k, 0, &costs, &rc_fn);
    assert_eq!(
        anchors.primary,
        [
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS22).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS11),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
                U16Cost::from(0u16)
            ),
            (
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
                U16Cost::from(0u16)
            )
        ]
        .map(SecondaryAnchor::from)
    );
}

#[test]
fn coordinates_inexact_1_symmetry_small() {
    let sequences = AlignmentSequences::new_complete(b"AC".to_vec(), b"AC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    let primary: HashSet<_> = anchors.primary.iter().cloned().collect();
    let secondaries = anchors
        .secondaries
        .map(|secondary| secondary.iter().cloned().collect::<HashSet<_>>());

    let primary_r = primary
        .iter()
        .cloned()
        .map(|a| PrimaryAnchor::new(a.range().swap(), a.cost()))
        .sorted_unstable()
        .collect::<HashSet<_>>();
    let secondaries_r = [
        &secondaries[TsKind::swap_index(0)],
        &secondaries[TsKind::swap_index(1)],
        &secondaries[TsKind::swap_index(2)],
        &secondaries[TsKind::swap_index(3)],
    ]
    .map(|secondary| {
        secondary
            .iter()
            .cloned()
            .map(|a| SecondaryAnchor::new(a.range().swap(), a.cost()))
            .sorted_unstable()
            .collect::<HashSet<_>>()
    });

    assert_eq!(
        primary,
        primary_r,
        "\nP forward:\n{}\nP reverse:\n{}",
        primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        primary_r
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
    );
    for (ts_kind, (sf, sr)) in TsKind::iter().zip(secondaries.iter().zip(secondaries_r.iter())) {
        assert_eq!(
            sf,
            sr,
            "\nS{} forward:\n{}\nS{} reverse:\n{}",
            ts_kind.digits(),
            sf.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
            ts_kind.digits(),
            sr.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        );
    }
}

#[test]
fn coordinates_inexact_1_symmetry() {
    let sequences_f = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let sequences_r = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors_f = Anchors::new_inexact(&sequences_f, k, 1, &costs, &rc_fn);
    let anchors_r = Anchors::new_inexact(&sequences_r, k, 1, &costs, &rc_fn);
    let primary: HashSet<_> = anchors_f.primary.iter().cloned().collect();
    let secondaries = anchors_f
        .secondaries
        .map(|secondary| secondary.iter().cloned().collect::<HashSet<_>>());

    let primary_r = anchors_r
        .primary
        .iter()
        .cloned()
        .map(|a| PrimaryAnchor::new(a.range().swap(), a.cost()))
        .sorted_unstable()
        .collect::<HashSet<_>>();
    let secondaries_r = [
        &anchors_r.secondaries[TsKind::swap_index(0)],
        &anchors_r.secondaries[TsKind::swap_index(1)],
        &anchors_r.secondaries[TsKind::swap_index(2)],
        &anchors_r.secondaries[TsKind::swap_index(3)],
    ]
    .map(|secondary| secondary.iter().cloned().collect::<HashSet<_>>());

    assert_eq!(
        primary,
        primary_r,
        "\nP forward:\n{}\nP reverse:\n{}",
        primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        primary_r
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
    );
    for (ts_kind, (sf, sr)) in TsKind::iter().zip(secondaries.iter().zip(secondaries_r.iter())) {
        assert_eq!(
            sf,
            sr,
            "\nS{} forward:\n{}\nS{} reverse:\n{}",
            ts_kind.digits(),
            sf.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
            ts_kind.digits(),
            sr.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        );
    }
}

#[test]
fn coordinates_inexact_1() {
    fn assert_slice_eq<Item: Eq + Debug + ToString>(
        context: impl Display,
        actual: &[Item],
        expected: &[Item],
    ) {
        assert_eq!(
            actual,
            expected,
            "\n{context}: actual --- expected:\n{}",
            actual
                .iter()
                .zip_longest(expected.iter())
                .fold(String::new(), |acc, a| {
                    let (l, r) = a.left_and_right();
                    let l = l
                        .map(ToString::to_string)
                        .unwrap_or_else(|| "None".to_string());
                    let r = r
                        .map(ToString::to_string)
                        .unwrap_or_else(|| "None".to_string());
                    acc + &l + if l == r { " == " } else { " != " } + &r + "\n"
                }),
        );
    }

    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    let mut expected_primary = [
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 1..3),
            U16Cost::from(2u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(3..4, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..4, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(3..4, 1..3),
            U16Cost::from(3u16),
        ),
    ]
    .map(PrimaryAnchor::from);
    expected_primary.sort_unstable();
    assert_slice_eq("P", &anchors.primary, &expected_primary);

    assert_slice_eq("S11", anchors.secondary_anchor_vec(TsKind::TS11), &[]);

    // TGTG ACGT
    let mut expected_secondary_12 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(2u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(1, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_12.sort_unstable();
    assert_slice_eq(
        "S12",
        anchors.secondary_anchor_vec(TsKind::TS12),
        &expected_secondary_12,
    );

    // TGCA ACAC
    let mut expected_secondary_21 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(2u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_21.sort_unstable();
    assert_slice_eq(
        "S21",
        anchors.secondary_anchor_vec(TsKind::TS21),
        &expected_secondary_21,
    );

    // TGCA ACGT
    let mut expected_secondary_22 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(1, 0, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_22.sort_unstable();
    assert_slice_eq(
        "S22",
        anchors.secondary_anchor_vec(TsKind::TS22),
        &expected_secondary_22,
    );
}

#[test]
fn coordinates_inexact_1_rev() {
    fn assert_slice_eq<Item: Eq + Debug + ToString>(
        context: impl Display,
        actual: &[Item],
        expected: &[Item],
    ) {
        assert_eq!(
            actual,
            expected,
            "\n{context}: actual --- expected:\n{}",
            actual
                .iter()
                .zip_longest(expected.iter())
                .fold(String::new(), |acc, a| {
                    let (l, r) = a.left_and_right();
                    let l = l
                        .map(ToString::to_string)
                        .unwrap_or_else(|| "None".to_string());
                    let r = r
                        .map(ToString::to_string)
                        .unwrap_or_else(|| "None".to_string());
                    acc + &l + if l == r { " == " } else { " != " } + &r + "\n"
                }),
        );
    }

    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();
    let costs = AlignmentCosts::new(
        costs.clone(),
        costs,
        U16Cost::zero().into(),
        TsLimits::new_unlimited(),
    );

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    let mut expected_primary = [
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..3, 1..3),
            U16Cost::from(2u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(3..4, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..4, 0..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(3..4, 1..3),
            U16Cost::from(3u16),
        ),
    ]
    .map(PrimaryAnchor::from)
    .map(PrimaryAnchor::swap);
    expected_primary.sort_unstable();
    assert_slice_eq("P", &anchors.primary, &expected_primary);

    assert_slice_eq("S22", anchors.secondary_anchor_vec(TsKind::TS22), &[]);

    // TGTG ACGT
    let mut expected_secondary_12 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(2u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(1, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_12.sort_unstable();
    assert_slice_eq(
        "S21",
        anchors.secondary_anchor_vec(TsKind::TS21),
        &expected_secondary_12,
    );

    // TGCA ACAC
    let mut expected_secondary_21 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(2u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_21.sort_unstable();
    assert_slice_eq(
        "S12",
        anchors.secondary_anchor_vec(TsKind::TS12),
        &expected_secondary_21,
    );

    // TGCA ACGT
    let mut expected_secondary_22 = [
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 0, 3..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..3),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 1, 2..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 0..1),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 2, 1..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 3, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 0..2),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(4, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 2, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 1..3),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(3, 0, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(2, 1, 2..4),
            U16Cost::from(3u16),
        ),
        (
            AnySecondaryAlignmentRange::new_from_ranges(1, 0, 2..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(SecondaryAnchor::from);
    expected_secondary_22.sort_unstable();
    assert_slice_eq(
        "S11",
        anchors.secondary_anchor_vec(TsKind::TS11),
        &expected_secondary_22,
    );
}

/// Computes the expected result of the reverse lookup:
/// every anchor whose start is among the coordinates gets the minimum cost of these coordinates.
fn expected_reverse_lookup<Coordinates: Copy + Eq>(
    anchor_starts: impl IntoIterator<Item = (AnchorIndex, Coordinates)>,
    coordinates: &[(Coordinates, U16Cost)],
) -> Vec<(AnchorIndex, U16Cost)> {
    anchor_starts
        .into_iter()
        .filter_map(|(index, start)| {
            coordinates
                .iter()
                .filter(|(coordinates, _)| *coordinates == start)
                .map(|(_, cost)| *cost)
                .min()
                .map(|cost| (index, cost))
        })
        .sorted()
        .collect()
}

#[test]
fn primary_reverse_lookup_duplicate_coordinates_and_starts() {
    let mut primary = [
        // Anchor before the first coordinate.
        (PrimaryAlignmentRange::new_from_ranges(0..2, 0..2), 0),
        // Three anchors with the same start but different ends, and only two coordinates.
        (PrimaryAlignmentRange::new_from_ranges(1..3, 1..3), 1),
        (PrimaryAlignmentRange::new_from_ranges(1..2, 1..3), 1),
        (PrimaryAlignmentRange::new_from_ranges(1..3, 1..2), 2),
        // Anchor without matching coordinate.
        (PrimaryAlignmentRange::new_from_ranges(2..4, 3..5), 0),
        // Anchor with a single matching coordinate.
        (PrimaryAlignmentRange::new_from_ranges(3..5, 3..5), 0),
        // Two anchors with the same start and three coordinates.
        (PrimaryAlignmentRange::new_from_ranges(4..6, 4..6), 0),
        (PrimaryAlignmentRange::new_from_ranges(4..6, 4..7), 1),
        // Anchor after the last coordinate.
        (PrimaryAlignmentRange::new_from_ranges(7..9, 7..9), 0),
    ]
    .map(|(range, cost)| PrimaryAnchor::new(range, U16Cost::from(cost as u16)))
    .to_vec();
    primary.sort_unstable();
    let anchors = Anchors {
        primary,
        secondaries: Default::default(),
    };

    let mut coordinates = [
        ((0, 1), 3),
        ((1, 1), 2),
        ((1, 1), 4),
        ((1, 2), 9),
        ((3, 3), 5),
        ((4, 4), 6),
        ((4, 4), 1),
        ((4, 4), 8),
        ((5, 5), 0),
    ]
    .map(|((a, b), cost)| {
        (
            PrimaryAlignmentCoordinates::new(a, b),
            U16Cost::from(cost as u16),
        )
    })
    .to_vec();
    // The reverse lookup requires the coordinates to be sorted by coordinates and then by cost.
    coordinates.sort_unstable();

    let actual: Vec<_> = anchors
        .primary_anchor_to_index_iter(coordinates.iter().copied())
        .sorted()
        .collect();
    let expected = expected_reverse_lookup(
        anchors
            .enumerate_primaries()
            .map(|(index, anchor)| (index, anchor.start())),
        &coordinates,
    );
    assert_eq!(expected.len(), 6);
    assert_eq!(actual, expected);
}

#[test]
fn secondary_reverse_lookup_duplicate_coordinates_and_starts() {
    for ts_kind in TsKind::iter() {
        let mut secondary = [
            // Anchor before the first coordinate.
            ((10, 8, 0..2), 0),
            // Three anchors with the same start but different ends, and only two coordinates.
            ((8, 6, 8..10), 0),
            ((8, 7, 8..10), 1),
            ((8, 6, 8..9), 1),
            // Anchor without matching coordinate.
            ((8, 6, 9..11), 0),
            // Anchor with a single matching coordinate.
            ((7, 5, 3..5), 0),
            // Nine anchors with the same start and three coordinates (the reproduction from the e2e tests).
            ((5, 2, 11..14), 0),
            ((5, 2, 11..13), 1),
            ((5, 2, 11..12), 2),
            ((5, 3, 11..14), 1),
            ((5, 3, 11..13), 0),
            ((5, 3, 11..12), 1),
            ((5, 4, 11..14), 2),
            ((5, 4, 11..13), 1),
            ((5, 4, 11..12), 1),
            // Anchor after the last coordinate.
            ((2, 0, 0..2), 0),
        ]
        .map(|((ancestor_offset, ancestor_limit, descendant), cost)| {
            SecondaryAnchor::new_from_ranges(
                ancestor_offset,
                ancestor_limit,
                descendant,
                U16Cost::from(cost as u16),
            )
        })
        .to_vec();
        secondary.sort_unstable();
        let mut secondaries: [Vec<_>; 4] = Default::default();
        secondaries[ts_kind.index()] = secondary;
        let anchors = Anchors {
            primary: Vec::new(),
            secondaries,
        };

        let mut coordinates = [
            ((9, 0), 3),
            ((8, 8), 7),
            ((8, 8), 2),
            ((7, 3), 5),
            ((5, 11), 10),
            ((5, 11), 4),
            ((5, 11), 10),
            ((5, 12), 1),
            ((3, 0), 0),
        ]
        .map(|((ancestor, descendant), cost)| {
            (
                AnySecondaryAlignmentCoordinates::new(ancestor, descendant),
                U16Cost::from(cost as u16),
            )
        })
        .to_vec();
        // The reverse lookup requires the coordinates to be sorted by coordinates and then by cost.
        coordinates.sort_unstable();

        let actual: Vec<_> = anchors
            .secondary_anchor_to_index_iter(coordinates.iter().copied(), ts_kind)
            .sorted()
            .collect();
        let expected = expected_reverse_lookup(
            anchors
                .enumerate_secondaries(ts_kind)
                .map(|(index, anchor)| (index, anchor.start())),
            &coordinates,
        );
        assert_eq!(expected.len(), 13);
        assert_eq!(actual, expected, "S{}", ts_kind.digits());

        // Other TS kinds have no anchors.
        for other_ts_kind in TsKind::iter().filter(|other| *other != ts_kind) {
            assert_eq!(
                anchors
                    .secondary_anchor_to_index_iter(coordinates.iter().copied(), other_ts_kind)
                    .count(),
                0,
            );
        }
    }
}

#[test]
fn reverse_lookup_without_costs() {
    let mut primary = [
        PrimaryAlignmentRange::new_from_ranges(1..3, 1..3),
        PrimaryAlignmentRange::new_from_ranges(1..2, 1..3),
        PrimaryAlignmentRange::new_from_ranges(1..3, 1..2),
        PrimaryAlignmentRange::new_from_ranges(2..4, 2..4),
    ]
    .map(|range| PrimaryAnchor::new(range, U16Cost::zero()))
    .to_vec();
    primary.sort_unstable();
    let anchors = Anchors {
        primary,
        secondaries: Default::default(),
    };
    let coordinates = [(1, 1), (1, 1), (3, 3)].map(|(a, b)| PrimaryAlignmentCoordinates::new(a, b));

    let actual: Vec<_> = anchors
        .primary_anchor_to_index_iter(coordinates)
        .sorted()
        .collect();
    // Each anchor is reported exactly once, even though its start coordinates are duplicated.
    assert_eq!(actual, [0usize, 1, 2].map(AnchorIndex::from));
}
