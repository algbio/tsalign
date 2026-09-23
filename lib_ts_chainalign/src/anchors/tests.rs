use std::{
    collections::HashSet,
    fmt::{Debug, Display},
};

use generic_a_star::cost::U16Cost;
use itertools::Itertools;

use crate::{
    alignment::{
        coordinates::range::{AnySecondaryAlignmentRange, PrimaryAlignmentRange},
        sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    anchors::{Anchors, PrimaryAnchor, SecondaryAnchor},
    costs::GapAffineCosts,
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
fn test_coordinates() {
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
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0)
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
fn test_coordinates_rev() {
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
            (AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4), 0),
            (AnySecondaryAlignmentRange::new_from_ranges(4, 2, 2..4), 0)
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
fn test_coordinates_inexact_0() {
    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

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
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
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
fn test_coordinates_rev_inexact_0() {
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

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
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
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
fn test_coordinates_inexact_1_symmetry_small() {
    let sequences = AlignmentSequences::new_complete(b"AC".to_vec(), b"AC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

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
fn test_coordinates_inexact_1_symmetry() {
    let sequences_f = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let sequences_r = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

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
fn test_coordinates_inexact_1() {
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
fn test_coordinates_rev_inexact_1() {
    // TODO: the expected values in this test are not complete.
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    let expected_primary = [
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..0, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 2..4),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..1, 0..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..1, 1..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..1, 2..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..1, 3..4),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..2, 2..3),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 2..2),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 2..3),
            U16Cost::from(2u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 4..4),
            U16Cost::from(3u16),
        ),
    ]
    .map(PrimaryAnchor::from);
    assert_eq!(
        anchors.primary,
        expected_primary,
        "\nactual:\n{}\nexpected:\n{}",
        anchors
            .primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        expected_primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS22).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[
            (
                AnySecondaryAlignmentRange::new_from_ranges(2, 0, 2..4),
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
