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
    let primary_r = anchors
        .primary
        .iter()
        .cloned()
        .map(|a| PrimaryAnchor::new(a.range().swap(), a.cost()))
        .sorted_unstable()
        .collect::<Vec<_>>();
    assert_eq!(
        anchors.primary,
        primary_r,
        "\nforward:\n{}\nreverse:\n{}",
        anchors
            .primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        primary_r
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
    );
    todo!("secondaries");
    /*for (sf, sr) in anchors_f
        .secondaries
        .into_iter()
        .rev()
        .zip(anchors_r.secondaries)
    {
        assert_eq!(
            sf,
            sr,
            "\nforward:\n{}\nreverse:\n{}",
            sf.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
            sr.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        );
    }*/
}

#[test]
fn test_coordinates_inexact_1_symmetry() {
    let sequences_f = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let sequences_r = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

    let anchors_f = Anchors::new_inexact(&sequences_f, k, 1, &costs, &rc_fn);
    let anchors_r = Anchors::new_inexact(&sequences_r, k, 1, &costs, &rc_fn);
    let primary_r = anchors_r
        .primary
        .into_iter()
        .map(|a| PrimaryAnchor::new(a.range().swap(), a.cost()))
        .sorted_unstable()
        .collect::<Vec<_>>();
    assert_eq!(
        anchors_f.primary,
        primary_r,
        "\nforward:\n{}\nreverse:\n{}",
        anchors_f
            .primary
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        primary_r
            .iter()
            .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
    );
    for (sf, sr) in anchors_f
        .secondaries
        .into_iter()
        .rev()
        .zip(anchors_r.secondaries)
    {
        assert_eq!(
            sf,
            sr,
            "\nforward:\n{}\nreverse:\n{}",
            sf.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
            sr.iter()
                .fold(String::new(), |acc, a| acc + &a.to_string() + "\n"),
        );
    }
}

#[test]
fn test_coordinates_inexact_1() {
    // TODO: the expected values in this test are not complete.
    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    let expected_primary = [
        (
            PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(0..1, 1..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 0..0),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..4, 0..2),
            U16Cost::from(0u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(1..2, 1..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(2..3, 1..1),
            U16Cost::from(3u16),
        ),
        (
            PrimaryAlignmentRange::new_from_ranges(3..4, 1..1),
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
            PrimaryAlignmentRange::new_from_ranges(4..4, 2..3),
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
