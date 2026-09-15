use generic_a_star::cost::U16Cost;

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
fn test_coordinates_inexact_1() {
    // TODO: the expected values in this test are not complete.
    let sequences = AlignmentSequences::new_complete(b"ACAC".to_vec(), b"ACGT".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
    assert_eq!(
        anchors.primary,
        [
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 0..2),
                U16Cost::from(0u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 0..1),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(0..2, 1..2),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(1..3, 0..1),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(2..4, 0..2),
                U16Cost::from(0u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(2..4, 0..1),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(2..4, 1..2),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(1..3, 1..2),
                U16Cost::from(3u16)
            ),
            (
                PrimaryAlignmentRange::new_from_ranges(1..3, 1..3),
                U16Cost::from(2u16)
            ),
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
fn test_coordinates_rev_inexact_1() {
    // TODO: the expected values in this test are not complete.
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;
    let costs = GapAffineCosts::new(2u16, 3, 1).into_costs::<U16Cost>();

    let anchors = Anchors::new_inexact(&sequences, k, 1, &costs, &rc_fn);
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
