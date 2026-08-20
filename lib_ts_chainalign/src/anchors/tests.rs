use crate::{
    alignment::{sequences::AlignmentSequences, ts_kind::TsKind},
    anchors::{Anchors, PrimaryAnchor, SecondaryAnchor},
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
        [(0, 0, 0), (2, 0, 0)].map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS11).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[(2, 2, 0), (4, 2, 0)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[(4, 0, 0), (4, 2, 0)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS22),
        &[(4, 0, 0), (3, 1, 0), (2, 2, 0)].map(SecondaryAnchor::from)
    );
}

#[test]
fn test_coordinates_rev() {
    let sequences = AlignmentSequences::new_complete(b"ACGT".to_vec(), b"ACAC".to_vec());
    let k = 2;

    let anchors = Anchors::new_exact(&sequences, k, &rc_fn);
    assert_eq!(
        anchors.primary,
        [(0, 0, 0), (0, 2, 0)].map(PrimaryAnchor::from)
    );
    assert!(anchors.secondary_anchor_vec(TsKind::TS22).is_empty());
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS21),
        &[(2, 2, 0), (4, 2, 0)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS12),
        &[(4, 0, 0), (4, 2, 0)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_anchor_vec(TsKind::TS11),
        &[(4, 0, 0), (3, 1, 0), (2, 2, 0)].map(SecondaryAnchor::from)
    );
}
