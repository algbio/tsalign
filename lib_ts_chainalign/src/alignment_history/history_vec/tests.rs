use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation::{self, *},
    history_vec::{AlignmentHistory, UnsignedIntAlignmentHistoryVec},
    test_util::{OPERATIONS, contains_anchor, longest_suffix_newest_first, reachable_histories},
};

/// Extends the empty history with the given operations, panicking if an extension is not possible.
fn history_from<AlignmentHistoryVec: AlignmentHistory>(
    operations: &[AlignmentHistoryOperation],
    anchor_k: u8,
    max_anchor_mutations: u8,
) -> AlignmentHistoryVec {
    operations
        .iter()
        .fold(AlignmentHistoryVec::default(), |history, operation| {
            history
                .try_extend(*operation, anchor_k, max_anchor_mutations)
                .unwrap()
        })
}

#[test]
fn required_capacity() {
    fn required_capacity(anchor_k: u8, max_anchor_mutations: u8) -> usize {
        UnsignedIntAlignmentHistoryVec::<u8>::required_capacity(anchor_k, max_anchor_mutations)
    }

    assert_eq!(required_capacity(8, 0), 7);
    assert_eq!(required_capacity(8, 1), 7);
    assert_eq!(required_capacity(8, 2), 8);
    assert_eq!(required_capacity(8, 3), 8);
    assert_eq!(required_capacity(8, 4), 9);

    assert_eq!(required_capacity(6, 0), 5);
    assert_eq!(required_capacity(6, 1), 5);
    assert_eq!(required_capacity(6, 2), 6);
    assert_eq!(required_capacity(6, 3), 6);
    assert_eq!(required_capacity(6, 4), 7);

    // The number of gap pairs is limited by the number of characters.
    assert_eq!(required_capacity(1, 0), 0);
    assert_eq!(required_capacity(1, 4), 0);
    assert_eq!(required_capacity(2, 4), 2);
    assert_eq!(required_capacity(3, 6), 4);
    assert_eq!(required_capacity(3, 10), 4);
}

#[test]
fn has_enough_capacity() {
    // `u8` can store three operations.
    assert!(UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
        4, 1
    ));
    assert!(UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
        3, 2
    ));
    assert!(!UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
        4, 2
    ));
    assert!(!UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
        5, 0
    ));
}

/// The maximum length of a reachable history is exactly the required capacity.
#[test]
fn required_capacity_is_tight() {
    for anchor_k in 1..=7 {
        for max_anchor_mutations in 0..=5 {
            let max_len = reachable_histories::<UnsignedIntAlignmentHistoryVec<u64>>(
                anchor_k,
                max_anchor_mutations,
            )
            .iter()
            .map(AlignmentHistory::len)
            .max()
            .unwrap();
            assert_eq!(
                max_len,
                UnsignedIntAlignmentHistoryVec::<u64>::required_capacity(
                    anchor_k,
                    max_anchor_mutations
                ),
                "anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
            );
        }
    }
}

#[test]
fn operations_are_newest_first() {
    let history: UnsignedIntAlignmentHistoryVec<u64> =
        history_from(&[Match, Substitution, GapInA], 10, 3);
    assert_eq!(
        history.operations().collect::<Vec<_>>(),
        vec![GapInA, Substitution, Match]
    );
    assert_eq!(history.len(), 3);
}

#[test]
fn drops_everything_up_to_the_oldest_mutation() {
    let operations = [Substitution, Match, Match, GapInA];

    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&operations, 10, 1);
    assert_eq!(
        history.operations().collect::<Vec<_>>(),
        vec![GapInA, Match, Match]
    );

    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&operations, 10, 2);
    assert_eq!(
        history.operations().collect::<Vec<_>>(),
        vec![GapInA, Match, Match, Substitution]
    );

    // With zero allowed mutations, only the trailing matches are kept.
    let history: UnsignedIntAlignmentHistoryVec<u64> =
        history_from(&[Match, Substitution, Match, Match], 10, 0);
    assert_eq!(history.operations().collect::<Vec<_>>(), vec![Match, Match]);
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[Match, Match, GapInB], 10, 0);
    assert!(history.is_empty());
}

#[test]
fn rejects_anchors() {
    // Three matches form an exact anchor with `anchor_k = 3`.
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[Match, Match], 3, 0);
    assert!(history.try_extend(Match, 3, 0).is_none());
    assert!(history.try_extend(Substitution, 3, 0).is_some());

    // With one allowed mutation, a substitution does not help.
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[Match, Match], 3, 1);
    assert!(history.try_extend(Substitution, 3, 1).is_none());
    // A gap in a covers a character of sequence b, so the subalignment covers three characters of sequence b.
    assert!(history.try_extend(GapInA, 3, 1).is_none());
    // Two gaps in a and b are enough.
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[GapInA, GapInB], 3, 1);
    assert!(history.try_extend(Match, 3, 1).is_some());

    // An anchor may consist only of gaps.
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[GapInB, GapInB], 3, 3);
    assert!(history.try_extend(GapInB, 3, 3).is_none());
    let history: UnsignedIntAlignmentHistoryVec<u64> = history_from(&[GapInB, GapInB], 3, 2);
    assert!(history.try_extend(GapInB, 3, 2).is_some());
}

/// Exhaustively compares [`AlignmentHistory::try_extend`] with a brute force check for anchors in the whole alignment.
fn check_try_extend_against_brute_force<AlignmentHistoryVec: AlignmentHistory>(
    anchor_k: u8,
    max_anchor_mutations: u8,
    max_alignment_len: usize,
) {
    fn recurse<AlignmentHistoryVec: AlignmentHistory>(
        alignment: &mut Vec<AlignmentHistoryOperation>,
        history: AlignmentHistoryVec,
        anchor_k: u8,
        max_anchor_mutations: u8,
        max_alignment_len: usize,
    ) {
        if alignment.len() == max_alignment_len {
            return;
        }

        for operation in OPERATIONS {
            alignment.push(operation);
            let extension = history.try_extend(operation, anchor_k, max_anchor_mutations);
            let has_anchor =
                contains_anchor(alignment, anchor_k.into(), max_anchor_mutations.into());
            assert_eq!(
                extension.is_none(),
                has_anchor,
                "alignment = {alignment:?}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
            );

            if let Some(extension) = extension {
                assert_eq!(
                    extension.operations().collect::<Vec<_>>(),
                    longest_suffix_newest_first(alignment, max_anchor_mutations.into()),
                    "alignment = {alignment:?}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
                );
                assert!(
                    extension.len()
                        <= AlignmentHistoryVec::required_capacity(anchor_k, max_anchor_mutations)
                );

                recurse(
                    alignment,
                    extension,
                    anchor_k,
                    max_anchor_mutations,
                    max_alignment_len,
                );
            }
            alignment.pop();
        }
    }

    recurse(
        &mut Vec::new(),
        AlignmentHistoryVec::default(),
        anchor_k,
        max_anchor_mutations,
        max_alignment_len,
    );
}

#[test]
fn try_extend_agrees_with_brute_force() {
    for anchor_k in 1..=5 {
        for max_anchor_mutations in 0..anchor_k {
            check_try_extend_against_brute_force::<UnsignedIntAlignmentHistoryVec<u64>>(
                anchor_k,
                max_anchor_mutations,
                7,
            );
        }
    }

    check_try_extend_against_brute_force::<UnsignedIntAlignmentHistoryVec<u64>>(7, 2, 9);
    check_try_extend_against_brute_force::<UnsignedIntAlignmentHistoryVec<u64>>(3, 3, 7);
}

#[test]
fn try_extend_agrees_with_brute_force_small_int() {
    for (anchor_k, max_anchor_mutations) in [(2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1)] {
        assert!(UnsignedIntAlignmentHistoryVec::<u8>::has_enough_capacity(
            anchor_k,
            max_anchor_mutations
        ));
        check_try_extend_against_brute_force::<UnsignedIntAlignmentHistoryVec<u8>>(
            anchor_k,
            max_anchor_mutations,
            7,
        );
    }
}
