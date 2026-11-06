use std::collections::HashSet;

use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation::{self, *},
    history_vec::{AlignmentHistory, UnsignedIntAlignmentHistoryVec},
    profile::{HistoryProfile, bytewise_at_most},
    test_util::{OPERATIONS, Rng, reachable_histories},
};

type History = UnsignedIntAlignmentHistoryVec<u64>;

fn profile_of(history: &History, max_anchor_mutations: u8) -> HistoryProfile {
    HistoryProfile::new(history.operations(), max_anchor_mutations).unwrap()
}

/// Returns the profile entries `(len_a(j), len_b(j))` for `j` in `0..=max_anchor_mutations`.
fn entries(profile: &HistoryProfile, max_anchor_mutations: u8) -> Vec<(u8, u8)> {
    (0..=max_anchor_mutations)
        .map(|j| (profile.len_a(j), profile.len_b(j)))
        .collect()
}

fn profile_from_entries(entries: &[(u8, u8)]) -> HistoryProfile {
    let mut len_a = [0; 8];
    let mut len_b = [0; 8];
    for (j, (a, b)) in entries.iter().enumerate() {
        len_a[j] = *a;
        len_b[j] = *b;
    }
    HistoryProfile {
        len_a: u64::from_le_bytes(len_a),
        len_b: u64::from_le_bytes(len_b),
    }
}

/// The extension rule from the module documentation, applied directly to the profile entries.
fn extend_entries(entries: &[(u8, u8)], operation: AlignmentHistoryOperation) -> Vec<(u8, u8)> {
    let len_a = u8::try_from(operation.len_a()).unwrap();
    let len_b = u8::try_from(operation.len_b()).unwrap();

    if operation.mutations() == 0 {
        entries
            .iter()
            .map(|(a, b)| (a + len_a, b + len_b))
            .collect()
    } else {
        std::iter::once((0, 0))
            .chain(
                entries[..entries.len() - 1]
                    .iter()
                    .map(|(a, b)| (a + len_a, b + len_b)),
            )
            .collect()
    }
}

#[test]
fn empty_history_has_zero_profile() {
    for max_anchor_mutations in 0..=HistoryProfile::MAX_SUPPORTED_ANCHOR_MUTATIONS {
        assert_eq!(
            HistoryProfile::new(std::iter::empty(), max_anchor_mutations),
            Some(HistoryProfile::default())
        );
    }
}

#[test]
fn example() {
    // From oldest to newest, the alignment is `M M S M I M M`.
    let history = [Match, Match, GapInA, Match, Substitution, Match, Match];

    let profile = HistoryProfile::new(history, 2).unwrap();
    assert_eq!(entries(&profile, 2), vec![(2, 2), (3, 4), (6, 7)]);

    let profile = HistoryProfile::new(history, 1).unwrap();
    assert_eq!(entries(&profile, 1), vec![(2, 2), (3, 4)]);

    let profile = HistoryProfile::new(history, 0).unwrap();
    assert_eq!(entries(&profile, 0), vec![(2, 2)]);

    // Starting with a mutation.
    let profile = HistoryProfile::new([GapInB, Match, GapInB], 2).unwrap();
    assert_eq!(entries(&profile, 2), vec![(0, 0), (2, 1), (3, 1)]);
}

#[test]
fn ignores_operations_beyond_max_anchor_mutations() {
    assert_eq!(
        HistoryProfile::new([Match, Substitution, Match, GapInA, Match, Match], 1),
        HistoryProfile::new([Match, Substitution, Match], 1),
    );
    assert_ne!(
        HistoryProfile::new([Match, Substitution, Match, GapInA, Match, Match], 2),
        HistoryProfile::new([Match, Substitution, Match], 2),
    );
}

#[test]
fn unsupported_max_anchor_mutations() {
    assert!(HistoryProfile::new([Match], HistoryProfile::MAX_SUPPORTED_ANCHOR_MUTATIONS).is_some());
    assert!(
        HistoryProfile::new([Match], HistoryProfile::MAX_SUPPORTED_ANCHOR_MUTATIONS + 1).is_none()
    );
}

#[test]
fn maximum_supported_values() {
    let history = vec![Match; HistoryProfile::MAX_SUPPORTED_LEN];
    let profile = HistoryProfile::new(history, 7).unwrap();
    assert_eq!(
        entries(&profile, 7),
        vec![(127, 127); 8],
        "all entries cover the whole history"
    );
    assert!(HistoryProfile::default().is_at_most(&profile));
    assert!(!profile.is_at_most(&HistoryProfile::default()));
    assert!(profile.is_at_most(&profile));
}

#[test]
#[should_panic]
fn too_long_history() {
    HistoryProfile::new(vec![Match; HistoryProfile::MAX_SUPPORTED_LEN + 1], 0);
}

#[test]
fn bytewise_at_most_agrees_with_naive() {
    fn naive(x: u64, y: u64) -> bool {
        x.to_le_bytes()
            .into_iter()
            .zip(y.to_le_bytes())
            .all(|(x, y)| x <= y)
    }

    let mut rng = Rng::new(0);
    for _ in 0..100_000 {
        // Draw bytes from a small range to get many equal and comparable bytes.
        let max = [2, 4, 128][rng.below(3)];
        let x = u64::from_le_bytes(std::array::from_fn(|_| rng.below(max) as u8));
        let y = u64::from_le_bytes(std::array::from_fn(|_| rng.below(max) as u8));
        assert_eq!(
            bytewise_at_most(x, y),
            naive(x, y),
            "x = {x:#x}, y = {y:#x}"
        );
        assert_eq!(
            bytewise_at_most(y, x),
            naive(y, x),
            "x = {x:#x}, y = {y:#x}"
        );
        assert!(bytewise_at_most(x, x));
    }

    for (x, y) in [
        (0, 0),
        (0, 0x7f7f_7f7f_7f7f_7f7f),
        (0x7f, 0x7e),
        (0x7e, 0x7f),
    ] {
        assert_eq!(
            bytewise_at_most(x, y),
            naive(x, y),
            "x = {x:#x}, y = {y:#x}"
        );
    }
}

#[test]
fn is_at_most_is_a_partial_order() {
    let (anchor_k, max_anchor_mutations) = (5, 2);
    let profiles: Vec<_> = reachable_histories::<History>(anchor_k, max_anchor_mutations)
        .iter()
        .map(|history| profile_of(history, max_anchor_mutations))
        .collect();

    for p in &profiles {
        assert!(p.is_at_most(p));
        for q in &profiles {
            if p.is_at_most(q) && q.is_at_most(p) {
                assert_eq!(p, q);
            }
            // Check against the entries.
            let naive = entries(p, max_anchor_mutations)
                .into_iter()
                .zip(entries(q, max_anchor_mutations))
                .all(|((pa, pb), (qa, qb))| pa <= qa && pb <= qb);
            assert_eq!(p.is_at_most(q), naive);
        }
    }

    // Transitivity on a sample.
    let mut rng = Rng::new(1);
    for _ in 0..100_000 {
        let p = &profiles[rng.below(profiles.len())];
        let q = &profiles[rng.below(profiles.len())];
        let r = &profiles[rng.below(profiles.len())];
        if p.is_at_most(q) && q.is_at_most(r) {
            assert!(p.is_at_most(r));
        }
    }
}

/// The profile determines the history, so no information is lost by using profiles instead of histories.
#[test]
fn profile_is_injective() {
    for anchor_k in 1..=7 {
        for max_anchor_mutations in 0..anchor_k.min(5) {
            let histories = reachable_histories::<History>(anchor_k, max_anchor_mutations);
            let profiles: HashSet<_> = histories
                .iter()
                .map(|history| profile_of(history, max_anchor_mutations))
                .collect();
            assert_eq!(
                histories.len(),
                profiles.len(),
                "anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
            );
        }
    }
}

/// Checks the validity and extension properties of the module documentation against [`AlignmentHistory::try_extend`].
#[test]
fn extension_rule_and_validity() {
    for anchor_k in 1..=7 {
        for max_anchor_mutations in 0..anchor_k.min(6) {
            for history in reachable_histories::<History>(anchor_k, max_anchor_mutations) {
                let profile_entries = entries(
                    &profile_of(&history, max_anchor_mutations),
                    max_anchor_mutations,
                );
                // The history itself is valid.
                let (last_a, last_b) = *profile_entries.last().unwrap();
                assert!(last_a < anchor_k && last_b < anchor_k);

                for operation in OPERATIONS {
                    let extended_entries = extend_entries(&profile_entries, operation);
                    let (last_a, last_b) = *extended_entries.last().unwrap();
                    let valid_by_profile = last_a < anchor_k && last_b < anchor_k;

                    match history.try_extend(operation, anchor_k, max_anchor_mutations) {
                        Some(extension) => {
                            assert!(valid_by_profile, "history = {history}, {operation}");
                            assert_eq!(
                                profile_of(&extension, max_anchor_mutations),
                                profile_from_entries(&extended_entries),
                                "history = {history}, {operation}"
                            );
                        }
                        None => assert!(!valid_by_profile, "history = {history}, {operation}"),
                    }
                }
            }
        }
    }
}

/// The central property: if the profile of `A` is at most the profile of `B`, then every valid continuation of `B` is a valid continuation of `A`.
#[test]
fn domination_preserves_valid_continuations() {
    fn check(a: History, b: History, anchor_k: u8, max_anchor_mutations: u8, remaining_len: usize) {
        assert!(
            profile_of(&a, max_anchor_mutations).is_at_most(&profile_of(&b, max_anchor_mutations))
        );
        if remaining_len == 0 {
            return;
        }

        for operation in OPERATIONS {
            if let Some(b_extension) = b.try_extend(operation, anchor_k, max_anchor_mutations) {
                let a_extension = a
                    .try_extend(operation, anchor_k, max_anchor_mutations)
                    .unwrap_or_else(|| {
                        panic!("{b} extends with {operation}, but the dominating {a} does not")
                    });
                check(
                    a_extension,
                    b_extension,
                    anchor_k,
                    max_anchor_mutations,
                    remaining_len - 1,
                );
            }
        }
    }

    for (anchor_k, max_anchor_mutations, continuation_len) in [
        (1, 0, 4),
        (2, 0, 5),
        (2, 1, 5),
        (3, 1, 4),
        (4, 1, 4),
        (4, 2, 3),
        (5, 2, 3),
        (4, 3, 3),
    ] {
        let histories = reachable_histories::<History>(anchor_k, max_anchor_mutations);
        let profiles: Vec<_> = histories
            .iter()
            .map(|history| profile_of(history, max_anchor_mutations))
            .collect();
        let mut strictly_dominating_pairs = 0;

        for (a, a_profile) in histories.iter().zip(&profiles) {
            for (b, b_profile) in histories.iter().zip(&profiles) {
                if a_profile.is_at_most(b_profile) {
                    if a != b {
                        strictly_dominating_pairs += 1;
                    }
                    check(*a, *b, anchor_k, max_anchor_mutations, continuation_len);
                }
            }
        }

        if anchor_k > 1 {
            assert!(
                strictly_dominating_pairs > 0,
                "anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}"
            );
        }
    }
}
