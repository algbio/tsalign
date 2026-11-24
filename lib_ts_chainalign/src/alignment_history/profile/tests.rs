use std::collections::HashMap;

use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation::{self, *},
    profile::HistoryProfile,
};

const OPERATIONS: [AlignmentHistoryOperation; 4] = [Match, Substitution, GapInA, GapInB];

/// Returns true if the given alignment contains a subalignment with at most `max_anchor_mutations` mutations
/// that covers at least `anchor_k` characters of one of the sequences.
fn brute_force_has_anchor(
    operations: &[AlignmentHistoryOperation],
    anchor_k: usize,
    max_anchor_mutations: usize,
) -> bool {
    (0..operations.len()).any(|start| {
        let mut len_a = 0;
        let mut len_b = 0;
        let mut mutations = 0;
        for operation in &operations[start..] {
            len_a += operation.len_a();
            len_b += operation.len_b();
            mutations += operation.mutations();
            if mutations > max_anchor_mutations {
                return false;
            }
            if len_a >= anchor_k || len_b >= anchor_k {
                return true;
            }
        }
        false
    })
}

/// Returns the longest suffix of the given alignment with at most `max_mutations` mutations.
fn longest_suffix(
    operations: &[AlignmentHistoryOperation],
    max_mutations: usize,
) -> &[AlignmentHistoryOperation] {
    let mut mutations = 0;
    let mut start = operations.len();
    while start > 0 {
        mutations += operations[start - 1].mutations();
        if mutations > max_mutations {
            break;
        }
        start -= 1;
    }
    &operations[start..]
}

/// Returns `(len_a(0..=M), len_b(0..=M))` computed from the definition.
fn brute_force_profile(
    operations: &[AlignmentHistoryOperation],
    max_index: usize,
) -> (Vec<usize>, Vec<usize>) {
    (0..=max_index)
        .map(|j| {
            let suffix = longest_suffix(operations, j);
            (
                suffix
                    .iter()
                    .map(|operation| operation.len_a())
                    .sum::<usize>(),
                suffix
                    .iter()
                    .map(|operation| operation.len_b())
                    .sum::<usize>(),
            )
        })
        .unzip()
}

struct Checker<'a> {
    anchor_k: u8,
    max_anchor_mutations: u8,
    max_len: usize,
    /// Maps each reached profile to the longest suffix with at most `max_anchor_mutations` mutations.
    profile_to_suffix: &'a mut HashMap<HistoryProfile, Vec<AlignmentHistoryOperation>>,
    checked_extensions: usize,
}

impl Checker<'_> {
    fn check(&mut self, operations: &mut Vec<AlignmentHistoryOperation>, profile: HistoryProfile) {
        let anchor_k = usize::from(self.anchor_k);
        let max_anchor_mutations = usize::from(self.max_anchor_mutations);
        let max_index = max_anchor_mutations.min(2 * anchor_k - 1);

        // The profile matches its definition.
        let (len_a, len_b) = brute_force_profile(operations, max_index);
        assert_eq!(
            profile.len_a().collect::<Vec<_>>(),
            len_a,
            "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
        );
        assert_eq!(
            profile.len_b().collect::<Vec<_>>(),
            len_b,
            "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
        );

        // The profile is in bijection with the longest suffix with at most `max_anchor_mutations` mutations.
        let suffix = longest_suffix(operations, max_anchor_mutations).to_vec();
        if let Some(previous_suffix) = self.profile_to_suffix.get(&profile) {
            assert_eq!(
                previous_suffix, &suffix,
                "profile {profile} k = {anchor_k} m = {max_anchor_mutations}"
            );
        } else {
            self.profile_to_suffix.insert(profile, suffix);
        }

        if operations.len() == self.max_len {
            return;
        }

        for operation in OPERATIONS {
            operations.push(operation);
            let expected_valid =
                !brute_force_has_anchor(operations, anchor_k, max_anchor_mutations);
            let extension = profile.try_extend(operation, self.anchor_k, self.max_anchor_mutations);
            assert_eq!(
                extension.is_some(),
                expected_valid,
                "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
            );
            self.checked_extensions += 1;

            if let Some(extension) = extension {
                self.check(operations, extension);
            }
            operations.pop();
        }
    }
}

/// Check [`HistoryProfile::try_extend`] against a brute-force anchor check on all short alignments,
/// and check that the profile is computed according to its definition.
#[test]
fn try_extend_matches_brute_force() {
    for anchor_k in 1..=4u8 {
        for max_anchor_mutations in 0..=8u8 {
            let mut profile_to_suffix = HashMap::new();
            let max_len = usize::from(2 * anchor_k + 3).min(10);
            let mut checker = Checker {
                anchor_k,
                max_anchor_mutations,
                max_len,
                profile_to_suffix: &mut profile_to_suffix,
                checked_extensions: 0,
            };
            checker.check(
                &mut Vec::new(),
                HistoryProfile::new_empty(anchor_k, max_anchor_mutations),
            );
            assert!(checker.checked_extensions > 0);
        }
    }
}

/// Same as [`try_extend_matches_brute_force`], but with larger `anchor_k` and fewer, longer alignments.
#[test]
fn try_extend_matches_brute_force_larger_k() {
    for (anchor_k, max_anchor_mutations) in [(6u8, 1u8), (7, 2), (8, 3)] {
        let mut profile_to_suffix = HashMap::new();
        let mut checker = Checker {
            anchor_k,
            max_anchor_mutations,
            max_len: 10,
            profile_to_suffix: &mut profile_to_suffix,
            checked_extensions: 0,
        };
        checker.check(
            &mut Vec::new(),
            HistoryProfile::new_empty(anchor_k, max_anchor_mutations),
        );
    }
}

#[test]
fn extension_rules() {
    let profile = HistoryProfile::new_empty(10, 2);
    assert_eq!(profile.to_string(), "a[0,0,0]b[0,0,0]");

    let profile = profile.try_extend(Match, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[1,1,1]b[1,1,1]");
    let profile = profile.try_extend(GapInA, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[0,1,1]b[0,2,2]");
    let profile = profile.try_extend(Match, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[1,2,2]b[1,3,3]");
    let profile = profile.try_extend(GapInB, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[0,2,3]b[0,1,3]");
    let profile = profile.try_extend(Substitution, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[0,1,3]b[0,1,2]");
    let profile = profile.try_extend(Substitution, 10, 2).unwrap();
    assert_eq!(profile.to_string(), "a[0,1,2]b[0,1,2]");
}

#[test]
fn match_run_limit() {
    let mut profile = HistoryProfile::new_empty(5, 0);
    for _ in 0..4 {
        profile = profile.try_extend(Match, 5, 0).unwrap();
    }
    assert!(profile.try_extend(Match, 5, 0).is_none());
    assert!(profile.try_extend(Substitution, 5, 0).is_some());
}

/// All parameters that were supported by the previous alignment history window (`anchor_k + ceil(max_anchor_mutations / 2) <= 61`)
/// are supported, and the documented range is exact.
#[test]
fn supported_range() {
    for anchor_k in 1..=u8::MAX {
        for max_anchor_mutations in 0..=u8::MAX {
            let k = u32::from(anchor_k);
            let m = u32::from(max_anchor_mutations);
            let supported = HistoryProfile::is_supported(anchor_k, max_anchor_mutations);
            assert_eq!(supported, k + m.min(2 * k - 1) <= 127, "k = {k}, m = {m}");
            if k + m.div_ceil(2) <= 61 {
                assert!(supported, "k = {k}, m = {m}");
            }
            if supported {
                // The largest encodings do not overflow.
                let mut profile = HistoryProfile::new_empty(anchor_k, max_anchor_mutations);
                for _ in 1..k {
                    profile = profile
                        .try_extend(Match, anchor_k, max_anchor_mutations)
                        .unwrap();
                }
                assert!(
                    profile
                        .try_extend(Match, anchor_k, max_anchor_mutations)
                        .is_none()
                );
            }
        }
    }
    assert!(!HistoryProfile::is_supported(0, 0));
}
