use std::collections::HashMap;

use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation::{self, *},
    profile::{HistoryProfile, HistoryProfileParameters, MAX_PROFILE_ANCHOR_K, MAX_PROFILE_INDEX},
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
fn brute_force_profile(operations: &[AlignmentHistoryOperation], max_index: usize) -> Entries {
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

/// The number of stored profile entries per sequence minus one, i.e. `M` in the module documentation.
fn max_index(anchor_k: u8, max_anchor_mutations: u8) -> usize {
    usize::from(max_anchor_mutations).min(2 * usize::from(anchor_k) - 1)
}

/// The entries `(len_a(0..=M), len_b(0..=M))` of a profile.
type Entries = (Vec<usize>, Vec<usize>);

/// Returns the entries `(len_a(0..=M), len_b(0..=M))` of the profile.
fn entries(profile: HistoryProfile, parameters: HistoryProfileParameters) -> Entries {
    (
        profile.len_a(parameters).collect(),
        profile.len_b(parameters).collect(),
    )
}

struct Checker<'a> {
    anchor_k: u8,
    max_anchor_mutations: u8,
    parameters: HistoryProfileParameters,
    max_len: usize,
    /// Maps each reached profile to the longest suffix with at most `max_anchor_mutations` mutations.
    profile_to_suffix: &'a mut HashMap<HistoryProfile, Vec<AlignmentHistoryOperation>>,
    checked_extensions: usize,
}

impl Checker<'_> {
    fn check(&mut self, operations: &mut Vec<AlignmentHistoryOperation>, profile: HistoryProfile) {
        let anchor_k = usize::from(self.anchor_k);
        let max_anchor_mutations = usize::from(self.max_anchor_mutations);
        let max_index = max_index(self.anchor_k, self.max_anchor_mutations);

        // The profile matches its definition.
        assert_eq!(
            entries(profile, self.parameters),
            brute_force_profile(operations, max_index),
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
            let extension = profile.try_extend(operation, self.parameters);
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
            if !HistoryProfile::is_supported(anchor_k, max_anchor_mutations) {
                continue;
            }
            let mut profile_to_suffix = HashMap::new();
            let max_len = usize::from(2 * anchor_k + 3).min(10);
            let mut checker = Checker {
                anchor_k,
                max_anchor_mutations,
                parameters: HistoryProfileParameters::new(anchor_k, max_anchor_mutations),
                max_len,
                profile_to_suffix: &mut profile_to_suffix,
                checked_extensions: 0,
            };
            checker.check(&mut Vec::new(), HistoryProfile::empty());
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
            parameters: HistoryProfileParameters::new(anchor_k, max_anchor_mutations),
            max_len: 10,
            profile_to_suffix: &mut profile_to_suffix,
            checked_extensions: 0,
        };
        checker.check(&mut Vec::new(), HistoryProfile::empty());
    }
}

#[test]
fn extension_rules() {
    let parameters = HistoryProfileParameters::new(10, 2);
    let check = |profile: HistoryProfile, len_a: [usize; 3], len_b: [usize; 3]| {
        assert_eq!(
            entries(profile, parameters),
            (len_a.to_vec(), len_b.to_vec()),
            "{profile}"
        );
    };

    let profile = HistoryProfile::empty();
    check(profile, [0, 0, 0], [0, 0, 0]);
    assert_eq!(profile.to_string(), "a[0,0,0,0]b[0,0,0,0]");

    let profile = profile.try_extend(Match, parameters).unwrap();
    check(profile, [1, 1, 1], [1, 1, 1]);
    assert_eq!(profile.to_string(), "a[1,1,1,0]b[1,1,1,0]");
    let profile = profile.try_extend(GapInA, parameters).unwrap();
    check(profile, [0, 1, 1], [0, 2, 2]);
    let profile = profile.try_extend(Match, parameters).unwrap();
    check(profile, [1, 2, 2], [1, 3, 3]);
    let profile = profile.try_extend(GapInB, parameters).unwrap();
    check(profile, [0, 2, 3], [0, 1, 3]);
    let profile = profile.try_extend(Substitution, parameters).unwrap();
    check(profile, [0, 1, 3], [0, 1, 2]);
    let profile = profile.try_extend(Substitution, parameters).unwrap();
    check(profile, [0, 1, 2], [0, 1, 2]);
    assert_eq!(profile.to_string(), "a[0,1,2,0]b[0,1,2,0]");
}

/// With `M = 3`, a mutation must not shift entry `M` of sequence a into entry `0` of sequence b.
#[test]
fn mutation_does_not_leak_between_sequences() {
    let parameters = HistoryProfileParameters::new(20, 3);
    let mut profile = HistoryProfile::empty();
    for operation in [Match, GapInA, Match, Match, GapInA, GapInA] {
        profile = profile.try_extend(operation, parameters).unwrap();
    }
    assert_eq!(
        entries(profile, parameters),
        (vec![0, 0, 2, 3], vec![0, 1, 4, 6])
    );
    let profile = profile.try_extend(GapInB, parameters).unwrap();
    assert_eq!(
        entries(profile, parameters),
        (vec![0, 1, 1, 3], vec![0, 0, 1, 4])
    );
}

#[test]
fn match_run_limit() {
    let parameters = HistoryProfileParameters::new(5, 0);
    let mut profile = HistoryProfile::empty();
    for _ in 0..4 {
        profile = profile.try_extend(Match, parameters).unwrap();
    }
    assert!(profile.try_extend(Match, parameters).is_none());
    assert!(profile.try_extend(Substitution, parameters).is_some());
}

/// The documented range is exact, and the largest entries do not overflow.
#[test]
fn supported_range() {
    for anchor_k in 0..=u8::MAX {
        for max_anchor_mutations in 0..=u8::MAX {
            let k = u32::from(anchor_k);
            let m = u32::from(max_anchor_mutations);
            let supported = HistoryProfile::is_supported(anchor_k, max_anchor_mutations);
            assert_eq!(
                supported,
                (1..=127).contains(&k) && m.min(2 * k - 1) <= 3,
                "k = {k}, m = {m}"
            );
            if supported {
                let parameters = HistoryProfileParameters::new(anchor_k, max_anchor_mutations);
                for operation in [Match, GapInA, GapInB] {
                    let mut profile = HistoryProfile::empty();
                    for _ in 1..k {
                        profile = profile.try_extend(Match, parameters).unwrap();
                    }
                    // With `M > 0`, a gap after `k - 1` matches creates an inexact anchor-like suffix.
                    if operation == Match || m > 0 {
                        assert!(profile.try_extend(operation, parameters).is_none());
                    } else {
                        assert!(profile.try_extend(operation, parameters).is_some());
                    }
                }
            }
        }
    }
    // All supported configurations of the windowed-min-mutations lower bound.
    for anchor_k in 1..=MAX_PROFILE_ANCHOR_K {
        for max_anchor_mutations in 0..=MAX_PROFILE_INDEX.min(anchor_k - 1) {
            assert!(HistoryProfile::is_supported(anchor_k, max_anchor_mutations));
        }
    }
}

#[test]
#[should_panic(expected = "do not support")]
fn unsupported_parameters_panic() {
    HistoryProfileParameters::new(20, 4);
}

/// A small deterministic pseudo random number generator (splitmix64).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Returns true with probability `numerator / 1000`.
    fn chance(&mut self, numerator: u64) -> bool {
        self.next() % 1000 < numerator
    }

    fn operation(&mut self, match_probability: u64) -> AlignmentHistoryOperation {
        if self.chance(match_probability) {
            Match
        } else {
            OPERATIONS[1 + (self.next() % 3) as usize]
        }
    }
}

/// All supported parameters, including some with `max_anchor_mutations >= anchor_k`.
fn all_supported_parameters() -> impl Iterator<Item = (u8, u8)> {
    (1..=MAX_PROFILE_ANCHOR_K).flat_map(|anchor_k| {
        (0..=8)
            .filter(move |&max_anchor_mutations| {
                HistoryProfile::is_supported(anchor_k, max_anchor_mutations)
            })
            .map(move |max_anchor_mutations| (anchor_k, max_anchor_mutations))
    })
}

/// Follow random alignments over all supported parameters and compare the profile against its definition.
///
/// Extensions that are rejected are checked to create an anchor-like suffix and are not applied,
/// as in the searches that use the profiles.
/// Returns the reached profiles together with their entries.
fn random_walk(
    anchor_k: u8,
    max_anchor_mutations: u8,
    rng: &mut Rng,
    steps: usize,
) -> Vec<(HistoryProfile, Entries)> {
    let parameters = HistoryProfileParameters::new(anchor_k, max_anchor_mutations);
    let max_index = max_index(anchor_k, max_anchor_mutations);
    let k = usize::from(anchor_k);
    let match_probability = [500, 900, 990][(rng.next() % 3) as usize];

    let mut operations = Vec::new();
    let mut profile = HistoryProfile::empty();
    let mut reached = vec![(profile, entries(profile, parameters))];
    for _ in 0..steps {
        let operation = rng.operation(match_probability);
        operations.push(operation);
        let (len_a, len_b) = brute_force_profile(&operations, max_index);
        let expected_valid = len_a[max_index] < k && len_b[max_index] < k;

        match profile.try_extend(operation, parameters) {
            Some(extension) => {
                assert!(
                    expected_valid,
                    "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
                );
                assert_eq!(
                    entries(extension, parameters),
                    (len_a, len_b),
                    "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
                );
                profile = extension;
                reached.push((profile, entries(profile, parameters)));
            }
            None => {
                assert!(
                    !expected_valid,
                    "{operations:?} k = {anchor_k} m = {max_anchor_mutations}"
                );
                operations.pop();
            }
        }
    }
    reached
}

#[test]
fn random_alignments_match_definition() {
    let mut rng = Rng(0x05ee_d224);
    for (anchor_k, max_anchor_mutations) in all_supported_parameters() {
        for _ in 0..3 {
            random_walk(anchor_k, max_anchor_mutations, &mut rng, 300);
        }
    }
}

/// [`HistoryProfile::is_componentwise_le`] equals the componentwise comparison of the entries,
/// and the comparison is preserved by extensions.
#[test]
fn componentwise_le() {
    let mut rng = Rng(0x00d0_0123);
    for (anchor_k, max_anchor_mutations) in all_supported_parameters() {
        let parameters = HistoryProfileParameters::new(anchor_k, max_anchor_mutations);
        let reached = random_walk(anchor_k, max_anchor_mutations, &mut rng, 200);
        let mut less_or_equal_pairs = 0;
        for _ in 0..500 {
            let (p, (p_a, p_b)) = &reached[(rng.next() % reached.len() as u64) as usize];
            let (q, (q_a, q_b)) = &reached[(rng.next() % reached.len() as u64) as usize];
            let expected = p_a.iter().zip(q_a).all(|(p, q)| p <= q)
                && p_b.iter().zip(q_b).all(|(p, q)| p <= q);
            assert_eq!(
                p.is_componentwise_le(q),
                expected,
                "{p} <= {q} k = {anchor_k} m = {max_anchor_mutations}"
            );
            assert_eq!(p.is_componentwise_le(q) && q.is_componentwise_le(p), p == q);

            if expected {
                less_or_equal_pairs += 1;
                for operation in OPERATIONS {
                    if let Some(q_extension) = q.try_extend(operation, parameters) {
                        let p_extension = p.try_extend(operation, parameters).unwrap_or_else(|| {
                            panic!("{p} <= {q} {operation} k = {anchor_k} m = {max_anchor_mutations}")
                        });
                        assert!(p_extension.is_componentwise_le(&q_extension));
                    }
                }
            }
        }
        assert!(less_or_equal_pairs > 0);
    }
}

/// Entries up to `anchor_k = 127` are compared correctly and fill a byte without overflowing into the next one.
#[test]
fn largest_entries() {
    let parameters = HistoryProfileParameters::new(127, 3);
    let mut profile = HistoryProfile::empty();
    for _ in 0..125 {
        profile = profile.try_extend(Match, parameters).unwrap();
    }
    let mutated = profile.try_extend(GapInA, parameters).unwrap();
    assert_eq!(
        entries(mutated, parameters),
        (vec![0, 125, 125, 125], vec![0, 126, 126, 126])
    );
    assert!(mutated.try_extend(GapInA, parameters).is_none());
    assert!(mutated.try_extend(GapInB, parameters).is_some());

    let large = profile.try_extend(Match, parameters).unwrap();
    assert_eq!(entries(large, parameters), (vec![126; 4], vec![126; 4]));
    assert!(large.try_extend(Match, parameters).is_none());
    assert!(large.try_extend(GapInA, parameters).is_none());
    assert!(large.try_extend(GapInB, parameters).is_none());
    assert!(large.try_extend(Substitution, parameters).is_none());

    let empty = HistoryProfile::empty();
    assert!(empty.is_componentwise_le(&large));
    assert!(!large.is_componentwise_le(&empty));
    assert!(large.is_componentwise_le(&large));
    assert!(mutated.is_componentwise_le(&large));
    assert!(!large.is_componentwise_le(&mutated));
}
