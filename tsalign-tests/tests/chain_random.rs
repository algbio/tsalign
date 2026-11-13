//! Randomised end-to-end tests of tschainalign.
//!
//! The random generator is seeded deterministically, so the tests are reproducible.
//! Set `TSALIGN_E2E_SEED` to use a different base seed and `TSALIGN_E2E_RANDOM_CASES` to change the number of cases per test.
//! Failing cases are reported with their sequences and stored in `target/tsalign-e2e/failures`.
//! See [`chain_common`] for the checks.

use chain_common::{Case, run_cases};

mod chain_common;

const DEFAULT_BASE_SEED: u64 = 0x5eed_7a11_9e7e_0001;
const DEFAULT_CASES: usize = 6;

/// SplitMix64, a tiny deterministic random generator.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in `[low, high]`.
    fn range(&mut self, low: usize, high: usize) -> usize {
        assert!(low <= high);
        low + (self.next_u64() % (high - low + 1) as u64) as usize
    }

    fn chance(&mut self, probability: f64) -> bool {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 <= probability
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[self.range(0, 3)]
    }

    fn other_base(&mut self, base: u8) -> u8 {
        loop {
            let other = self.base();
            if other != base {
                return other;
            }
        }
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            c => panic!("Unsupported character {c}"),
        })
        .collect()
}

fn uniform(rng: &mut Rng, length: usize) -> Vec<u8> {
    (0..length).map(|_| rng.base()).collect()
}

/// A sequence that is repetitive in the biological sense.
fn repetitive(rng: &mut Rng, length: usize) -> Vec<u8> {
    let mut sequence = match rng.range(0, 3) {
        // Tandem repeat (micro- or minisatellite) with a few noisy copies, embedded in random sequence.
        0 => {
            let unit = {
                let length = rng.range(1, 6);
                uniform(rng, length)
            };
            let flank = rng.range(0, length / 4);
            let mut sequence = uniform(rng, flank);
            while sequence.len() < length - flank.min(length / 8) {
                sequence.extend(unit.iter().map(|&c| {
                    if rng.chance(0.03) {
                        rng.other_base(c)
                    } else {
                        c
                    }
                }));
            }
            sequence
        }
        // Interspersed repeat: a mobile-element-like segment inserted several times with divergence.
        1 => {
            let element = {
                let length = rng.range(12, 25);
                uniform(rng, length)
            };
            let mut sequence = Vec::new();
            while sequence.len() < length {
                sequence.extend({
                    let length = rng.range(0, 12);
                    uniform(rng, length)
                });
                sequence.extend(element.iter().map(|&c| {
                    if rng.chance(0.05) {
                        rng.other_base(c)
                    } else {
                        c
                    }
                }));
            }
            sequence
        }
        // Low complexity: strongly biased composition (e.g. AT-rich).
        2 => {
            let alphabet = if rng.chance(0.5) { b"AT" } else { b"GC" };
            (0..length)
                .map(|_| {
                    if rng.chance(0.85) {
                        alphabet[rng.range(0, 1)]
                    } else {
                        rng.base()
                    }
                })
                .collect()
        }
        // Mix of tandem repeats with different units.
        _ => {
            let mut sequence = Vec::new();
            while sequence.len() < length {
                let unit = {
                    let length = rng.range(2, 4);
                    uniform(rng, length)
                };
                for _ in 0..rng.range(2, 6) {
                    sequence.extend_from_slice(&unit);
                }
            }
            sequence
        }
    };

    // Fill up with random characters if a generator produced too little, and cut to length.
    while sequence.len() < length {
        sequence.push(rng.base());
    }
    sequence.truncate(length);
    sequence
}

/// Apply random substitutions, insertions and deletions.
fn mutate(rng: &mut Rng, sequence: &[u8], rate: f64) -> Vec<u8> {
    let mut result = Vec::with_capacity(sequence.len() + 8);
    for &c in sequence {
        if rng.chance(rate) {
            match rng.range(0, 3) {
                0 | 1 => result.push(rng.other_base(c)),
                2 => {
                    result.push(c);
                    for _ in 0..rng.range(1, 3) {
                        result.push(rng.base());
                    }
                }
                _ => {
                    // Deletion: skip c.
                }
            }
        } else {
            result.push(c);
        }
    }
    if result.is_empty() {
        result.push(rng.base());
    }
    result
}

/// Replace a segment of `descendant` by the reverse complement of a nearby segment of `ancestor`.
/// The ancestor may be the descendant itself (intra-sequence template switch).
fn insert_template_switch(rng: &mut Rng, ancestor: &[u8], descendant: &mut Vec<u8>) {
    let length = rng.range(10, 25);
    if descendant.len() < length + 2 || ancestor.len() < length {
        return;
    }
    let descendant_start = rng.range(1, descendant.len() - length - 1);
    // Ancestor segment around the same position, shifted by a few characters.
    let shift = rng.range(0, 30) as isize - 15;
    let ancestor_start =
        (descendant_start as isize + shift).clamp(0, (ancestor.len() - length) as isize) as usize;
    let inner = reverse_complement(&ancestor[ancestor_start..ancestor_start + length]);
    // Vary the length of the replaced descendant segment a little to also create length differences.
    let replaced_length = (length as isize + rng.range(0, 4) as isize - 2)
        .clamp(1, (descendant.len() - descendant_start) as isize)
        as usize;
    descendant.splice(descendant_start..descendant_start + replaced_length, inner);
}

fn base_seed() -> u64 {
    std::env::var("TSALIGN_E2E_SEED")
        .ok()
        .and_then(|seed| seed.parse().ok())
        .unwrap_or(DEFAULT_BASE_SEED)
}

fn case_count() -> usize {
    std::env::var("TSALIGN_E2E_RANDOM_CASES")
        .ok()
        .and_then(|cases| cases.parse().ok())
        .unwrap_or(DEFAULT_CASES)
}

/// Generate the cases of one category.
/// Each case has its own seed, such that a case can be reproduced independently of the case count.
fn generate(
    category: &str,
    category_id: u64,
    generator: impl Fn(&mut Rng) -> (Vec<u8>, Vec<u8>),
) -> Vec<Case> {
    let base_seed = base_seed();
    (0..case_count())
        .map(|index| {
            let seed = base_seed
                .wrapping_add(category_id.wrapping_mul(0x1_0000_0000))
                .wrapping_add(index as u64);
            let mut rng = Rng::new(seed);
            let (reference, query) = generator(&mut rng);
            // Randomly swap the roles, such that template switches also go from the query into the reference.
            let (reference, query) = if rng.chance(0.5) {
                (query, reference)
            } else {
                (reference, query)
            };
            Case::new(format!("{category} (seed {seed})"), reference, query)
        })
        .collect()
}

/// Length of a random sequence: mostly up to 80, sometimes up to 100.
fn random_length(rng: &mut Rng) -> usize {
    if rng.chance(0.8) {
        rng.range(20, 80)
    } else {
        rng.range(80, 100)
    }
}

fn independent(rng: &mut Rng, sequence: fn(&mut Rng, usize) -> Vec<u8>) -> (Vec<u8>, Vec<u8>) {
    let reference_length = random_length(rng);
    let reference = sequence(rng, reference_length);
    let query_length = random_length(rng);
    let query = sequence(rng, query_length);
    (reference, query)
}

fn mutated(rng: &mut Rng, sequence: fn(&mut Rng, usize) -> Vec<u8>) -> (Vec<u8>, Vec<u8>) {
    let length = random_length(rng);
    let reference = sequence(rng, length);
    let rate = [0.02, 0.03, 0.05, 0.1][rng.range(0, 3)];
    let query = mutate(rng, &reference, rate);
    (reference, query)
}

fn with_template_switches(
    rng: &mut Rng,
    sequence: fn(&mut Rng, usize) -> Vec<u8>,
) -> (Vec<u8>, Vec<u8>) {
    let length = rng.range(35, 100);
    let reference = sequence(rng, length);
    let mut query = reference.clone();
    for _ in 0..rng.range(1, 2) {
        if rng.chance(0.7) {
            insert_template_switch(rng, &reference, &mut query);
        } else {
            let ancestor = query.clone();
            insert_template_switch(rng, &ancestor, &mut query);
        }
    }
    let rate = [0.0, 0.01, 0.03][rng.range(0, 2)];
    let query = mutate(rng, &query, rate);
    (reference, query)
}

#[test]
fn uniform_independent() {
    run_cases(&generate("uniform independent", 1, |rng| {
        independent(rng, uniform)
    }));
}

#[test]
fn uniform_mutated() {
    run_cases(&generate("uniform mutated", 2, |rng| mutated(rng, uniform)));
}

#[test]
fn uniform_template_switches() {
    run_cases(&generate("uniform with template switches", 3, |rng| {
        with_template_switches(rng, uniform)
    }));
}

#[test]
fn repetitive_independent() {
    run_cases(&generate("repetitive independent", 4, |rng| {
        independent(rng, repetitive)
    }));
}

#[test]
fn repetitive_mutated() {
    run_cases(&generate("repetitive mutated", 5, |rng| {
        mutated(rng, repetitive)
    }));
}

#[test]
fn repetitive_template_switches() {
    run_cases(&generate("repetitive with template switches", 6, |rng| {
        with_template_switches(rng, repetitive)
    }));
}

#[test]
fn generators_are_deterministic_and_valid() {
    let a = generate("check", 7, |rng| with_template_switches(rng, repetitive));
    let b = generate("check", 7, |rng| with_template_switches(rng, repetitive));
    for (a, b) in a.iter().zip(&b) {
        assert_eq!(a.reference, b.reference);
        assert_eq!(a.query, b.query);
        assert!(!a.reference.is_empty() && !a.query.is_empty());
        assert!(a.max_length() <= 110, "{}", a.max_length());
        assert!(
            a.reference
                .iter()
                .chain(&a.query)
                .all(|c| b"ACGT".contains(c))
        );
    }
}
