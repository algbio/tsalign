//! Shared helpers of the integration tests of lib_tsshow.
//!
//! - Loading the alignment fixtures from `tests/fixtures` (see `tests/fixtures/generate.sh`).
//! - Snapshot (golden file) comparison against `tests/snapshots`; set `TSSHOW_BLESS=1` to (re)write the snapshots.
//! - Building alignment results from sequences and alignments, and generating random alignments.
//!
//! Not every test binary uses every helper.
#![allow(dead_code)]

use std::{
    env, fs,
    path::{Path, PathBuf},
};

use compact_genome::implementation::{
    alphabets::dna_alphabet_or_n::DnaAlphabetOrN, vec_sequence::VectorGenome,
};
use compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use generic_a_star::AStarResult;
use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_geometry::AlignmentRange,
    alignment_result::{AlignmentResult, AlignmentType, TSMUncertaintyRange},
    costs::U64Cost,
};
use lib_tsshow::{
    plain_text::show_template_switches,
    svg::{SvgConfig, UncertaintyRangeMode, create_ts_svg},
};

pub type Result = AlignmentResult<U64Cost>;

pub fn crate_directory() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn repository_directory() -> PathBuf {
    crate_directory().parent().unwrap().to_owned()
}

pub fn fixture_directory() -> PathBuf {
    crate_directory().join("tests").join("fixtures")
}

pub fn snapshot_directory() -> PathBuf {
    crate_directory().join("tests").join("snapshots")
}

/// Loads a toml output of `tsalign align`.
pub fn load_result(path: impl AsRef<Path>) -> Result {
    let path = path.as_ref();
    let text = fs::read_to_string(path).unwrap_or_else(|error| panic!("reading {path:?}: {error}"));
    toml::from_str(&text).unwrap_or_else(|error| panic!("parsing {path:?}: {error}"))
}

/// An alignment with template switches together with the corresponding alignment without template switches.
pub struct Fixture {
    pub name: String,
    pub result: Result,
    pub no_ts_result: Option<Result>,
}

impl Fixture {
    /// Loads `<name>.toml` and `<name>.no_ts.toml` from the fixture directory.
    pub fn load(name: &str) -> Self {
        let directory = fixture_directory();
        let no_ts_path = directory.join(format!("{name}.no_ts.toml"));
        Self {
            name: name.to_owned(),
            result: load_result(directory.join(format!("{name}.toml"))),
            no_ts_result: no_ts_path.exists().then(|| load_result(no_ts_path)),
        }
    }

    /// Loads a pair of alignments from `test_files` in the repository root.
    pub fn load_test_files(name: &str, ts_file: &str, no_ts_file: &str) -> Self {
        let directory = repository_directory().join("test_files");
        Self {
            name: name.to_owned(),
            result: load_result(directory.join(ts_file)),
            no_ts_result: Some(load_result(directory.join(no_ts_file))),
        }
    }
}

/// The names of all alignment fixtures (without the `.no_ts` ones).
pub fn fixture_names() -> Vec<String> {
    let mut names: Vec<_> = fs::read_dir(fixture_directory())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter_map(|name| name.strip_suffix(".toml").map(str::to_owned))
        .filter(|name| !name.ends_with(".no_ts"))
        .collect();
    names.sort();
    names
}

/// The alignments from `test_files` that the old show-bug scripts and PRs #219/#220 were about.
///
/// Returns the name, the TS alignment and the no-TS alignment.
pub fn test_files_fixtures() -> Vec<Fixture> {
    [
        ("ari_200", "ari_200.toml", "ari_200_nots.toml"),
        ("ari_868", "ari_868.toml", "ari_868_nots.toml"),
        ("underscore", "underscore.toml", "underscore_no_ts.toml"),
        (
            "twin_heli_MDC1-AS1_10",
            "twin_heli_MDC1-AS1_10.toml",
            "twin_heli_MDC1-AS1_10_no_ts.toml",
        ),
        (
            "twin_heli_linc01237_114",
            "twin_heli_linc01237_114.toml",
            "twin_heli_linc01237_114_no_ts.toml",
        ),
    ]
    .into_iter()
    .map(|(name, ts, no_ts)| Fixture::load_test_files(name, ts, no_ts))
    .collect()
}

pub fn render_plain_text(result: &Result, no_ts_result: &Option<Result>) -> String {
    let mut output = Vec::new();
    show_template_switches(&mut output, result, no_ts_result).unwrap();
    String::from_utf8(output).unwrap()
}

pub fn render_svg(
    result: &Result,
    no_ts_result: &Option<Result>,
    config: &SvgConfig,
) -> lib_tsshow::error::Result<String> {
    let mut output = Vec::new();
    create_ts_svg(&mut output, result, no_ts_result, config)?;
    Ok(String::from_utf8(output).unwrap())
}

pub fn svg_config(
    render_arrows: bool,
    render_more_complement: bool,
    restrict_context: Option<usize>,
    uncertainty_range_mode: UncertaintyRangeMode,
) -> SvgConfig {
    SvgConfig {
        render_arrows,
        render_more_complement,
        restrict_context,
        uncertainty_range_mode,
    }
}

pub fn default_svg_config() -> SvgConfig {
    svg_config(false, false, None, UncertaintyRangeMode::None)
}

/// Replaces the content of the `<style>` element (the embedded fonts) by a short placeholder.
///
/// The embedded fonts are large, and their bytes depend on the font subsetting library.
/// They are tested separately.
pub fn strip_embedded_fonts(svg: &str) -> String {
    let start = svg.find("<style").expect("SVG has a style element");
    let content_start = start + svg[start..].find('>').unwrap() + 1;
    let content_end = content_start + svg[content_start..].find("</style>").unwrap();
    let families: Vec<_> = svg[content_start..content_end]
        .match_indices("font-family: \"")
        .map(|(index, prefix)| {
            let family_start = content_start + index + prefix.len();
            let family_end = family_start + svg[family_start..].find('"').unwrap();
            &svg[family_start..family_end]
        })
        .collect();
    format!(
        "{}[embedded fonts: {}]{}",
        &svg[..content_start],
        families.join(", "),
        &svg[content_end..]
    )
}

/// Compares `actual` against the snapshot `tests/snapshots/<name>`.
///
/// If the environment variable `TSSHOW_BLESS` is set to `1`, then the snapshot is (re)written instead.
/// On a mismatch, the actual output is written to the test's temporary directory for inspection.
pub fn assert_snapshot(name: &str, actual: &str) {
    let path = snapshot_directory().join(name);

    if env::var("TSSHOW_BLESS").is_ok_and(|value| value == "1") {
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, actual).unwrap();
        return;
    }

    let expected = fs::read_to_string(&path).unwrap_or_else(|error| {
        panic!(
            "Snapshot {path:?} cannot be read ({error}).\n\
             Create it with `TSSHOW_BLESS=1 cargo test -p lib_tsshow` and review it."
        )
    });

    if expected != actual {
        let actual_path = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
            .join("tsshow-snapshots")
            .join(name);
        fs::create_dir_all(actual_path.parent().unwrap()).unwrap();
        fs::write(&actual_path, actual).unwrap();

        let (line_number, (expected_line, actual_line)) = expected
            .lines()
            .chain(std::iter::repeat("<end of file>"))
            .zip(actual.lines().chain(std::iter::repeat("<end of file>")))
            .enumerate()
            .find(|(_, (expected, actual))| expected != actual)
            .unwrap();
        panic!(
            "Snapshot {path:?} differs from the actual output, written to {actual_path:?}.\n\
             First difference in line {}:\n  expected: {expected_line}\n  actual:   {actual_line}\n\
             If the change is intended, re-bless the snapshots with `TSSHOW_BLESS=1 cargo test -p lib_tsshow` \
             and review the diff.",
            line_number + 1,
        );
    }
}

/// Builds an alignment result from two DNA sequences and an alignment of them.
pub fn build_result(
    reference: &str,
    query: &str,
    alignment: Vec<(usize, AlignmentType)>,
    range: AlignmentRange,
) -> Result {
    let reference_genome = VectorGenome::<DnaAlphabetOrN>::from_slice_u8(reference.as_bytes())
        .expect("reference is DNA");
    let query_genome =
        VectorGenome::<DnaAlphabetOrN>::from_slice_u8(query.as_bytes()).expect("query is DNA");
    AlignmentResult::new_with_target(
        alignment,
        reference_genome.as_genome_subsequence(),
        query_genome.as_genome_subsequence(),
        "reference",
        "query",
        range,
        AStarResult::FoundTarget {
            identifier: (),
            cost: U64Cost::from(0u64),
        },
        0.0,
        0,
        0,
        0,
        reference.len(),
        query.len(),
    )
}

/// A small deterministic pseudo random number generator (xorshift64*), so that the tests need no extra dependency.
pub struct Random(u64);

impl Random {
    pub fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// A uniformly random number in `0..limit`.
    pub fn below(&mut self, limit: usize) -> usize {
        assert!(limit > 0);
        (self.next_u64() % limit as u64) as usize
    }

    /// A uniformly random number in `range`.
    pub fn range(&mut self, range: std::ops::RangeInclusive<usize>) -> usize {
        range.start() + self.below(range.end() - range.start() + 1)
    }

    pub fn chance(&mut self, numerator: usize, denominator: usize) -> bool {
        self.below(denominator) < numerator
    }

    pub fn dna(&mut self, length: usize) -> String {
        (0..length)
            .map(|_| b"ACGT"[self.below(4)] as char)
            .collect()
    }
}

/// A random alignment with template switches, see [`random_alignment`].
#[derive(Debug, Clone)]
pub struct RandomAlignment {
    pub reference: String,
    pub query: String,
    pub alignment: Vec<(usize, AlignmentType)>,
    pub range: AlignmentRange,
}

impl RandomAlignment {
    pub fn result(&self) -> Result {
        build_result(
            &self.reference,
            &self.query,
            self.alignment.clone(),
            self.range.clone(),
        )
    }

    /// A rust expression that rebuilds this alignment, for reproducing failures in a unit test.
    pub fn reproduction(&self) -> String {
        format!(
            "reference: {:?}\nquery: {:?}\nrange: {:?}\nalignment: {:?}",
            self.reference, self.query, self.range, self.alignment
        )
    }
}

/// Generates a random valid alignment with at least one template switch.
///
/// The alignment respects the invariants of alignments created by the aligners:
/// the primary alignment advances through `offset..limit` of both sequences,
/// every template switch inner aligns at least one character of the descendant
/// to a character of the ancestor that lies within the ancestor,
/// and the anti-descendant gap of a template switch does not leave the anti-descendant.
///
/// `forward` allows forward template switches, `uncertainty` adds random TSM uncertainty ranges,
/// and `ranges` allows the alignment to start and end within the sequences.
pub fn random_alignment(
    random: &mut Random,
    forward: bool,
    uncertainty: bool,
    ranges: bool,
) -> RandomAlignment {
    loop {
        if let Some(alignment) = try_random_alignment(random, forward, uncertainty, ranges) {
            return alignment;
        }
    }
}

/// Limits the TSM uncertainty ranges to what an aligner could report:
/// the switch points can only be shifted over matches, and the inner must keep at least one character.
fn clamp_uncertainty_ranges(alignment: &[(usize, AlignmentType)]) -> Vec<(usize, AlignmentType)> {
    let mut flat: Vec<AlignmentType> = alignment
        .iter()
        .flat_map(|(multiplicity, alignment_type)| {
            std::iter::repeat_n(*alignment_type, *multiplicity)
        })
        .collect();
    let count_run = |range: &mut dyn Iterator<Item = &AlignmentType>, alignment_type| {
        range.take_while(|t| **t == alignment_type).count() as u16
    };

    for index in 0..flat.len() {
        let AlignmentType::TemplateSwitchEntrance {
            uncertainty_range: Some(mut range),
            first_offset,
            descendant,
            ancestor,
            direction,
        } = flat[index]
        else {
            continue;
        };
        let exit = index
            + flat[index..]
                .iter()
                .position(|t| matches!(t, AlignmentType::TemplateSwitchExit { .. }))
                .unwrap();
        let inner = &flat[index + 1..exit];
        let inner_length = inner
            .iter()
            .filter(|t| **t != AlignmentType::SecondaryDeletion)
            .count() as u16;

        range.start_left_shift = range.start_left_shift.min(count_run(
            &mut flat[..index].iter().rev(),
            AlignmentType::PrimaryMatch,
        ));
        range.end_right_shift = range.end_right_shift.min(count_run(
            &mut flat[exit + 1..].iter(),
            AlignmentType::PrimaryMatch,
        ));
        range.start_right_shift = range
            .start_right_shift
            .min(count_run(&mut inner.iter(), AlignmentType::SecondaryMatch));
        range.end_left_shift = range.end_left_shift.min(count_run(
            &mut inner.iter().rev(),
            AlignmentType::SecondaryMatch,
        ));
        while range.start_right_shift + range.end_left_shift >= inner_length
            && range.start_right_shift + range.end_left_shift > 0
        {
            if range.start_right_shift > 0 {
                range.start_right_shift -= 1;
            } else {
                range.end_left_shift -= 1;
            }
        }

        flat[index] = AlignmentType::TemplateSwitchEntrance {
            uncertainty_range: Some(range),
            first_offset,
            descendant,
            ancestor,
            direction,
        };
    }

    let mut result = Vec::new();
    for alignment_type in flat {
        push(&mut result, alignment_type);
    }
    result
}

fn push(alignment: &mut Vec<(usize, AlignmentType)>, alignment_type: AlignmentType) {
    match alignment.last_mut() {
        Some((multiplicity, last)) if *last == alignment_type => *multiplicity += 1,
        _ => alignment.push((1, alignment_type)),
    }
}

fn try_random_alignment(
    random: &mut Random,
    forward: bool,
    uncertainty: bool,
    ranges: bool,
) -> Option<RandomAlignment> {
    let reference_length = random.range(4..=30);
    let reference = random.dna(reference_length);
    let query_length = random.range(4..=30);
    let query = random.dna(query_length);
    let (reference_offset, query_offset) = if ranges && random.chance(1, 3) {
        (random.range(0..=3), random.range(0..=3))
    } else {
        (0, 0)
    };
    let lengths = [reference.len(), query.len()];
    let mut position = [reference_offset, query_offset];
    let mut alignment: Vec<(usize, AlignmentType)> = Vec::new();
    let mut template_switch_count = 0;
    let steps = random.range(1..=20);

    for _ in 0..steps {
        let remaining = [lengths[0] - position[0], lengths[1] - position[1]];
        let choice = random.below(10);

        if choice < 3 && remaining.iter().all(|r| *r > 0) {
            // Template switch.
            let descendant_index = random.below(2);
            let ancestor_index = random.below(2);
            let anti_descendant_index = 1 - descendant_index;
            let direction = if forward && random.chance(1, 3) {
                TemplateSwitchDirection::Forward
            } else {
                TemplateSwitchDirection::Reverse
            };

            // The inner alignment: at least one match, and the descendant must not leave its sequence.
            let mut inner = Vec::new();
            let mut descendant_length = 0;
            let mut ancestor_length = 0;
            let inner_steps = random.range(1..=6);
            for step in 0..inner_steps {
                let alignment_type = if step == 0 {
                    AlignmentType::SecondaryMatch
                } else {
                    match random.below(6) {
                        0 => AlignmentType::SecondaryInsertion,
                        1 => AlignmentType::SecondaryDeletion,
                        2 => AlignmentType::SecondarySubstitution,
                        _ => AlignmentType::SecondaryMatch,
                    }
                };
                let consumes_descendant = alignment_type != AlignmentType::SecondaryDeletion;
                let consumes_ancestor = alignment_type != AlignmentType::SecondaryInsertion;
                if consumes_descendant && descendant_length + 1 > remaining[descendant_index] {
                    break;
                }
                if consumes_ancestor && ancestor_length + 1 > lengths[ancestor_index] {
                    break;
                }
                descendant_length += usize::from(consumes_descendant);
                ancestor_length += usize::from(consumes_ancestor);
                inner.push(alignment_type);
            }

            // SP2 is the first ancestor position of the inner. Forward inners go right, reverse inners go left.
            let sp2 = match direction {
                TemplateSwitchDirection::Forward => {
                    random.range(0..=lengths[ancestor_index] - ancestor_length)
                }
                TemplateSwitchDirection::Reverse => {
                    random.range(ancestor_length..=lengths[ancestor_index])
                }
            };
            let first_offset = sp2 as isize - position[ancestor_index] as isize;

            // The anti-descendant gap, i.e. SP4 - SP1 on the anti-descendant.
            // It may repeat characters (negative gap), but not beyond the start of the alignment.
            let anti_descendant_offset = [reference_offset, query_offset][anti_descendant_index];
            let anti_descendant_gap = random.range(
                anti_descendant_offset
                    ..=position[anti_descendant_index] + remaining[anti_descendant_index],
            ) as isize
                - position[anti_descendant_index] as isize;
            // Keep negative gaps small most of the time, like in real alignments.
            let anti_descendant_gap = if anti_descendant_gap < -6 && random.chance(3, 4) {
                anti_descendant_gap / 3
            } else {
                anti_descendant_gap
            };

            let uncertainty_range = if uncertainty && random.chance(1, 2) {
                Some(TSMUncertaintyRange {
                    start_left_shift: random.range(0..=2) as u16,
                    start_right_shift: random.range(0..=descendant_length.min(2)) as u16,
                    end_left_shift: random.range(0..=descendant_length.min(2)) as u16,
                    end_right_shift: random.range(0..=2) as u16,
                })
            } else if uncertainty {
                None
            } else {
                Some(TSMUncertaintyRange::new())
            };

            push(
                &mut alignment,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset,
                    uncertainty_range,
                    descendant: if descendant_index == 0 {
                        TemplateSwitchDescendant::Reference
                    } else {
                        TemplateSwitchDescendant::Query
                    },
                    ancestor: if ancestor_index == 0 {
                        TemplateSwitchAncestor::Reference
                    } else {
                        TemplateSwitchAncestor::Query
                    },
                    direction,
                },
            );
            for alignment_type in inner {
                push(&mut alignment, alignment_type);
            }
            push(
                &mut alignment,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap,
                },
            );

            position[descendant_index] += descendant_length;
            position[anti_descendant_index] =
                (position[anti_descendant_index] as isize + anti_descendant_gap) as usize;
            template_switch_count += 1;
        } else {
            let alignment_type = match choice {
                3 if remaining[1] > 0 => AlignmentType::PrimaryInsertion,
                4 if remaining[0] > 0 => AlignmentType::PrimaryDeletion,
                5 if remaining[0] > 0 && remaining[1] > 0 => AlignmentType::PrimarySubstitution,
                _ if remaining[0] > 0 && remaining[1] > 0 => AlignmentType::PrimaryMatch,
                _ => continue,
            };
            let length = random.range(1..=4);
            for _ in 0..length {
                let remaining = [lengths[0] - position[0], lengths[1] - position[1]];
                let (consumes_reference, consumes_query) = match alignment_type {
                    AlignmentType::PrimaryInsertion => (false, true),
                    AlignmentType::PrimaryDeletion => (true, false),
                    _ => (true, true),
                };
                if (consumes_reference && remaining[0] == 0)
                    || (consumes_query && remaining[1] == 0)
                {
                    break;
                }
                position[0] += usize::from(consumes_reference);
                position[1] += usize::from(consumes_query);
                push(&mut alignment, alignment_type);
            }
        }
    }

    if template_switch_count == 0 {
        return None;
    }

    // Without ranges, the alignment must cover the complete sequences.
    if !ranges {
        while position[0] < lengths[0] || position[1] < lengths[1] {
            let alignment_type = if position[0] < lengths[0] && position[1] < lengths[1] {
                AlignmentType::PrimaryMatch
            } else if position[0] < lengths[0] {
                AlignmentType::PrimaryDeletion
            } else {
                AlignmentType::PrimaryInsertion
            };
            position[0] += usize::from(alignment_type != AlignmentType::PrimaryInsertion);
            position[1] += usize::from(alignment_type != AlignmentType::PrimaryDeletion);
            push(&mut alignment, alignment_type);
        }
    }

    if uncertainty {
        alignment = clamp_uncertainty_ranges(&alignment);
    }

    Some(RandomAlignment {
        reference,
        query,
        alignment,
        range: AlignmentRange::new(reference_offset, query_offset, position[0], position[1]),
    })
}
