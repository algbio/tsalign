//! Property tests: random small alignments with template switches must render without panics,
//! and the TS arrangement must keep its invariants.
//!
//! The random alignments are deterministic. Set `TSSHOW_FUZZ_SEED` to use another seed
//! and `TSSHOW_FUZZ_CASES` to change the number of cases per test.

use std::{env, panic};

use common::{
    Random, RandomAlignment, random_alignment, render_plain_text, render_svg, svg_config,
};
use lib_tsalign::alignment_result::AlignmentResult;
use lib_tsshow::{
    svg::UncertaintyRangeMode,
    ts_arrangement::{TsArrangement, character::Char, source::SourceChar},
};

mod common;

thread_local! {
    static PANIC_LOCATION: std::cell::Cell<Option<String>> = const { std::cell::Cell::new(None) };
}

fn seed() -> u64 {
    env::var("TSSHOW_FUZZ_SEED")
        .map(|seed| seed.parse().expect("TSSHOW_FUZZ_SEED is a number"))
        .unwrap_or(279)
}

fn cases() -> usize {
    env::var("TSSHOW_FUZZ_CASES")
        .map(|cases| cases.parse().expect("TSSHOW_FUZZ_CASES is a number"))
        .unwrap_or(200)
}

/// Runs `check` on random alignments, and reports the smallest failing alignments.
/// Known bugs that are not fixed yet: panics whose location contains the first string
/// and whose message contains the second string are reported, but do not fail the tests.
///
/// They are reproduced by ignored tests in `plain_text.rs` and `svg_regressions.rs`.
/// Remove the entries when fixing the bugs.
const KNOWN_BUGS: &[(&str, &str)] = &[
    // Issue #284. Plain text: an inner or the downstream is aligned beyond the anti-primary range.
    (
        "mutlipair_alignment_renderer.rs",
        "translated_reference_sequence.len()",
    ),
    ("mutlipair_alignment_renderer.rs", "index out of bounds"),
    // Issue #285. SVG: overlapping forward TSs.
    ("ts_arrangement/source.rs", "Blank has no copy property"),
];

fn is_known_bug(message: &str) -> bool {
    let location = message.split(": ").next().unwrap();
    KNOWN_BUGS
        .iter()
        .any(|(file, text)| location.contains(file) && message.contains(text))
}

fn fuzz(
    name: &str,
    forward: bool,
    uncertainty: bool,
    ranges: bool,
    check: impl Fn(&RandomAlignment) + panic::RefUnwindSafe,
) {
    let mut random = Random::new(seed() ^ name.len() as u64);
    let mut failures = Vec::new();

    // The panic hook is global, so the fuzz tests must not run concurrently.
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _lock = LOCK.lock().unwrap_or_else(|poisoned| poisoned.into_inner());

    // Silence the panic messages of the failing cases, they are reported below with the location of the panic.
    let hook = panic::take_hook();
    panic::set_hook(Box::new(|info| {
        PANIC_LOCATION.set(info.location().map(ToString::to_string));
    }));
    for _ in 0..cases() {
        let alignment = random_alignment(&mut random, forward, uncertainty, ranges);
        if let Err(error) = panic::catch_unwind(|| check(&alignment)) {
            let message = error
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| error.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default();
            let location = PANIC_LOCATION.take().unwrap_or_default();
            failures.push((
                alignment,
                format!("{location}: {}", message.replace('\n', " ")),
            ));
        }
    }
    panic::set_hook(hook);

    failures.sort_by_key(|(alignment, _)| {
        alignment.reference.len() + alignment.query.len() + alignment.alignment.len()
    });
    // Group the failures by the location of the panic, and report the smallest failing alignment of each group.
    let mut groups: Vec<(&str, &RandomAlignment, &str)> = Vec::new();
    for (alignment, message) in &failures {
        let location = message.split(": ").next().unwrap();
        if !groups.iter().any(|(l, _, _)| *l == location) {
            groups.push((location, alignment, message));
        }
    }
    let (known, groups): (Vec<_>, Vec<_>) = groups
        .into_iter()
        .partition(|(_, _, message)| is_known_bug(message));
    for (_, alignment, message) in &known {
        println!(
            "Known bug: {message}\nSmallest failing alignment:\n{}",
            alignment.reproduction()
        );
    }
    let unknown_failures = failures
        .iter()
        .filter(|(_, message)| !is_known_bug(message))
        .count();

    let report: Vec<String> = groups
        .iter()
        .map(|(_, alignment, message)| {
            format!(
                "{message}\nSmallest failing alignment:\n{}",
                alignment.reproduction()
            )
        })
        .collect();
    assert!(
        unknown_failures == 0,
        "{name}: {} of {} random alignments failed (seed {}).\n\n{}",
        unknown_failures,
        cases(),
        seed(),
        report.join("\n\n"),
    );
}

#[test]
fn random_reverse_alignments_render_as_plain_text() {
    fuzz("plain_text", false, false, false, |alignment| {
        render_plain_text(&alignment.result(), &None);
    });
}

#[test]
fn random_reverse_alignments_render_as_svg() {
    fuzz("svg", false, false, false, |alignment| {
        render_svg(
            &alignment.result(),
            &None,
            &svg_config(true, false, None, UncertaintyRangeMode::None),
        )
        .unwrap();
    });
}

#[test]
fn random_reverse_alignments_with_options_render_as_svg() {
    fuzz("svg_options", false, true, false, |alignment| {
        let result = alignment.result();
        for config in [
            svg_config(true, true, Some(2), UncertaintyRangeMode::Full),
            svg_config(false, false, None, UncertaintyRangeMode::InnerOnly),
            svg_config(false, false, Some(0), UncertaintyRangeMode::None),
        ] {
            render_svg(&result, &None, &config).unwrap();
        }
    });
}

#[test]
fn random_alignments_with_ranges_render() {
    fuzz("ranges", false, false, true, |alignment| {
        let result = alignment.result();
        render_plain_text(&result, &None);
        render_svg(
            &result,
            &None,
            &svg_config(true, false, None, UncertaintyRangeMode::None),
        )
        .unwrap();
    });
}

#[test]
fn random_forward_alignments_render() {
    fuzz("forward", true, false, false, |alignment| {
        let result = alignment.result();
        render_plain_text(&result, &None);
        render_svg(
            &result,
            &None,
            &svg_config(true, false, None, UncertaintyRangeMode::None),
        )
        .unwrap();
    });
}

/// The source rows of the arrangement contain each character of their sequence exactly once (not counting copies),
/// in order, and all rows have the same width.
#[test]
fn random_arrangements_keep_their_invariants() {
    fuzz("arrangement", false, true, true, |alignment| {
        let AlignmentResult::WithTarget {
            alignment: alignment_steps,
            statistics,
        } = alignment.result()
        else {
            unreachable!()
        };
        for uncertainty_range_mode in [
            UncertaintyRangeMode::None,
            UncertaintyRangeMode::InnerOnly,
            UncertaintyRangeMode::Full,
        ] {
            let mut arrangement = TsArrangement::new(
                statistics.reference_offset,
                statistics.query_offset,
                alignment.reference.len(),
                alignment.query.len(),
                alignment_steps.iter_flat_cloned(),
                uncertainty_range_mode,
            )
            .unwrap();
            let visible = check_arrangement(&arrangement, alignment, true);
            arrangement.show_complete_complements_if_used();
            arrangement.remove_empty_columns();
            assert_eq!(
                check_arrangement(&arrangement, alignment, false),
                visible,
                "removing empty columns removed visible characters"
            );
        }
    });
}

/// Checks the invariants of the arrangement and returns the visible source columns of the reference and the query.
///
/// If `complete` is false, then hidden characters may have been removed.
fn check_arrangement(
    arrangement: &TsArrangement,
    alignment: &RandomAlignment,
    complete: bool,
) -> Vec<Vec<usize>> {
    let width = arrangement.width();
    assert_eq!(arrangement.reference().len(), width);
    assert_eq!(arrangement.query().len(), width);
    assert_eq!(arrangement.reference_complement().len(), width);
    assert_eq!(arrangement.query_complement().len(), width);
    for inner in arrangement.inners().iter_values() {
        assert_eq!(inner.sequence().len(), width, "inner has a different width");
        // The inner of a TS is a substring of its descendant, also if it consists of copies (#286).
        let inner_source_columns: Vec<usize> = inner
            .template_switch()
            .inner
            .iter()
            .map(|c| c.source_column().primitive())
            .collect();
        assert!(
            inner_source_columns.windows(2).all(|w| w[0] + 1 == w[1]),
            "the inner of a TS is not a substring of its descendant: {inner_source_columns:?}"
        );
    }

    let mut visible = Vec::new();
    for (row, length) in [
        (arrangement.reference(), alignment.reference.len()),
        (arrangement.query(), alignment.query.len()),
    ] {
        let source_columns: Vec<usize> = row
            .iter_values()
            .filter(|c| c.is_source_char())
            .map(|c| c.source_column().primitive())
            .collect();
        if complete {
            assert_eq!(
                source_columns,
                (0..length).collect::<Vec<_>>(),
                "the source row does not contain each source character once, in order"
            );
        } else {
            assert!(
                source_columns.is_sorted_by(|a, b| a < b)
                    && source_columns.iter().all(|c| *c < length),
                "the source row contains source characters out of order: {source_columns:?}"
            );
        }
        assert!(
            row.iter_values()
                .filter(|c| matches!(c, SourceChar::Separator))
                .count()
                <= 2,
            "more than two separators"
        );
        visible.push(
            row.iter_values()
                .filter(|c| c.is_source_char() && !c.is_hidden())
                .map(|c| c.source_column().primitive())
                .collect(),
        );
    }
    visible
}
