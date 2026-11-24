//! Tests of the plain text output (`tsalign show -t`, `Alignment.viz_template_switches()` in Python).

use common::{Fixture, build_result, fixture_names, render_plain_text};
use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_geometry::AlignmentRange,
    alignment_result::{AlignmentResult, AlignmentType, TSMUncertaintyRange},
};
use lib_tsshow::plain_text::show_template_switches;

mod common;

/// `twin_10_ts_only.fa`: the query is the reverse complement of the reference.
fn whole_sequence_reverse_ts() -> common::Result {
    build_result(
        "GGAGGAGGTT",
        "AACCTCCTCC",
        vec![
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    uncertainty_range: Some(TSMUncertaintyRange::new()),
                    descendant: TemplateSwitchDescendant::Query,
                    ancestor: TemplateSwitchAncestor::Reference,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (10, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 10,
                },
            ),
        ],
        AlignmentRange::new_complete(10, 10),
    )
}

/// The complete output goes to the given writer, including the headers of each template switch.
#[test]
fn output_goes_to_the_writer() {
    let fixture = Fixture::load("classic_both_descendants");
    let output = render_plain_text(&fixture.result, &fixture.no_ts_result);

    assert_eq!(output.matches("Switch process:").count(), 2, "{output}");
    assert_eq!(output.matches("No-ts alignment:").count(), 2, "{output}");
    assert_eq!(output.matches("Direction: reverse").count(), 2, "{output}");
    assert!(output.contains("Parent: "), "{output}");
    assert!(output.contains("Child: "), "{output}");
}

/// Each template switch is shown once, with its direction.
#[test]
fn every_template_switch_is_shown() {
    for name in fixture_names() {
        let fixture = Fixture::load(&name);
        let AlignmentResult::WithTarget { alignment, .. } = &fixture.result else {
            panic!("{name} has no alignment");
        };
        let directions: Vec<_> = alignment
            .iter_compact()
            .filter_map(|(_, alignment_type)| match alignment_type {
                AlignmentType::TemplateSwitchEntrance { direction, .. } => Some(*direction),
                _ => None,
            })
            .collect();
        let output = std::panic::catch_unwind(|| render_plain_text(&fixture.result, &None))
            .unwrap_or_else(|_| panic!("{name}: rendering panicked"));

        assert_eq!(
            output.matches("Switch process:").count(),
            directions.len(),
            "{name}:\n{output}"
        );
        let forward = directions
            .iter()
            .filter(|d| **d == TemplateSwitchDirection::Forward)
            .count();
        assert_eq!(
            output.matches("Direction: forward").count(),
            forward,
            "{name}:\n{output}"
        );
        assert_eq!(
            output.matches("Direction: reverse").count(),
            directions.len() - forward,
            "{name}:\n{output}"
        );
    }
}

/// A TS that spans both sequences completely (the Python smoke test of `viz_template_switches`).
#[test]
fn whole_sequence_template_switch() {
    let output = render_plain_text(&whole_sequence_reverse_ts(), &None);
    assert!(output.contains("Switch process:"), "{output}");
    assert!(output.contains("AACCTCCTCC"), "{output}");
}

/// An aborted alignment has nothing to show, which is not an error.
#[test]
fn alignment_without_target_shows_nothing() {
    let fixture = Fixture::load_test_files(
        "twin_ari_chrX_146823507_146823598",
        "twin_ari_chrX_146823507_146823598_no_ts.toml",
        "twin_ari_chrX_146823507_146823598_no_ts.toml",
    );
    assert!(matches!(
        fixture.result,
        AlignmentResult::WithoutTarget { .. }
    ));
    assert_eq!(render_plain_text(&fixture.result, &None), "");
}

/// An alignment without template switches has nothing to show.
#[test]
fn alignment_without_template_switches_shows_nothing() {
    let fixture = Fixture::load("classic_qr_reverse");
    let no_ts_result = fixture.no_ts_result.unwrap();
    assert_eq!(render_plain_text(&no_ts_result, &None), "");
}

/// A no-TS alignment without target is skipped.
#[test]
fn no_ts_alignment_without_target_is_skipped() {
    let fixture = Fixture::load_test_files(
        "twin_ari_chrX_146823507_146823598",
        "twin_ari_chrX_146823507_146823598.toml",
        "twin_ari_chrX_146823507_146823598_no_ts.toml",
    );
    let output = render_plain_text(&fixture.result, &fixture.no_ts_result);
    assert_eq!(output.matches("Switch process:").count(), 3, "{output}");
    assert!(!output.contains("No-ts alignment:"), "{output}");
}

/// A "no-TS" alignment that contains template switches is rejected with an error instead of a panic.
#[test]
fn no_ts_alignment_with_template_switches_is_an_error() {
    let fixture = Fixture::load("classic_qr_reverse");
    let mut output = Vec::new();
    let result = std::panic::catch_unwind(move || {
        show_template_switches(&mut output, &fixture.result, &Some(fixture.result.clone()))
    })
    .expect("show_template_switches panicked");
    let error = result.expect_err("the no-TS alignment with template switches was accepted");
    assert!(error.to_string().contains("template switches"), "{error}");
}

/// A no-TS alignment of other sequences is rejected with an error instead of a panic or a wrong output.
#[test]
fn no_ts_alignment_of_other_sequences_is_an_error() {
    let fixture = Fixture::load("classic_qr_reverse");
    let other = Fixture::load("classic_rq_reverse");
    let mut output = Vec::new();
    let result = std::panic::catch_unwind(move || {
        show_template_switches(&mut output, &fixture.result, &other.no_ts_result)
    })
    .expect("show_template_switches panicked");
    let error = result.expect_err("the no-TS alignment of other sequences was accepted");
    assert!(error.to_string().contains("sequences"), "{error}");
}
