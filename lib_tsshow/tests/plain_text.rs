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
    // The inner of a reverse TS is shown reversed, aligned to the reverse complement of the ancestor.
    assert!(output.contains("ParentR: CCTCCTCCAA"), "{output}");
    assert!(output.contains("Child2:  CCTCCTCCAA"), "{output}");
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

/// Found by the property tests: an upstream that starts with a gap in the anti-descendant, and a downstream that ends
/// with a gap in the anti-descendant, made the multipair alignment renderer panic (`multi_insert`).
#[test]
fn gaps_at_both_ends_of_the_anti_descendant() {
    let result = build_result(
        "CAG",
        "ATCCG",
        vec![
            (1, AlignmentType::PrimaryDeletion),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 3,
                    uncertainty_range: Some(TSMUncertaintyRange::new()),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (1, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 5,
                },
            ),
            (1, AlignmentType::PrimaryDeletion),
        ],
        AlignmentRange::new_complete(3, 5),
    );
    let output = render_plain_text(&result, &None);
    assert!(output.contains("ChildF:  -ATCCG-"), "{output}");
}

/// Two TSs with negative anti-descendant gaps (found by the property tests, issue #284).
/// The first has an empty downstream, so its upstream ends behind its downstream on the anti-primary.
fn multiple_negative_gaps_result() -> common::Result {
    let reverse = TemplateSwitchDirection::Reverse;
    let entrance = |first_offset, descendant, ancestor| AlignmentType::TemplateSwitchEntrance {
        first_offset,
        uncertainty_range: Some(TSMUncertaintyRange::new()),
        descendant,
        ancestor,
        direction: reverse,
    };
    build_result(
        "TGACTATGTACCCGAACTAAATCATCTAT",
        "TATTATTACGACCCGCCCGTGCGA",
        vec![
            (2, AlignmentType::PrimaryMatch),
            (2, AlignmentType::PrimaryDeletion),
            (4, AlignmentType::PrimaryMatch),
            (1, AlignmentType::PrimaryInsertion),
            (3, AlignmentType::PrimaryMatch),
            (4, AlignmentType::PrimaryInsertion),
            (4, AlignmentType::PrimaryMatch),
            (
                1,
                entrance(
                    -12,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Query,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondaryDeletion),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: -4,
                },
            ),
            (
                1,
                entrance(
                    -10,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Reference,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (2, AlignmentType::SecondarySubstitution),
            (2, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondaryDeletion),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: -4,
                },
            ),
            (8, AlignmentType::PrimaryMatch),
            (6, AlignmentType::PrimaryInsertion),
        ],
        AlignmentRange::new_complete(29, 24),
    )
}

/// The upstream (F1) of a TS may end behind the end of its downstream (F3) on the anti-primary.
/// The anti-primary must then be shown up to the end of F1 (issue #284).
#[test]
fn multiple_negative_gaps() {
    let output = render_plain_text(&multiple_negative_gaps_result(), &None);

    assert_eq!(output.matches("Switch process:").count(), 2, "{output}");
    // The first TS: F1 ends at query position 18, while F3 is empty and ends at query position 14.
    assert!(
        output.contains("Parent1: TGACTATG-TAC----CCGA\n"),
        "{output}"
    );
    assert!(
        output.contains("ChildF:  TA--TTATTACGACCCGCCC\n"),
        "{output}"
    );
    assert!(output.contains("Parent2:       -A"), "{output}");
    // The second TS: F3 starts at query position 10, before F1 ends at query position 14.
    assert!(output.contains("Parent3: TCATCTAT------\n"), "{output}");
    assert!(output.contains("ChildF:  ACCCGCCCGTGCGA\n"), "{output}");
}

/// The no-ts alignment is shown over the same anti-primary interval as the TS (issue #284).
#[test]
fn multiple_negative_gaps_with_no_ts_alignment() {
    let no_ts_result = build_result(
        "TGACTATGTACCCGAACTAAATCATCTAT",
        "TATTATTACGACCCGCCCGTGCGA",
        vec![
            (24, AlignmentType::PrimaryMatch),
            (5, AlignmentType::PrimaryDeletion),
        ],
        AlignmentRange::new_complete(29, 24),
    );
    let output = render_plain_text(&multiple_negative_gaps_result(), &Some(no_ts_result));

    assert_eq!(output.matches("No-ts alignment:").count(), 2, "{output}");
    // The first TS covers the query from 0 to 18, the second from 10 to 24.
    assert!(output.contains("Child:  TATTATTACGACCCGCCC\n"), "{output}");
    assert!(output.contains("Child:  ACCCGCCCGTGCGA\n"), "{output}");
}
