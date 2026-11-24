//! Regression tests of the SVG output for bugs found by the property tests (`fuzz.rs`).

use common::{build_result, default_svg_config, render_svg};
use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_geometry::AlignmentRange,
    alignment_result::{AlignmentType, TSMUncertaintyRange},
};
use lib_tsshow::svg::UncertaintyRangeMode;

mod common;

fn entrance(
    first_offset: isize,
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
    direction: TemplateSwitchDirection,
) -> AlignmentType {
    AlignmentType::TemplateSwitchEntrance {
        first_offset,
        uncertainty_range: Some(TSMUncertaintyRange::new()),
        descendant,
        ancestor,
        direction,
    }
}

fn exit(anti_descendant_gap: isize) -> AlignmentType {
    AlignmentType::TemplateSwitchExit {
        anti_descendant_gap,
    }
}

/// A forward TS with a substitution in its inner, whose ancestor character is hidden,
/// because it is the inner of the TS itself (descendant and ancestor are the reference).
#[test]
fn forward_substitution_against_hidden_ancestor() {
    let result = build_result(
        "GAAA",
        "TGAA",
        vec![
            (
                1,
                entrance(
                    0,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Reference,
                    TemplateSwitchDirection::Forward,
                ),
            ),
            (2, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (1, exit(1)),
            (
                1,
                entrance(
                    -1,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Reference,
                    TemplateSwitchDirection::Forward,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (1, AlignmentType::SecondaryDeletion),
            (1, AlignmentType::SecondaryMatch),
            (1, exit(-1)),
            (2, AlignmentType::PrimaryInsertion),
        ],
        AlignmentRange::new_complete(4, 4),
    );
    render_svg(&result, &None, &default_svg_config()).unwrap();
}

/// `--restrict-context-to` with a TS whose hidden inner was removed as an empty column before restricting the context.
#[test]
fn restricted_context_with_removed_hidden_inner() {
    let reverse = TemplateSwitchDirection::Reverse;
    let result = build_result(
        "GCTTACTATA",
        "TTATTGATGAGAAGC",
        vec![
            (7, AlignmentType::PrimaryMatch),
            (
                1,
                entrance(
                    -1,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Reference,
                    reverse,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (1, exit(3)),
            (
                1,
                entrance(
                    -5,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Query,
                    reverse,
                ),
            ),
            (4, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (1, exit(2)),
        ],
        AlignmentRange::new_complete(10, 15),
    );
    let mut config = default_svg_config();
    config.restrict_context = Some(2);
    render_svg(&result, &None, &config).unwrap();
}

fn entrance_with_uncertainty(
    first_offset: isize,
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
    [
        start_left_shift,
        start_right_shift,
        end_left_shift,
        end_right_shift,
    ]: [u16; 4],
) -> AlignmentType {
    AlignmentType::TemplateSwitchEntrance {
        first_offset,
        uncertainty_range: Some(TSMUncertaintyRange {
            start_left_shift,
            start_right_shift,
            end_left_shift,
            end_right_shift,
        }),
        descendant,
        ancestor,
        direction: TemplateSwitchDirection::Reverse,
    }
}

fn render_with_all_uncertainty_modes(result: &common::Result) {
    for mode in [UncertaintyRangeMode::InnerOnly, UncertaintyRangeMode::Full] {
        let mut config = default_svg_config();
        config.uncertainty_range_mode = mode;
        render_svg(result, &None, &config).unwrap();
    }
}

/// The uncertainty range of a TS inner that is at the right end of the arrangement
/// extends beyond the arrangement (index out of bounds in `inner.rs`).
#[test]
fn uncertainty_range_beyond_the_right_end() {
    let result = build_result(
        "CGTTCAACCTGACCGAT",
        "CCCGTCGCTT",
        vec![
            (2, AlignmentType::PrimaryMatch),
            (2, AlignmentType::PrimarySubstitution),
            (4, AlignmentType::PrimaryMatch),
            (
                1,
                entrance_with_uncertainty(
                    8,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Reference,
                    [2, 0, 1, 0],
                ),
            ),
            (2, AlignmentType::SecondaryMatch),
            (1, exit(-2)),
            (11, AlignmentType::PrimaryDeletion),
        ],
        AlignmentRange::new_complete(17, 10),
    );
    render_with_all_uncertainty_modes(&result);
}

/// The uncertainty range of a TS inner that is at the left end of the arrangement
/// extends beyond the arrangement (subtraction overflow of an arrangement column).
#[test]
fn uncertainty_range_beyond_the_left_end() {
    let result = build_result(
        "GAGACCGGGACAGG",
        "CGCCGCCCCTAT",
        vec![
            (
                1,
                entrance(
                    8,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Query,
                    TemplateSwitchDirection::Reverse,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (2, AlignmentType::SecondarySubstitution),
            (1, exit(12)),
            (
                1,
                entrance_with_uncertainty(
                    -9,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Reference,
                    [0, 0, 1, 2],
                ),
            ),
            (2, AlignmentType::SecondaryMatch),
            (1, exit(-2)),
            (4, AlignmentType::PrimaryMatch),
            (3, AlignmentType::PrimaryInsertion),
        ],
        AlignmentRange::new_complete(14, 12),
    );
    render_with_all_uncertainty_modes(&result);
}

/// Issue #285: overlapping forward TSs (found by the property tests).
#[test]
fn overlapping_forward_template_switches() {
    let forward = TemplateSwitchDirection::Forward;
    let result = build_result(
        "TACGAACATGAACGCT",
        "CAAGTCAGGGAACTCCCAGATATGGGTC",
        vec![
            (4, AlignmentType::PrimaryInsertion),
            (
                1,
                entrance(
                    1,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Reference,
                    forward,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondaryDeletion),
            (2, AlignmentType::SecondarySubstitution),
            (1, AlignmentType::SecondaryInsertion),
            (1, AlignmentType::SecondaryMatch),
            (1, exit(20)),
            (
                1,
                entrance(
                    -1,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Reference,
                    TemplateSwitchDirection::Reverse,
                ),
            ),
            (1, AlignmentType::SecondaryMatch),
            (1, exit(0)),
            (
                1,
                entrance(
                    -10,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Query,
                    forward,
                ),
            ),
            (3, AlignmentType::SecondaryMatch),
            (1, exit(2)),
            (9, AlignmentType::PrimaryDeletion),
        ],
        AlignmentRange::new_complete(16, 28),
    );
    render_svg(&result, &None, &default_svg_config()).unwrap();
}

/// Known bug: overlapping reverse TSs with uncertainty ranges (found by the property tests).
#[test]
#[ignore = "issue #286"]
fn overlapping_reverse_template_switches_with_uncertainty() {
    use TemplateSwitchAncestor as A;
    use TemplateSwitchDescendant as D;
    let reverse = TemplateSwitchDirection::Reverse;
    let mut alignment = vec![
        (4, AlignmentType::PrimaryDeletion),
        (
            1,
            entrance_with_uncertainty(18, D::Query, A::Query, [0, 0, 1, 0]),
        ),
        (1, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondaryInsertion),
        (2, AlignmentType::SecondaryMatch),
        (1, exit(1)),
    ];
    let without_uncertainty = |first_offset, descendant, ancestor| {
        let mut entrance = entrance(first_offset, descendant, ancestor, reverse);
        if let AlignmentType::TemplateSwitchEntrance {
            uncertainty_range, ..
        } = &mut entrance
        {
            *uncertainty_range = None;
        }
        entrance
    };
    alignment.extend([
        (1, without_uncertainty(9, D::Query, A::Reference)),
        (1, AlignmentType::SecondaryMatch),
        (1, exit(18)),
        (1, without_uncertainty(-12, D::Query, A::Reference)),
        (1, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondaryDeletion),
        (1, exit(-6)),
        (1, without_uncertainty(-7, D::Query, A::Reference)),
        (1, AlignmentType::SecondaryMatch),
        (1, exit(-2)),
        (1, without_uncertainty(7, D::Reference, A::Reference)),
        (1, AlignmentType::SecondaryMatch),
        (1, AlignmentType::SecondaryDeletion),
        (2, AlignmentType::SecondaryInsertion),
        (2, AlignmentType::SecondaryMatch),
        (1, exit(-6)),
    ]);
    let result = build_result(
        "CAGTCAAGTCGTGGATAGAGACAGAAT",
        "AACCGCCCTTGAGCTGCCAGAAAGTTC",
        alignment,
        AlignmentRange::new(1, 3, 21, 4),
    );
    render_with_all_uncertainty_modes(&result);
}
