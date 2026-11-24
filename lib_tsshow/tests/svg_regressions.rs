//! Regression tests of the SVG output for bugs found by the property tests (`fuzz.rs`).

use common::{build_result, default_svg_config, render_svg};
use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_geometry::AlignmentRange,
    alignment_result::{AlignmentType, TSMUncertaintyRange},
};

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
