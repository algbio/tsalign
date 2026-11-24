//! Structural checks of the SVG output that do not depend on the exact layout.
//!
//! For every fixture and several configurations, the SVG must be well-formed XML with a valid view box,
//! contain the expected labels, numbers and arrows for each template switch, and show the characters of the sequences:
//! - the rows of the no-TS alignment contain every character of their sequence exactly once, in order;
//! - the reference and query rows of the TS alignment contain every character of their sequence exactly once,
//!   except for the TS inners of which they are the descendant, which are shown in the inner rows instead;
//! - the inner rows contain the characters of the TS inners (reversed for reverse TSs).

use std::collections::BTreeMap;

use common::{Fixture, Result, fixture_names, render_svg, svg_config, test_files_fixtures};
use lib_tsalign::{
    aligner::{TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_result::{AlignmentResult, AlignmentType},
};
use lib_tsshow::svg::{SvgConfig, UncertaintyRangeMode};
use roxmltree::{Document, Node};

mod common;

/// The fill colours of copied characters (see `COPY_COLORS` in `svg.rs`).
const COPY_COLORS: &[&str] = &["#00CC00", "#009900", "#006600", "#003300"];
/// The width of a character of the sequences in the SVG.
const CHARACTER_WIDTH: f32 = 2.0;
const LABEL_COLOR: &str = "#555555";

#[test]
fn fixtures_are_structurally_valid() {
    let mut failures = Vec::new();
    for name in fixture_names() {
        let fixture = Fixture::load(&name);
        check_fixture(&fixture, &mut failures);
    }
    assert!(
        failures.is_empty(),
        "{} structural checks failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn test_files_are_structurally_valid() {
    let mut failures = Vec::new();
    for fixture in test_files_fixtures() {
        check_fixture(&fixture, &mut failures);
    }
    assert!(
        failures.is_empty(),
        "{} structural checks failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn configs() -> Vec<(&'static str, SvgConfig)> {
    vec![
        (
            "default",
            svg_config(false, false, None, UncertaintyRangeMode::None),
        ),
        (
            "arrows",
            svg_config(true, false, None, UncertaintyRangeMode::None),
        ),
        (
            "more_complement",
            svg_config(true, true, None, UncertaintyRangeMode::None),
        ),
        (
            "inner_only",
            svg_config(false, false, None, UncertaintyRangeMode::InnerOnly),
        ),
        (
            "full",
            svg_config(true, false, None, UncertaintyRangeMode::Full),
        ),
        (
            "context_5",
            svg_config(true, false, Some(5), UncertaintyRangeMode::None),
        ),
        (
            "context_0_full",
            svg_config(true, true, Some(0), UncertaintyRangeMode::Full),
        ),
    ]
}

fn check_fixture(fixture: &Fixture, failures: &mut Vec<String>) {
    for (config_name, config) in configs() {
        for with_no_ts in [true, false] {
            let no_ts_result = if with_no_ts {
                fixture.no_ts_result.clone()
            } else {
                None
            };
            let context = format!(
                "{} [{config_name}{}]",
                fixture.name,
                if with_no_ts { ", with no-TS" } else { "" }
            );
            let outcome = std::panic::catch_unwind(|| {
                let svg = render_svg(&fixture.result, &no_ts_result, &config)
                    .unwrap_or_else(|error| panic!("rendering failed: {error}"));
                check_svg(&svg, &fixture.result, &no_ts_result, &config)
            });
            match outcome {
                Ok(Ok(())) => {}
                Ok(Err(error)) => failures.push(format!("{context}: {error}")),
                Err(error) => failures.push(format!(
                    "{context}: panicked: {}",
                    error
                        .downcast_ref::<String>()
                        .cloned()
                        .or_else(|| error.downcast_ref::<&str>().map(|s| s.to_string()))
                        .unwrap_or_default()
                )),
            }
        }
    }
}

/// A template switch as seen from the alignment, computed independently of lib_tsshow.
#[derive(Debug)]
struct ExpectedTemplateSwitch {
    descendant: TemplateSwitchDescendant,
    direction: TemplateSwitchDirection,
    /// The descendant range of the inner.
    sp1: usize,
    sp4: usize,
}

fn expected_template_switches(result: &Result) -> Vec<ExpectedTemplateSwitch> {
    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        panic!("fixture has no alignment");
    };
    let mut reference = statistics.reference_offset;
    let mut query = statistics.query_offset;
    let mut current = None;
    let mut result = Vec::new();

    for alignment_type in alignment.iter_flat_cloned() {
        match alignment_type {
            AlignmentType::PrimaryMatch | AlignmentType::PrimarySubstitution => {
                reference += 1;
                query += 1;
            }
            AlignmentType::PrimaryInsertion => query += 1,
            AlignmentType::PrimaryDeletion => reference += 1,
            AlignmentType::TemplateSwitchEntrance {
                descendant,
                direction,
                ..
            } => {
                let sp1 = match descendant {
                    TemplateSwitchDescendant::Reference => reference,
                    TemplateSwitchDescendant::Query => query,
                };
                current = Some(ExpectedTemplateSwitch {
                    descendant,
                    direction,
                    sp1,
                    sp4: sp1,
                });
            }
            AlignmentType::SecondaryMatch
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryInsertion => current.as_mut().unwrap().sp4 += 1,
            AlignmentType::SecondaryDeletion => {}
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap,
            } => {
                let ts = current.take().unwrap();
                match ts.descendant {
                    TemplateSwitchDescendant::Reference => {
                        reference = ts.sp4;
                        query = query.checked_add_signed(anti_descendant_gap).unwrap();
                    }
                    TemplateSwitchDescendant::Query => {
                        query = ts.sp4;
                        reference = reference.checked_add_signed(anti_descendant_gap).unwrap();
                    }
                }
                result.push(ts);
            }
            other => panic!("unexpected alignment type {other}"),
        }
    }

    result
}

/// A row of characters in the SVG, by column.
#[derive(Debug, Default)]
struct Row {
    characters: BTreeMap<usize, (char, String)>,
}

impl Row {
    /// The characters of the row from left to right, without gaps, separators and copies, in upper case.
    fn source_characters(&self) -> String {
        self.characters
            .values()
            .filter(|(c, fill)| *c != '-' && *c != '|' && !COPY_COLORS.contains(&fill.as_str()))
            .map(|(c, _)| c.to_ascii_uppercase())
            .collect()
    }
}

/// The absolute translation of a node, summing up all `translate(x y)` transforms of its ancestors.
fn translation(node: Node) -> (f32, f32) {
    let mut x = 0.0;
    let mut y = 0.0;
    for ancestor in node.ancestors() {
        if let Some(transform) = ancestor.attribute("transform") {
            if let Some(arguments) = transform
                .strip_prefix("translate(")
                .and_then(|t| t.strip_suffix(')'))
            {
                let numbers: Vec<f32> = arguments
                    .split([' ', ','])
                    .map(|n| n.parse().unwrap())
                    .collect();
                x += numbers[0];
                y += numbers[1];
            }
        }
    }
    (x, y)
}

fn text(node: Node) -> String {
    node.text().unwrap_or_default().trim().to_owned()
}

fn is_subsequence(needle: &str, haystack: &str) -> bool {
    let mut haystack = haystack.chars();
    needle.chars().all(|c| haystack.any(|h| h == c))
}

fn check_svg(
    svg: &str,
    result: &Result,
    no_ts_result: &Option<Result>,
    config: &SvgConfig,
) -> std::result::Result<(), String> {
    let document = Document::parse(svg).map_err(|error| format!("invalid XML: {error}"))?;
    let root = document.root_element();
    if root.tag_name().name() != "svg" {
        return Err("root element is not <svg>".to_owned());
    }

    // View box.
    let view_box: Vec<f32> = root
        .attribute("viewBox")
        .ok_or("no viewBox")?
        .split_whitespace()
        .map(|n| n.parse().map_err(|_| format!("invalid viewBox number {n}")))
        .collect::<std::result::Result<_, _>>()?;
    if view_box.len() != 4
        || view_box.iter().any(|n| !n.is_finite())
        || view_box[2] <= 0.0
        || view_box[3] <= 0.0
    {
        return Err(format!("invalid viewBox {view_box:?}"));
    }

    // No coordinate may be NaN or infinite.
    for node in document.descendants() {
        for attribute in node.attributes() {
            if attribute.value().contains("NaN") || attribute.value().contains("inf") {
                return Err(format!(
                    "attribute {}={:?} is not finite",
                    attribute.name(),
                    attribute.value()
                ));
            }
        }
    }

    // Fonts.
    let style = document
        .descendants()
        .find(|n| n.has_tag_name("style"))
        .ok_or("no <style>")?;
    let style = style.text().unwrap_or_default();
    for family in ["Liberation Mono", "DejaVu Sans Mono", "DejaVu Sans"] {
        if style
            .matches(&format!("font-family: \"{family}\";"))
            .count()
            != 1
        {
            return Err(format!("font {family} is not embedded exactly once"));
        }
    }

    let statistics = result.statistics();
    let sequences = &statistics.sequences;
    let template_switches = expected_template_switches(result);

    // Labels and rows.
    let labels: Vec<(String, f32)> = document
        .descendants()
        .filter(|n| n.has_tag_name("text") && n.attribute("fill") == Some(LABEL_COLOR))
        .map(|n| {
            let (_, y) = translation(n);
            let y = y + n.attribute("y").unwrap().parse::<f32>().unwrap();
            (text(n), y)
        })
        .collect();
    let mut rows: Vec<(f32, Row)> = Vec::new();
    for group in document
        .descendants()
        .filter(|n| n.has_tag_name("g") && n.attribute("font-family") == Some("Liberation Mono"))
    {
        let (_, y) = translation(group);
        let mut row = Row::default();
        for character in group.children().filter(|n| n.has_tag_name("text")) {
            let x: f32 = character.attribute("x").unwrap().parse().unwrap();
            let column = (x / CHARACTER_WIDTH).round() as usize;
            let c = text(character);
            if c.chars().count() != 1 {
                return Err(format!("character element contains {c:?}"));
            }
            let fill = character.attribute("fill").unwrap_or_default().to_owned();
            if row
                .characters
                .insert(column, (c.chars().next().unwrap(), fill))
                .is_some()
            {
                return Err(format!("two characters in column {column} of row at y={y}"));
            }
        }
        rows.push((y, row));
    }

    let row_of = |label: &str, occurrence: usize| -> std::result::Result<&Row, String> {
        let (_, y) = labels
            .iter()
            .filter(|(l, _)| l == label)
            .nth(occurrence)
            .ok_or_else(|| format!("label {label:?} #{occurrence} is missing"))?;
        rows.iter()
            .find(|(row_y, _)| (row_y - y).abs() < 1e-3)
            .map(|(_, row)| row)
            .ok_or_else(|| format!("no row for label {label:?} #{occurrence}"))
    };

    let count_labels = |label: &str| labels.iter().filter(|(l, _)| l == label).count();
    let no_ts_rows = usize::from(no_ts_result.is_some());
    for (label, expected) in [
        ("Reference complement:", 1),
        ("Query complement:", 1),
        ("Reference:", 1 + no_ts_rows),
        ("Query:", 1 + no_ts_rows),
    ] {
        if count_labels(label) != expected {
            return Err(format!(
                "label {label:?} occurs {} times instead of {expected}",
                count_labels(label)
            ));
        }
    }

    let running_numbers: Vec<char> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        .chars()
        .take(template_switches.len())
        .collect();
    let inner_labels = labels
        .iter()
        .filter(|(l, _)| l.ends_with(" inner:"))
        .count();
    if inner_labels != template_switches.len() {
        return Err(format!(
            "{inner_labels} inner rows for {} template switches",
            template_switches.len()
        ));
    }

    // Numbers at the switch points.
    let texts: Vec<String> = document
        .descendants()
        .filter(|n| n.has_tag_name("text"))
        .map(text)
        .collect();
    for running_number in &running_numbers {
        if count_labels(&format!("TS-{running_number} inner:")) != 1 {
            return Err(format!("inner label of TS {running_number} is missing"));
        }
        for switch_point in 1..=4 {
            let number = format!("{running_number}{switch_point}");
            let count = texts.iter().filter(|t| **t == number).count();
            if count != 1 {
                return Err(format!("number {number} occurs {count} times"));
            }
        }
    }

    // Arrows.
    let arrows = document
        .descendants()
        .filter(|n| n.attribute("marker-end") == Some("url(#arrow_head_red)"))
        .count();
    let expected_arrows = if config.render_arrows {
        2 * template_switches.len()
    } else {
        0
    };
    if arrows != expected_arrows {
        return Err(format!("{arrows} arrows instead of {expected_arrows}"));
    }

    // Legend with the sequence names.
    for name in [&sequences.reference_name, &sequences.query_name] {
        if !texts.iter().any(|t| t == name.trim()) {
            return Err(format!("the legend does not contain the name {name:?}"));
        }
    }

    // Characters of the no-TS alignment.
    if no_ts_result.is_some() {
        for (label, sequence) in [
            ("Reference:", &sequences.reference),
            ("Query:", &sequences.query),
        ] {
            let row = row_of(label, 1)?;
            let characters = row.source_characters();
            let complete = config.restrict_context.is_none();
            if (complete && characters != *sequence)
                || (!complete && !sequence.contains(&characters))
            {
                return Err(format!(
                    "no-TS row {label} shows {characters:?} instead of {} {sequence:?}",
                    if complete { "" } else { "a substring of" }
                ));
            }
        }
    }

    // Characters of the TS alignment.
    for (label, sequence, descendant) in [
        (
            "Reference:",
            &sequences.reference,
            TemplateSwitchDescendant::Reference,
        ),
        ("Query:", &sequences.query, TemplateSwitchDescendant::Query),
    ] {
        let row = row_of(label, 0)?;
        let characters = row.source_characters();
        // The sequence without the inners it is the descendant of.
        let mut expected = String::new();
        let mut position = 0;
        let mut inners: Vec<_> = template_switches
            .iter()
            .filter(|ts| ts.descendant == descendant)
            .map(|ts| (ts.sp1, ts.sp4))
            .collect();
        inners.sort();
        for (sp1, sp4) in inners {
            if sp1 >= position {
                expected.push_str(&sequence[position..sp1]);
            }
            position = position.max(sp4);
        }
        expected.push_str(&sequence[position..]);

        let exact = config.restrict_context.is_none()
            && config.uncertainty_range_mode != UncertaintyRangeMode::Full;
        if exact && characters != expected {
            return Err(format!(
                "TS row {label} shows {characters:?} instead of {expected:?}"
            ));
        }
        if !exact && !is_subsequence(&characters, sequence) {
            return Err(format!(
                "TS row {label} shows {characters:?}, which is not a subsequence of {sequence:?}"
            ));
        }
    }

    // Characters of the inners.
    if config.restrict_context.is_none()
        && config.uncertainty_range_mode == UncertaintyRangeMode::None
    {
        for (running_number, ts) in running_numbers.iter().zip(&template_switches) {
            let row = row_of(&format!("TS-{running_number} inner:"), 0)?;
            let characters: String = row
                .characters
                .values()
                .filter(|(c, _)| *c != '-')
                .map(|(c, _)| c.to_ascii_uppercase())
                .collect();
            let descendant_sequence = match ts.descendant {
                TemplateSwitchDescendant::Reference => &sequences.reference,
                TemplateSwitchDescendant::Query => &sequences.query,
            };
            let inner = &descendant_sequence[ts.sp1..ts.sp4];
            let expected: String = match ts.direction {
                TemplateSwitchDirection::Forward => inner.to_owned(),
                TemplateSwitchDirection::Reverse => inner.chars().rev().collect(),
            };
            if characters != expected {
                return Err(format!(
                    "inner row of TS {running_number} shows {characters:?} instead of {expected:?}"
                ));
            }
        }
    }

    Ok(())
}
