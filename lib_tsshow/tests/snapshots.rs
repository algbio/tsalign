//! Golden-file (snapshot) tests of the rendered plain text and SVG output.
//!
//! The snapshots live in `tests/snapshots`. To update them after an intended change of the output, run
//! `TSSHOW_BLESS=1 cargo test -p lib_tsshow --test snapshots` and review the diff of the snapshots.
//! The embedded fonts are replaced by a placeholder in the SVG snapshots (see [`common::strip_embedded_fonts`]).

use common::{
    Fixture, assert_snapshot, default_svg_config, fixture_names, render_plain_text, render_svg,
    strip_embedded_fonts, svg_config, test_files_fixtures,
};
use lib_tsshow::svg::UncertaintyRangeMode;

mod common;

/// Snapshots of every fixture in `tests/fixtures`: plain text, and SVG with the default and with all options.
#[test]
fn fixture_snapshots() {
    let names = fixture_names();
    assert!(names.len() >= 20, "fixtures are missing: {names:?}");
    let mut failures = Vec::new();

    for name in names {
        let fixture = Fixture::load(&name);
        let result = std::panic::catch_unwind(|| fixture_snapshot(&fixture));
        if let Err(error) = result {
            failures.push(format!("{name}: {}", panic_message(&error)));
        }
    }

    assert!(
        failures.is_empty(),
        "{} fixtures do not match their snapshots:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn fixture_snapshot(fixture: &Fixture) {
    let Fixture {
        name,
        result,
        no_ts_result,
    } = fixture;

    assert_snapshot(
        &format!("{name}.txt"),
        &render_plain_text(result, no_ts_result),
    );
    assert_snapshot(
        &format!("{name}.svg"),
        &strip_embedded_fonts(&render_svg(result, no_ts_result, &default_svg_config()).unwrap()),
    );
    assert_snapshot(
        &format!("{name}.all_options.svg"),
        &strip_embedded_fonts(
            &render_svg(
                result,
                no_ts_result,
                &svg_config(true, true, Some(5), UncertaintyRangeMode::Full),
            )
            .unwrap(),
        ),
    );
}

/// Snapshots of the alignments in `test_files` with the options of the old show-bug scripts
/// (`test_files/ari_200_show_bug.sh`, `test_files/ari_868_show_bug.sh`: `-are -z 5`).
#[test]
fn test_files_snapshots() {
    let mut failures = Vec::new();

    for fixture in test_files_fixtures() {
        let result = std::panic::catch_unwind(|| {
            let Fixture {
                name,
                result,
                no_ts_result,
            } = &fixture;
            assert_snapshot(
                &format!("test_files/{name}.txt"),
                &render_plain_text(result, no_ts_result),
            );
            assert_snapshot(
                &format!("test_files/{name}.svg"),
                &strip_embedded_fonts(
                    &render_svg(
                        result,
                        no_ts_result,
                        &svg_config(true, false, Some(5), UncertaintyRangeMode::Full),
                    )
                    .unwrap(),
                ),
            );
        });
        if let Err(error) = result {
            failures.push(format!("{}: {}", fixture.name, panic_message(&error)));
        }
    }

    assert!(
        failures.is_empty(),
        "{} test_files alignments do not match their snapshots:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The SVG options that are not covered by the fixture snapshots, on one fixture each.
#[test]
fn svg_option_snapshots() {
    let fixture = Fixture::load("classic_both_descendants");
    let configs = [
        (
            "arrows",
            svg_config(true, false, None, UncertaintyRangeMode::None),
        ),
        (
            "more_complement",
            svg_config(false, true, None, UncertaintyRangeMode::None),
        ),
        (
            "context_0",
            svg_config(false, false, Some(0), UncertaintyRangeMode::None),
        ),
        (
            "context_3",
            svg_config(false, false, Some(3), UncertaintyRangeMode::None),
        ),
        (
            "uncertainty_inner_only",
            svg_config(false, false, None, UncertaintyRangeMode::InnerOnly),
        ),
        (
            "uncertainty_full",
            svg_config(false, false, None, UncertaintyRangeMode::Full),
        ),
    ];

    for (config_name, config) in configs {
        assert_snapshot(
            &format!("options/{}.{config_name}.svg", fixture.name),
            &strip_embedded_fonts(
                &render_svg(&fixture.result, &fixture.no_ts_result, &config).unwrap(),
            ),
        );
    }

    // Without the no-TS alignment.
    assert_snapshot(
        &format!("options/{}.without_no_ts.svg", fixture.name),
        &strip_embedded_fonts(&render_svg(&fixture.result, &None, &default_svg_config()).unwrap()),
    );
    assert_snapshot(
        &format!("options/{}.without_no_ts.txt", fixture.name),
        &render_plain_text(&fixture.result, &None),
    );
}

/// Rendering the same input twice gives the same output, also in the embedded fonts.
#[test]
fn rendering_is_deterministic() {
    let fixture = Fixture::load("classic_three_ts");
    let config = svg_config(true, true, Some(5), UncertaintyRangeMode::Full);

    let svg = render_svg(&fixture.result, &fixture.no_ts_result, &config).unwrap();
    assert_eq!(
        svg,
        render_svg(&fixture.result, &fixture.no_ts_result, &config).unwrap()
    );

    let text = render_plain_text(&fixture.result, &fixture.no_ts_result);
    assert_eq!(
        text,
        render_plain_text(&fixture.result, &fixture.no_ts_result)
    );
}

fn panic_message(error: &Box<dyn std::any::Any + Send>) -> String {
    error
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| error.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_else(|| "<unknown panic>".to_owned())
}
