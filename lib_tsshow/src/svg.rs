use std::{collections::HashMap, io::Write};

use arrows::{Arrow, ArrowEndpointDirection, add_arrow_defs};
use font::{CharacterData, sans_serif_mono, svg_string, typewriter};
use indexed_str::IndexedStr;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::AlignmentResult,
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info};
use numbers::{Number, NumberAlignment};
use svg::{
    Document,
    node::element::{Circle, Group, Text},
};

use crate::{
    error::{Error, Result},
    plain_text::mutlipair_alignment_renderer::Character,
    svg::font::svg_text,
    ts_arrangement::{
        TsArrangement,
        character::Char,
        complement::ComplementChar,
        index_types::ArrangementColumn,
        inner::InnerChar,
        row::TsArrangementRow,
        source::{SourceChar, TsSourceArrangement},
        template_switch::TemplateSwitch,
    },
};

mod arrows;
mod font;
mod indexed_str;
pub mod labelled_sequence;
mod numbers;

const SVG_PADDING: f32 = 10.0;
const COPY_COLORS: &[&str] = &["#00CC00", "#009900", "#006600", "#003300"];
const COMPLEMENT_SOURCE_HIDDEN_COLOR: &str = "grey";
const TS_RUNNING_NUMBER: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

struct SvgLocation {
    pub x: f32,
    pub y: f32,
}

impl SvgLocation {
    pub fn as_transform(&self) -> String {
        format!("translate({} {})", self.x, self.y)
    }
}

pub struct SvgConfig {
    pub render_arrows: bool,
    pub render_more_complement: bool,
    /// Restrict the context around the template switches to this many characters on each side.
    /// If `None`, the full sequences will be rendered.
    pub restrict_context: Option<usize>,
}

pub fn create_ts_svg(
    output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    config: &SvgConfig,
) -> Result<()> {
    info!("Creating template switch SVG");

    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        return Err(Error::AlignmentHasNoTarget);
    };
    debug!("Alignment: {}", alignment.cigar());

    let reference = &statistics.sequences.reference;
    let query = &statistics.sequences.query;
    let reference_c: String = statistics.sequences.reference_rc.chars().rev().collect();
    let query_c: String = statistics.sequences.query_rc.chars().rev().collect();

    debug!("Computing TS arrangement");
    let mut ts_arrangement = TsArrangement::new(
        statistics.reference_offset,
        statistics.query_offset,
        reference.len(),
        query.len(),
        alignment.iter_flat_cloned(),
    )?;

    if config.render_more_complement {
        ts_arrangement.show_complete_complements_if_used();
    }

    ts_arrangement.remove_empty_columns();

    let (reference_context_range, query_context_range) =
        if let Some(context_character_amount) = config.restrict_context {
            ts_arrangement.limit_context_to(context_character_amount)
        } else {
            (
                0usize.into()..reference.len().into(),
                0usize.into()..query.len().into(),
            )
        };

    let ts_arrangement = ts_arrangement;

    let no_ts_arrangement = no_ts_result
        .as_ref()
        .map(|no_ts_result| {
            debug!("Computing no-TS arrangement");
            let AlignmentResult::WithTarget {
                alignment,
                statistics,
            } = no_ts_result
            else {
                return Err(Error::NoTsAlignmentHasNoTarget);
            };

            let reference = &statistics.sequences.reference;
            let query = &statistics.sequences.query;
            Ok(TsSourceArrangement::new(
                statistics.reference_offset,
                statistics.query_offset,
                reference.len(),
                query.len(),
                alignment.iter_flat_cloned(),
                &mut Vec::new(),
            )
            .map(|mut no_ts_arrangement| {
                let min_start = no_ts_arrangement
                    .reference_source_to_arrangement_column(reference_context_range.start)
                    .min(
                        no_ts_arrangement
                            .query_source_to_arrangement_column(query_context_range.start),
                    );
                let limit = no_ts_arrangement.width();
                let max_end = no_ts_arrangement
                    .try_reference_source_to_arrangement_column(reference_context_range.end)
                    .unwrap_or(limit.into())
                    .max(
                        no_ts_arrangement
                            .try_query_source_to_arrangement_column(query_context_range.end)
                            .unwrap_or(limit.into()),
                    );

                no_ts_arrangement.remove_columns(
                    (0usize..min_start.primitive())
                        .chain(max_end.primitive()..limit)
                        .map(ArrangementColumn::new),
                );
                no_ts_arrangement
            }))
        })
        .transpose()?
        .transpose()?;

    debug!("Computing TS arrows and labels");
    let mut arrows = Vec::new();
    let mut numbers = Vec::new();

    for (
        inner_identifier,
        TemplateSwitch {
            index: ts_index,
            primary,
            secondary,
            sp1_reference,
            sp1_query,
            sp2_secondary,
            sp3_secondary,
            sp4_reference,
            sp4_query,
            ..
        },
    ) in ts_arrangement.template_switches()
    {
        /*let anti_primary_sp4_minus_one = match primary {
            TemplateSwitchPrimary::Reference => sp4_query.checked_sub(1).map(|column| {
                ts_arrangement.query_arrangement_char_to_arrangement_column(column) + 1usize
            }),
            TemplateSwitchPrimary::Query => sp4_reference.checked_sub(1).map(|column| {
                ts_arrangement.reference_arrangement_char_to_arrangement_column(column) + 1usize
            }),
        };

        let (anti_primary_sp4, anti_primary_row) = match primary {
            TemplateSwitchPrimary::Reference => (
                ts_arrangement
                    .query()
                    .iter()
                    .take(
                        usize::from(
                            ts_arrangement.query_arrangement_char_to_arrangement_column(*sp4_query),
                        ) + 1usize,
                    )
                    .skip(anti_primary_sp4_minus_one.unwrap_or(0.into()).into())
                    .find(|(_, c)| !c.is_blank())
                    .unwrap()
                    .0,
                TsArrangementRow::Query,
            ),
            TemplateSwitchPrimary::Query => (
                ts_arrangement
                    .reference()
                    .iter()
                    .take(
                        usize::from(
                            ts_arrangement
                                .reference_arrangement_char_to_arrangement_column(*sp4_reference),
                        ) + 1usize,
                    )
                    .skip(anti_primary_sp4_minus_one.unwrap_or(0.into()).into())
                    .find(|(_, c)| !c.is_blank())
                    .unwrap()
                    .0,
                TsArrangementRow::Reference,
            ),
        };

        if let Some(anti_primary_sp4_minus_one) = anti_primary_sp4_minus_one {
            if (anti_primary_sp4 - anti_primary_sp4_minus_one) > ArrangementColumn::ZERO {
                arrows.push(Arrow::new_skip(
                    anti_primary_sp4_minus_one,
                    anti_primary_sp4,
                    anti_primary_row,
                ));
            }
        }*/

        let primary_sp4_minus_one = match primary {
            TemplateSwitchPrimary::Reference => sp4_reference.checked_sub(1).map(|column| {
                ts_arrangement.reference_arrangement_char_to_arrangement_column(column) + 1usize
            }),
            TemplateSwitchPrimary::Query => sp4_query.checked_sub(1).map(|column| {
                ts_arrangement.query_arrangement_char_to_arrangement_column(column) + 1usize
            }),
        };

        let (primary_sp1, primary_sp4, primary_row) = match primary {
            TemplateSwitchPrimary::Reference => (
                ts_arrangement.reference_arrangement_char_to_arrangement_column(*sp1_reference),
                ts_arrangement
                    .reference()
                    .iter()
                    .take(
                        usize::from(
                            ts_arrangement
                                .reference_arrangement_char_to_arrangement_column(*sp4_reference),
                        ) + 1usize,
                    )
                    .skip(primary_sp4_minus_one.unwrap_or(0.into()).into())
                    .find(|(_, c)| !c.is_blank())
                    .unwrap()
                    .0,
                TsArrangementRow::Reference,
            ),
            TemplateSwitchPrimary::Query => (
                ts_arrangement.query_arrangement_char_to_arrangement_column(*sp1_query),
                ts_arrangement
                    .query()
                    .iter()
                    .take(
                        usize::from(
                            ts_arrangement.query_arrangement_char_to_arrangement_column(*sp4_query),
                        ) + 1usize,
                    )
                    .skip(primary_sp4_minus_one.unwrap_or(0.into()).into())
                    .find(|(_, c)| !c.is_blank())
                    .unwrap()
                    .0,
                TsArrangementRow::Query,
            ),
        };

        let forward = sp2_secondary < sp3_secondary;
        let (secondary_limit, secondary_offset, inner_row) = match secondary {
            TemplateSwitchSecondary::Reference => (
                ts_arrangement.inner_last_non_blank_column(inner_identifier) + 1usize,
                ts_arrangement.inner_first_non_blank_column(inner_identifier),
                TsArrangementRow::Inner {
                    index: inner_identifier,
                },
            ),
            TemplateSwitchSecondary::Query => (
                ts_arrangement.inner_last_non_blank_column(inner_identifier) + 1usize,
                ts_arrangement.inner_first_non_blank_column(inner_identifier),
                TsArrangementRow::Inner {
                    index: inner_identifier,
                },
            ),
        };

        let running_number = TS_RUNNING_NUMBER.chars().nth(*ts_index).unwrap();
        let number1 = Number::new(
            format!("{running_number}1"),
            primary_sp1,
            primary_row,
            NumberAlignment::Left,
            0.5,
        );
        let number2 = Number::new(
            format!("{running_number}2"),
            if forward {
                secondary_offset
            } else {
                secondary_limit
            },
            inner_row,
            if forward {
                NumberAlignment::Right
            } else {
                NumberAlignment::Left
            },
            0.5,
        );
        let number3 = Number::new(
            format!("{running_number}3"),
            if forward {
                secondary_limit
            } else {
                secondary_offset
            },
            inner_row,
            if forward {
                NumberAlignment::Left
            } else {
                NumberAlignment::Right
            },
            0.5,
        );
        let number4 = Number::new(
            format!("{running_number}4"),
            primary_sp4,
            primary_row,
            NumberAlignment::Right,
            0.5,
        );

        if config.render_arrows {
            arrows.push(Arrow::new_curved(
                primary_sp1,
                number1.width(),
                primary_row,
                ArrowEndpointDirection::Forward,
                if forward {
                    secondary_offset
                } else {
                    secondary_limit
                },
                number2.width(),
                inner_row,
                if forward {
                    ArrowEndpointDirection::Backward
                } else {
                    ArrowEndpointDirection::Forward
                },
            ));

            arrows.push(Arrow::new_curved(
                if forward {
                    secondary_limit
                } else {
                    secondary_offset
                },
                number3.width(),
                inner_row,
                if forward {
                    ArrowEndpointDirection::Forward
                } else {
                    ArrowEndpointDirection::Backward
                },
                primary_sp4,
                number4.width(),
                primary_row,
                ArrowEndpointDirection::Backward,
            ));
        }

        numbers.extend([number1, number2, number3, number4]);
    }

    debug!("Rendering TS arrangement");
    let mut ts_group = Group::new();
    let mut label_group = Group::new();

    let reference = IndexedStr::new(reference);
    let query = IndexedStr::new(query);
    let reference_c = IndexedStr::new(&reference_c);
    let query_c = IndexedStr::new(&query_c);

    let mut rows = HashMap::new();
    let mut y = 0.0;
    let mut ts_label_group_width = 0.0f32;

    // Reference complement inners.
    for (identifier, reference_inner) in ts_arrangement.reference_complement_inners().rev() {
        ts_group = ts_group.add(svg_string(
            reference_inner
                .sequence()
                .iter_values()
                .map(render_inner_char(
                    match reference_inner.template_switch().primary {
                        TemplateSwitchPrimary::Reference => &reference,
                        TemplateSwitchPrimary::Query => &query,
                    },
                )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));

        let label = format!(
            "TS-{} inner:",
            TS_RUNNING_NUMBER
                .chars()
                .nth(reference_inner.template_switch().index)
                .unwrap()
        );
        ts_label_group_width = ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        label_group = label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));

        rows.insert(TsArrangementRow::Inner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    // Reference complement.
    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference_complement()
            .iter_values()
            .map(render_complement_char(&reference_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));

    let label = "Reference complement:";
    ts_label_group_width = ts_label_group_width
        .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    label_group = label_group.add(svg_string(
        label.chars().map(render_label_char),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));

    rows.insert(TsArrangementRow::ReferenceComplement, y);
    y += typewriter::FONT.character_height;

    // Reference inners.
    for (identifier, reference_inner) in ts_arrangement.reference_inners().rev() {
        ts_group = ts_group.add(svg_string(
            reference_inner
                .sequence()
                .iter_values()
                .map(render_inner_char(
                    match reference_inner.template_switch().primary {
                        TemplateSwitchPrimary::Reference => &reference,
                        TemplateSwitchPrimary::Query => &query,
                    },
                )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));

        let label = format!(
            "TS-{} inner:",
            TS_RUNNING_NUMBER
                .chars()
                .nth(reference_inner.template_switch().index)
                .unwrap()
        );
        ts_label_group_width = ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        label_group = label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));

        rows.insert(TsArrangementRow::Inner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    // Reference.
    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference()
            .iter_values()
            .map(render_source_char(&reference)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));

    let label = "Reference:";
    ts_label_group_width = ts_label_group_width
        .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    label_group = label_group.add(svg_string(
        label.chars().map(render_label_char),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));

    rows.insert(TsArrangementRow::Reference, y);
    y += typewriter::FONT.character_height;

    // Query.
    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query()
            .iter_values()
            .map(render_source_char(&query)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));

    let label = "Query:";
    ts_label_group_width = ts_label_group_width
        .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    label_group = label_group.add(svg_string(
        label.chars().map(render_label_char),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));

    rows.insert(TsArrangementRow::Query, y);
    y += typewriter::FONT.character_height;

    // Query inners.
    for (identifier, query_inner) in ts_arrangement.query_inners() {
        ts_group = ts_group.add(svg_string(
            query_inner.sequence().iter_values().map(render_inner_char(
                match query_inner.template_switch().primary {
                    TemplateSwitchPrimary::Reference => &reference,
                    TemplateSwitchPrimary::Query => &query,
                },
            )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));

        let label = format!(
            "TS-{} inner:",
            TS_RUNNING_NUMBER
                .chars()
                .nth(query_inner.template_switch().index)
                .unwrap()
        );
        ts_label_group_width = ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        label_group = label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));

        rows.insert(TsArrangementRow::Inner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    // Query complement.
    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query_complement()
            .iter_values()
            .map(render_complement_char(&query_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));

    let label = "Query complement:";
    ts_label_group_width = ts_label_group_width
        .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    label_group = label_group.add(svg_string(
        label.chars().map(render_label_char),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));

    rows.insert(TsArrangementRow::QueryComplement, y);
    y += typewriter::FONT.character_height;

    // Query complement inners.
    for (identifier, query_inner) in ts_arrangement.query_complement_inners() {
        ts_group = ts_group.add(svg_string(
            query_inner.sequence().iter_values().map(render_inner_char(
                match query_inner.template_switch().primary {
                    TemplateSwitchPrimary::Reference => &reference,
                    TemplateSwitchPrimary::Query => &query,
                },
            )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));

        let label = format!(
            "TS-{} inner:",
            TS_RUNNING_NUMBER
                .chars()
                .nth(query_inner.template_switch().index)
                .unwrap()
        );
        ts_label_group_width = ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        label_group = label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));

        rows.insert(TsArrangementRow::Inner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    debug!("Rendering arrows");
    for arrow in &arrows {
        ts_group = ts_group.add(arrow.render(|row| *rows.get(row).unwrap()));
    }

    debug!("Rendering numbers");
    for number in &numbers {
        ts_group = ts_group.add(number.render(|row| *rows.get(row).unwrap()));
    }

    let ts_label_group_width = ts_label_group_width + sans_serif_mono::FONT.character_width;
    let ts_label_group =
        label_group.set("transform", SvgLocation { x: 0.0, y: 0.0 }.as_transform());

    let ts_group_width = ts_arrangement.width() as f32 * typewriter::FONT.character_width;
    let ts_group_height = y;

    let mut body_group = Group::new().set(
        "transform",
        SvgLocation {
            x: SVG_PADDING,
            y: SVG_PADDING,
        }
        .as_transform(),
    );
    body_group = body_group.add(ts_label_group);

    let (body_group_width, body_group_height, label_group_width) = if let Some(no_ts_arrangement) =
        no_ts_arrangement
    {
        debug!("Rendering no-TS arrangement");

        let vertical_spacer_height = typewriter::FONT.character_height;
        let mut no_ts_group = Group::new();
        let mut no_ts_label_group = Group::new();

        let mut y = 0.0;
        let mut no_ts_label_group_width = 0.0f32;

        no_ts_group = no_ts_group.add(svg_string(
            no_ts_arrangement
                .reference()
                .iter_values()
                .map(render_source_char(&reference)),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        let label = "Reference:";
        no_ts_label_group_width = no_ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        no_ts_label_group = no_ts_label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));
        y += typewriter::FONT.character_height;

        no_ts_group = no_ts_group.add(svg_string(
            no_ts_arrangement
                .query()
                .iter_values()
                .map(render_source_char(&query)),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        let label = "Query:";
        no_ts_label_group_width = no_ts_label_group_width
            .max(label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
        no_ts_label_group = no_ts_label_group.add(svg_string(
            label.chars().map(render_label_char),
            &SvgLocation { x: 0.0, y },
            &sans_serif_mono::FONT,
        ));
        y += typewriter::FONT.character_height;

        let no_ts_label_group_width =
            no_ts_label_group_width + sans_serif_mono::FONT.character_width;
        let no_ts_group_width = no_ts_arrangement.width() as f32 * typewriter::FONT.character_width;
        let no_ts_group_height = y;
        let label_group_width = (ts_label_group_width).max(no_ts_label_group_width);

        let no_ts_label_group = no_ts_label_group.set(
            "transform",
            SvgLocation {
                x: 0.0,
                y: ts_group_height + vertical_spacer_height,
            }
            .as_transform(),
        );
        let no_ts_group = no_ts_group.set(
            "transform",
            SvgLocation {
                x: label_group_width,
                y: ts_group_height + vertical_spacer_height,
            }
            .as_transform(),
        );

        body_group = body_group.add(no_ts_label_group).add(no_ts_group);

        (
            label_group_width + ts_group_width.max(no_ts_group_width),
            ts_group_height + vertical_spacer_height + no_ts_group_height,
            label_group_width,
        )
    } else {
        (
            ts_label_group_width + ts_group_width,
            ts_group_height,
            ts_label_group_width,
        )
    };

    let ts_group = ts_group.set(
        "transform",
        SvgLocation {
            x: label_group_width,
            y: 0.0,
        }
        .as_transform(),
    );
    body_group = body_group.add(ts_group);

    let (legend_group, legend_width, legend_height) = legend(
        &statistics.sequences.reference_name,
        &statistics.sequences.query_name,
        0.6,
    );
    let vertical_spacer_height = typewriter::FONT.character_height;
    body_group = body_group.add(
        legend_group.set(
            "transform",
            SvgLocation {
                x: 0.0,
                y: body_group_height + vertical_spacer_height,
            }
            .as_transform(),
        ),
    );
    let body_group_width = body_group_width.max(legend_width);
    let body_group_height = body_group_height + vertical_spacer_height + legend_height;

    debug!("Creating SVG root");
    let mut svg = Document::new();
    svg = add_arrow_defs(svg);
    svg = svg.add(Circle::new().set("r", 1e5).set("fill", "white"));
    svg = svg.set(
        "viewBox",
        (
            0,
            0,
            body_group_width + 2.0 * SVG_PADDING,
            body_group_height + 2.0 * SVG_PADDING,
        ),
    );
    svg = svg.add(body_group);

    debug!("Storing SVG");
    svg::write(output, &svg)?;

    Ok(())
}

fn render_source_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&SourceChar) -> Character<CharacterData> {
    |source_char| match source_char {
        SourceChar::Source {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Separator => Character::new_char_with_default('|'),
        SourceChar::Hidden { .. } | SourceChar::Spacer | SourceChar::Blank => {
            Character::new_char_with_default(' ')
        }
    }
}

fn render_complement_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&ComplementChar) -> Character<CharacterData> {
    |complement_char| match complement_char {
        ComplementChar::Complement {
            column,
            lower_case,
            source_hidden,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(
                c,
                if *source_hidden {
                    CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
                } else {
                    Default::default()
                },
            )
        }
        ComplementChar::Gap { source_hidden } => Character::new_char(
            '-',
            if *source_hidden {
                CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
            } else {
                Default::default()
            },
        ),
        ComplementChar::Hidden { .. } | ComplementChar::Blank => {
            Character::new_char(' ', Default::default())
        }
    }
}

fn render_inner_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&InnerChar) -> Character<CharacterData> {
    |inner_char| match inner_char {
        InnerChar::Inner {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };
            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Blank => Character::new_char(' ', Default::default()),
    }
}

fn render_label_char(c: char) -> Character<CharacterData> {
    Character::new_char(c, CharacterData::new_colored("#555555"))
}

fn copy_color(copy_depth: &Option<usize>) -> impl ToString {
    if let Some(copy_depth) = copy_depth {
        COPY_COLORS[copy_depth % COPY_COLORS.len()]
    } else {
        "black"
    }
}

/// Returns the SVG legend along with its width and height.
fn legend(reference_name: &str, query_name: &str, scale: f32) -> (Group, f32, f32) {
    let mut result = Group::new();

    let headline = "Legend:";
    result = result.add(svg_string(
        headline
            .chars()
            .map(Character::<CharacterData>::new_char_with_default),
        &SvgLocation { x: 0.0, y: 0.0 },
        &sans_serif_mono::FONT,
    ));
    let headline_width = headline.chars().count() as f32 * sans_serif_mono::FONT.character_width;
    let headline_height = sans_serif_mono::FONT.character_height * 1.4;

    let mut label_width = 0.0f32;
    let mut explanation_width = 0.0f32;
    let mut y = headline_height;

    // Labels.
    let reference_label = "Reference";
    result = result.add(svg_string(
        reference_label
            .chars()
            .map(Character::<CharacterData>::new_char_with_default),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));
    label_width = label_width
        .max(reference_label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    y += typewriter::FONT
        .character_height
        .max(sans_serif_mono::FONT.character_height);

    let query_label = "Query";
    result = result.add(svg_string(
        query_label
            .chars()
            .map(Character::<CharacterData>::new_char_with_default),
        &SvgLocation { x: 0.0, y },
        &sans_serif_mono::FONT,
    ));
    label_width =
        label_width.max(query_label.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    y += typewriter::FONT
        .character_height
        .max(sans_serif_mono::FONT.character_height);

    let copy_label = "GREEN CHARACTERS";
    result = result.add(svg_string(
        copy_label
            .chars()
            .map(|c| Character::new_char(c, CharacterData::new_colored(COPY_COLORS[0]))),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    label_width =
        label_width.max(copy_label.chars().count() as f32 * typewriter::FONT.character_width);

    // Explanations.
    let label_width = label_width + typewriter::FONT.character_width;
    y = headline_height;

    result = result.add(svg_text(reference_name, &SvgLocation { x: label_width, y }));
    explanation_width = explanation_width
        .max(reference_name.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    y += typewriter::FONT
        .character_height
        .max(sans_serif_mono::FONT.character_height);

    result = result.add(svg_text(query_name, &SvgLocation { x: label_width, y }));
    explanation_width = explanation_width
        .max(query_name.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    y += typewriter::FONT
        .character_height
        .max(sans_serif_mono::FONT.character_height);

    let copy_explanation = "Repeated characters due to a TS with SP4 < SP1";
    result = result.add(svg_string(
        copy_explanation
            .chars()
            .map(Character::<CharacterData>::new_char_with_default),
        &SvgLocation { x: label_width, y },
        &sans_serif_mono::FONT,
    ));
    explanation_width = explanation_width
        .max(copy_explanation.chars().count() as f32 * sans_serif_mono::FONT.character_width);
    y += typewriter::FONT
        .character_height
        .max(sans_serif_mono::FONT.character_height);

    result = result.set("transform", format!("scale({scale})"));
    let legend_width = headline_width.max(label_width + explanation_width) * scale;
    let legend_height = y * scale;

    (Group::new().add(result), legend_width, legend_height)
}

pub fn create_error_svg(output: impl Write, error: Error) -> Result<()> {
    let svg = Document::new()
        .set("viewBox", (0, 0, 1920, 1080))
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(
            Text::new(error.to_string())
                .set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform())
                .set("font", "Sans serif"),
        );

    svg::write(output, &svg)?;
    Ok(())
}
