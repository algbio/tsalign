use std::{io::Write, iter};

use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{
            AlignmentResult,
            a_star_sequences::SequencePair,
            alignment::stream::{AlignmentStream, AlignmentStreamCoordinates},
        },
        template_switch_distance::{
            AlignmentType, TemplateSwitchAncestor, TemplateSwitchDescendant,
        },
    },
    costs::U64Cost,
};
use log::{debug, info, trace, warn};
use mutlipair_alignment_renderer::MultipairAlignmentRenderer;
use parse_template_switches::{STREAM_PADDING, TSShow};

use crate::error::Result;

pub mod mutlipair_alignment_renderer;
mod parse_template_switches;

pub fn show_template_switches(
    mut output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
) -> Result<()> {
    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        warn!("Alignment was aborted early, no template switches present");
        return Ok(());
    };

    info!("CIGAR: {} (Cost: {:.0})", result.cigar(), statistics.cost);
    if let Some(no_ts_result) = no_ts_result.as_ref() {
        info!(
            "No-ts CIGAR: {} (Cost: {:.0})",
            no_ts_result.cigar(),
            no_ts_result.statistics().cost
        )
    }

    debug!("Collecting template switches");
    let template_switches = parse_template_switches::parse(
        alignment,
        statistics.reference_offset,
        statistics.query_offset,
    )?;
    info!("Found {} template switches", template_switches.len());

    for (index, template_switch) in template_switches.iter().enumerate() {
        info!("Showing template switch {}", index + 1);
        show_template_switch(
            &mut output,
            template_switch,
            &statistics.sequences,
            no_ts_result,
            statistics.reference_offset,
            statistics.query_offset,
        );
    }

    Ok(())
}

fn show_template_switch(
    mut output: impl Write,
    template_switch: &TSShow<AlignmentType>,
    sequences: &SequencePair,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    reference_offset: usize,
    query_offset: usize,
) {
    trace!("Showing template switch\n{template_switch:?}");

    let forward = template_switch.sp2_ancestor_offset < template_switch.sp3_ancestor_offset;
    let reference = &sequences.reference;
    let reference_c: String = sequences.reference_rc.chars().rev().collect();
    let query = &sequences.query;
    let query_c: String = sequences.query_rc.chars().rev().collect();

    let (
        primary_label,
        primary_name,
        primary,
        primary_c,
        primary_coordinate_picker,
        anti_primary_label,
        anti_primary_name,
        anti_primary,
        anti_primary_c,
        anti_primary_coordinate_picker,
        invert_alignment,
    ) = match template_switch.descendant {
        TemplateSwitchDescendant::Reference => (
            "Parent".to_string(),
            &sequences.reference_name,
            reference,
            &reference_c,
            Box::new(|alignment_coordinates: &AlignmentStreamCoordinates| {
                alignment_coordinates.reference()
            }) as Box<dyn Fn(&AlignmentStreamCoordinates) -> usize>,
            "Child".to_string(),
            &sequences.query_name,
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentStreamCoordinates| {
                alignment_coordinates.query()
            }) as Box<dyn Fn(&AlignmentStreamCoordinates) -> usize>,
            true,
        ),
        TemplateSwitchDescendant::Query => (
            "Child".to_string(),
            &sequences.query_name,
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentStreamCoordinates| {
                alignment_coordinates.query()
            }) as Box<dyn Fn(&AlignmentStreamCoordinates) -> usize>,
            "Parent".to_string(),
            &sequences.reference_name,
            reference,
            &reference_c,
            Box::new(|alignment_coordinates: &AlignmentStreamCoordinates| {
                alignment_coordinates.reference()
            }) as Box<dyn Fn(&AlignmentStreamCoordinates) -> usize>,
            false,
        ),
    };
    let primary_equals_secondary = (template_switch.descendant
        == TemplateSwitchDescendant::Reference
        && template_switch.ancestor == TemplateSwitchAncestor::Reference)
        || (template_switch.descendant == TemplateSwitchDescendant::Query
            && template_switch.ancestor == TemplateSwitchAncestor::Query);

    let primary_forward_label = format!("{primary_label}F");
    let primary_reverse_label = format!("{primary_label}R");
    let anti_primary_forward_label = format!("{anti_primary_label}F");
    let anti_primary_reverse_label = format!("{anti_primary_label}R");
    let f1_label = format!("{primary_label}1");
    let f2_label = format!("{primary_label}2");
    let f3_label = format!("{primary_label}3");

    let primary_offset = primary_coordinate_picker(&template_switch.upstream_offset);
    let primary_limit = primary_coordinate_picker(&template_switch.downstream_limit);

    let anti_primary_f1_offset = anti_primary_coordinate_picker(&template_switch.upstream_offset);
    let anti_primary_f3_offset = anti_primary_coordinate_picker(&template_switch.sp4_offset);
    let anti_primary_offset = anti_primary_f1_offset.min(anti_primary_f3_offset);
    let anti_primary_f1_limit = anti_primary_coordinate_picker(&template_switch.sp1_offset);
    let anti_primary_f3_limit = anti_primary_coordinate_picker(&template_switch.downstream_limit);
    let anti_primary_limit = anti_primary_f1_limit.max(anti_primary_f3_limit);

    let primary_sp1 = primary_coordinate_picker(&template_switch.sp1_offset);
    let anti_primary_sp1 = anti_primary_coordinate_picker(&template_switch.sp1_offset);
    let primary_sp4 = primary_coordinate_picker(&template_switch.sp4_offset);
    let anti_primary_sp4 = anti_primary_coordinate_picker(&template_switch.sp4_offset);

    let (ts_inner, ts_inner_alignment) = if forward {
        (
            primary[primary_coordinate_picker(&template_switch.sp1_offset)
                ..primary_coordinate_picker(&template_switch.sp4_offset)]
                .to_string(),
            template_switch.template_switch.clone(),
        )
    } else {
        (
            primary[primary_coordinate_picker(&template_switch.sp1_offset)
                ..primary_coordinate_picker(&template_switch.sp4_offset)]
                .chars()
                .rev()
                .collect(),
            template_switch.template_switch.reverse(),
        )
    };

    debug!("Primary offset: {primary_offset}");
    debug!("SP1 primary offset: {primary_sp1}");
    debug!("SP4 primary offset: {primary_sp4}");
    debug!("Primary limit: {primary_limit}");
    debug!("Primary len: {}", primary.len());
    debug!("Anti-primary offset: {anti_primary_offset}");
    debug!("SP1 anti-primary offset: {anti_primary_sp1}");
    debug!("SP4 anti-primary offset: {anti_primary_sp4}");
    debug!("Anti-primary limit: {anti_primary_f3_limit}");
    debug!("Anti-primary len: {}", anti_primary.len());
    debug!(
        "SP2 secondary offset: {}",
        template_switch.sp2_ancestor_offset
    );
    debug!(
        "SP3 secondary offset: {}",
        template_switch.sp3_ancestor_offset
    );

    if primary_equals_secondary {
        let primary_extended_offset = primary_offset.min(
            template_switch
                .sp3_ancestor_offset
                .min(template_switch.sp2_ancestor_offset)
                .saturating_sub(STREAM_PADDING),
        );
        let primary_extended_limit = primary_limit.max(
            primary.chars().count().min(
                template_switch
                    .sp2_ancestor_offset
                    .max(template_switch.sp3_ancestor_offset)
                    + STREAM_PADDING,
            ),
        );

        debug!("Primary extended offset: {primary_extended_offset}");
        debug!("Primary extended limit: {primary_extended_limit}");

        debug!("Creating outside renderer");
        let mut outside_renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_forward_label.clone(),
            anti_primary[anti_primary_offset..anti_primary_limit].chars(),
        );

        debug!("Adding F1");
        let f1_sequence =
            &primary[primary_offset..primary_coordinate_picker(&template_switch.sp1_offset)];
        trace!("Current renderer content:\n{}", {
            let mut out = Vec::new();
            outside_renderer
                .render_without_names(&mut out, [&anti_primary_forward_label])
                .unwrap();
            String::from_utf8(out).unwrap()
        });
        trace!("f1_sequence: {f1_sequence}");
        trace!(
            "f1_alignment: {} ({:?})",
            template_switch.upstream.cigar(),
            template_switch.upstream
        );

        outside_renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_f1_offset
                .checked_sub(anti_primary_offset)
                .unwrap(),
            f1_label.clone(),
            f1_sequence.chars(),
            template_switch.upstream.iter_flat_cloned(),
            true,
            invert_alignment,
        );
        debug!("Adding F3");
        outside_renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_f3_offset
                .checked_sub(anti_primary_offset)
                .unwrap(),
            f3_label.clone(),
            primary[primary_coordinate_picker(&template_switch.sp4_offset)..primary_limit].chars(),
            template_switch.downstream.iter_flat_cloned(),
            true,
            invert_alignment,
        );

        debug!("Creating inside renderer");
        let mut inside_renderer = if forward {
            MultipairAlignmentRenderer::new_without_data(
                primary_forward_label.clone(),
                primary[primary_extended_offset..primary_extended_limit].chars(),
            )
        } else {
            MultipairAlignmentRenderer::new_without_data(
                primary_reverse_label.clone(),
                primary_c[primary_extended_offset..primary_extended_limit].chars(),
            )
        };
        debug!("Adding F2");
        inside_renderer.add_aligned_sequence_without_data(
            if forward {
                &primary_forward_label
            } else {
                &primary_reverse_label
            },
            template_switch
                .sp3_ancestor_offset
                .min(template_switch.sp2_ancestor_offset)
                - primary_extended_offset,
            f2_label.clone(),
            ts_inner.chars(),
            ts_inner_alignment.iter_flat_cloned(),
            true,
            false,
        );

        println!("{anti_primary_label}: {anti_primary_name}");
        println!("{primary_label}: {primary_name}");
        println!("Direction: {}", if forward { "forward" } else { "reverse" });
        println!();
        println!("Switch process:");

        debug!("Rendering");
        outside_renderer
            .render(
                &mut output,
                [&f1_label, &f3_label, &anti_primary_forward_label],
            )
            .unwrap();
        println!();
        inside_renderer
            .render(
                &mut output,
                [
                    if forward {
                        &primary_forward_label
                    } else {
                        &primary_reverse_label
                    },
                    &f2_label,
                ],
            )
            .unwrap();
    } else {
        let anti_primary_extended_offset = anti_primary_offset.min(
            template_switch
                .sp3_ancestor_offset
                .min(template_switch.sp2_ancestor_offset)
                .saturating_sub(STREAM_PADDING),
        );
        let anti_primary_extended_limit = anti_primary_f3_limit.max(
            anti_primary.chars().count().min(
                template_switch
                    .sp2_ancestor_offset
                    .max(template_switch.sp3_ancestor_offset)
                    + STREAM_PADDING,
            ),
        );

        debug!("Anti-primary extended offset: {anti_primary_extended_offset}");
        debug!("Anti-primary extended limit: {anti_primary_extended_limit}");

        debug!("Creating renderer");
        let mut renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_forward_label.clone(),
            anti_primary[anti_primary_extended_offset..anti_primary_extended_limit].chars(),
        );

        if !forward {
            debug!("Adding complement primary");
            renderer.add_aligned_sequence_without_data(
                &anti_primary_forward_label,
                0,
                anti_primary_reverse_label.clone(),
                anti_primary_c[anti_primary_extended_offset..anti_primary_extended_limit].chars(),
                iter::repeat_n(
                    AlignmentType::PrimaryMatch,
                    anti_primary_extended_limit - anti_primary_extended_offset,
                ),
                false,
                false,
            );
        }

        debug!("Adding F1");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_offset - anti_primary_extended_offset,
            f1_label.clone(),
            primary[primary_offset..primary_coordinate_picker(&template_switch.sp1_offset)].chars(),
            template_switch.upstream.iter_flat_cloned(),
            true,
            invert_alignment,
        );
        debug!("Adding F3");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_coordinate_picker(&template_switch.sp4_offset)
                - anti_primary_extended_offset,
            f3_label.clone(),
            primary[primary_coordinate_picker(&template_switch.sp4_offset)..primary_limit].chars(),
            template_switch.downstream.iter_flat_cloned(),
            true,
            invert_alignment,
        );

        debug!("Adding F2");
        renderer.add_aligned_sequence_without_data(
            if forward {
                &anti_primary_forward_label
            } else {
                &anti_primary_reverse_label
            },
            template_switch
                .sp3_ancestor_offset
                .min(template_switch.sp2_ancestor_offset)
                - anti_primary_extended_offset,
            f2_label.clone(),
            ts_inner.chars(),
            ts_inner_alignment.iter_flat_cloned(),
            true,
            false,
        );

        println!("{anti_primary_label}: {anti_primary_name}");
        println!("{primary_label}: {primary_name}");
        println!("Direction: {}", if forward { "forward" } else { "reverse" });
        println!();
        println!("Switch process:");

        debug!("Rendering");
        renderer
            .render(
                &mut output,
                if forward {
                    vec![&f1_label, &f3_label, &anti_primary_forward_label, &f2_label]
                } else {
                    vec![
                        &f1_label,
                        &f3_label,
                        &anti_primary_forward_label,
                        &anti_primary_reverse_label,
                        &f2_label,
                    ]
                },
            )
            .unwrap();
    }

    if let Some(no_ts_result) = no_ts_result {
        let AlignmentResult::WithTarget {
            alignment: no_ts_alignment,
            ..
        } = no_ts_result
        else {
            warn!("No-ts alignment was aborted early, unable to render");
            return;
        };

        assert!(
            no_ts_alignment
                .iter_compact()
                .all(|(_, alignment_type)| !matches!(
                    alignment_type,
                    AlignmentType::TemplateSwitchEntrance { .. }
                )),
            "No-ts alignment must not contain template switches."
        );

        // Find subsequence of no-ts alignment that matches ts alignment interval.
        let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
        for alignment_type in no_ts_alignment.iter_flat_cloned() {
            if anti_primary_coordinate_picker(&stream.head_coordinates()) >= anti_primary_f3_limit {
                break;
            } else {
                stream.push(1, alignment_type);
            }
        }
        assert_eq!(
            anti_primary_coordinate_picker(&stream.head_coordinates()),
            anti_primary_f3_limit
        );

        while anti_primary_coordinate_picker(&stream.tail_coordinates()) < anti_primary_offset {
            stream.pop_one();
        }
        assert_eq!(
            anti_primary_coordinate_picker(&stream.tail_coordinates()),
            anti_primary_offset
        );

        debug!("Creating no-ts renderer");
        let mut renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_label.clone(),
            anti_primary[anti_primary_offset..anti_primary_f3_limit].chars(),
        );

        debug!("Adding primary");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_label,
            0,
            primary_label.clone(),
            primary[primary_coordinate_picker(&stream.tail_coordinates())
                ..primary_coordinate_picker(&stream.head_coordinates())]
                .chars(),
            stream.stream_iter_flat(),
            true,
            invert_alignment,
        );

        println!();
        println!("No-ts alignment:");

        debug!("Rendering");
        renderer
            .render(&mut output, [&anti_primary_label, &primary_label])
            .unwrap();
    } else {
        debug!("No no-ts alignment given, skipping");
    }
}
