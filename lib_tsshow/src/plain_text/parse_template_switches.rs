use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::{
        Alignment,
        iter::CompactAlignmentIterCloned,
        stream::{AlignmentStream, AlignmentStreamCoordinates},
    },
    template_switch_distance::{
        AlignmentType, TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection,
    },
};
use log::{debug, trace};

use crate::error::Result;

pub const STREAM_DEFAULT_LENGTH: usize = 20;
pub const STREAM_PADDING: usize = 10;

#[derive(Debug)]
pub struct TSShow<AlignmentType> {
    pub upstream_offset: AlignmentStreamCoordinates,
    pub downstream_limit: AlignmentStreamCoordinates,
    pub sp1_offset: AlignmentStreamCoordinates,
    pub sp2_ancestor_offset: usize,
    pub sp3_ancestor_offset: usize,
    pub sp4_offset: AlignmentStreamCoordinates,
    pub descendant: TemplateSwitchDescendant,
    pub ancestor: TemplateSwitchAncestor,
    pub upstream: Alignment<AlignmentType>,
    pub template_switch: Alignment<AlignmentType>,
    pub downstream: Alignment<AlignmentType>,
}

pub fn parse(
    alignment: &Alignment<AlignmentType>,
    reference_offset: usize,
    query_offset: usize,
) -> Result<Vec<TSShow<AlignmentType>>> {
    let mut template_switches = Vec::new();
    let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
    let mut alignment = alignment.iter_compact_cloned();

    while let Some((multiplicity, alignment_type)) = alignment.peek_front_cloned() {
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            template_switches.push(parse_template_switch(&mut alignment, &mut stream)?);
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push(multiplicity, alignment_type);
            alignment.next();
        }
    }

    Ok(template_switches)
}

fn parse_template_switch(
    alignment: &mut CompactAlignmentIterCloned<AlignmentType>,
    stream: &mut AlignmentStream,
) -> Result<TSShow<AlignmentType>> {
    let (
        multiplicity,
        alignment_type @ AlignmentType::TemplateSwitchEntrance {
            descendant,
            ancestor,
            direction,
            first_offset,
            ..
        },
    ) = alignment.peek_front_cloned().unwrap()
    else {
        unreachable!("Function is only called with a template switch entrance")
    };
    debug!("Parsing TS with first_offset {first_offset}");

    let sp1_offset = stream.head_coordinates();
    let upstream = stream.clone();
    let mut template_switch = Vec::new();

    stream.push(multiplicity, alignment_type);
    alignment.next().unwrap();

    let sp2_ancestor_offset: usize = match ancestor {
        TemplateSwitchAncestor::Reference => (sp1_offset.reference() as isize + first_offset)
            .try_into()
            .unwrap(),
        TemplateSwitchAncestor::Query => (sp1_offset.query() as isize + first_offset)
            .try_into()
            .unwrap(),
    };
    let mut sp3_ancestor_offset = sp2_ancestor_offset;

    while let Some((multiplicity, alignment_type)) = alignment.next() {
        debug!("alignment type: {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            panic!("Found template switch entrance within template switch");
        } else if let AlignmentType::TemplateSwitchExit { .. } = alignment_type {
            trace!("Found TS exit");
            stream.push(multiplicity, alignment_type);

            trace!("Creating upstream");
            let mut upstream = upstream;
            upstream.pop(
                STREAM_DEFAULT_LENGTH.max(
                    sp1_offset
                        .reference()
                        .max(sp1_offset.query())
                        .saturating_sub(sp2_ancestor_offset.min(sp3_ancestor_offset))
                        + STREAM_PADDING,
                ),
            );
            let upstream_offset = upstream.tail_coordinates();
            let upstream = upstream.stream_alignment();

            let template_switch = template_switch.into();

            trace!("Creating downstream");
            stream.clear();
            let sp4_offset = stream.head_coordinates();
            let downstream = parse_downstream(
                alignment,
                stream,
                STREAM_DEFAULT_LENGTH.max(
                    sp3_ancestor_offset.max(sp2_ancestor_offset).saturating_sub(
                        sp4_offset.reference().min(sp4_offset.query()) + STREAM_PADDING,
                    ),
                ),
            );
            let downstream_limit = stream.head_coordinates();

            trace!("Returning TSShow");
            return Ok(TSShow {
                upstream_offset,
                downstream_limit,
                sp1_offset,
                sp2_ancestor_offset,
                sp3_ancestor_offset,
                sp4_offset,
                descendant,
                ancestor,
                upstream,
                template_switch,
                downstream,
            });
        } else {
            template_switch.push((multiplicity, alignment_type));
            stream.push(multiplicity, alignment_type);

            if matches!(
                alignment_type,
                AlignmentType::SecondaryDeletion
                    | AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch
            ) {
                match direction {
                    TemplateSwitchDirection::Forward => {
                        sp3_ancestor_offset += multiplicity;
                        assert!(sp3_ancestor_offset > sp2_ancestor_offset);
                    }
                    TemplateSwitchDirection::Reverse => {
                        sp3_ancestor_offset -= multiplicity;
                        assert!(sp3_ancestor_offset < sp2_ancestor_offset);
                    }
                }
            }
        }
    }

    panic!("Found template switch without exit")
}

fn parse_downstream(
    alignment: &mut CompactAlignmentIterCloned<AlignmentType>,
    stream: &mut AlignmentStream,
    requested_length: usize,
) -> Alignment<AlignmentType> {
    debug!("Parsing downstream");

    stream.clear();

    while let Some((_, alignment_type)) = alignment.peek_front_cloned() {
        trace!("alignment_type: {alignment_type}");
        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            break;
        } else if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
            panic!("Found template switch exit without matching entrance");
        } else {
            stream.push_until_full(
                &mut alignment.peek_front_multiplicity_mut().unwrap(),
                alignment_type,
                requested_length,
            );

            if stream.is_full(requested_length) {
                trace!("Stream is full");
                break;
            }
        }
    }

    stream.stream_alignment()
}
