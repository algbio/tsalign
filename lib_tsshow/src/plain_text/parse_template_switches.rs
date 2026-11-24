use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_result::{
        AlignmentType,
        alignment::{
            Alignment,
            iter::CompactAlignmentIterCloned,
            stream::{AlignmentStream, AlignmentStreamCoordinates},
        },
    },
};
use log::{debug, trace};

use crate::error::Result;

pub const STREAM_DEFAULT_LENGTH: usize = 20;
pub const STREAM_PADDING: usize = 10;

#[derive(Debug)]
pub struct TSShow {
    pub upstream_offset: AlignmentStreamCoordinates,
    pub downstream_limit: AlignmentStreamCoordinates,
    pub sp1_offset: AlignmentStreamCoordinates,
    pub sp2_ancestor_offset: usize,
    pub sp3_ancestor_offset: usize,
    pub sp4_offset: AlignmentStreamCoordinates,
    pub descendant: TemplateSwitchDescendant,
    pub ancestor: TemplateSwitchAncestor,
    pub upstream: Alignment,
    pub template_switch: Alignment,
    pub downstream: Alignment,
}

pub fn parse(
    alignment: &Alignment,
    reference_offset: usize,
    query_offset: usize,
) -> Result<Vec<TSShow>> {
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
    alignment: &mut CompactAlignmentIterCloned<'_>,
    stream: &mut AlignmentStream,
) -> Result<TSShow> {
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
    alignment: &mut CompactAlignmentIterCloned<'_>,
    stream: &mut AlignmentStream,
    requested_length: usize,
) -> Alignment {
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

#[cfg(test)]
mod tests {
    use lib_tsalign::{
        aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
        alignment_result::{AlignmentType, TSMUncertaintyRange, alignment::Alignment},
    };

    use super::parse;

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

    #[test]
    fn reverse_template_switch() {
        let alignment = Alignment::from(vec![
            (4, AlignmentType::PrimaryMatch),
            (
                1,
                entrance(
                    3,
                    TemplateSwitchDescendant::Query,
                    TemplateSwitchAncestor::Reference,
                    TemplateSwitchDirection::Reverse,
                ),
            ),
            (2, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (1, exit(3)),
            (3, AlignmentType::PrimaryMatch),
        ]);
        let template_switches = parse(&alignment, 0, 0).unwrap();
        assert_eq!(template_switches.len(), 1);
        let ts = &template_switches[0];

        assert_eq!(ts.descendant, TemplateSwitchDescendant::Query);
        assert_eq!(ts.ancestor, TemplateSwitchAncestor::Reference);
        assert_eq!((ts.sp1_offset.reference(), ts.sp1_offset.query()), (4, 4));
        // SP2 = SP1 on the ancestor + first offset, and the reverse inner goes left.
        assert_eq!(ts.sp2_ancestor_offset, 7);
        assert_eq!(ts.sp3_ancestor_offset, 4);
        assert_eq!((ts.sp4_offset.reference(), ts.sp4_offset.query()), (7, 7));
        assert_eq!(
            (ts.upstream_offset.reference(), ts.upstream_offset.query()),
            (0, 0)
        );
        assert_eq!(
            (ts.downstream_limit.reference(), ts.downstream_limit.query()),
            (10, 10)
        );
        assert_eq!(ts.upstream.cigar(), "4=");
        assert_eq!(ts.template_switch.cigar(), "2=1X");
        assert_eq!(ts.downstream.cigar(), "3=");
    }

    #[test]
    fn forward_template_switch_with_offsets() {
        let alignment = Alignment::from(vec![
            (2, AlignmentType::PrimaryMatch),
            (
                1,
                entrance(
                    -2,
                    TemplateSwitchDescendant::Reference,
                    TemplateSwitchAncestor::Query,
                    TemplateSwitchDirection::Forward,
                ),
            ),
            (3, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondaryDeletion),
            (1, exit(-1)),
            (2, AlignmentType::PrimaryMatch),
        ]);
        let template_switches = parse(&alignment, 5, 10).unwrap();
        assert_eq!(template_switches.len(), 1);
        let ts = &template_switches[0];

        assert_eq!((ts.sp1_offset.reference(), ts.sp1_offset.query()), (7, 12));
        assert_eq!(ts.sp2_ancestor_offset, 10);
        assert_eq!(ts.sp3_ancestor_offset, 14);
        // The descendant advances by the inner without the deletion, the anti-descendant by the gap.
        assert_eq!((ts.sp4_offset.reference(), ts.sp4_offset.query()), (10, 11));
    }

    #[test]
    fn several_template_switches() {
        let ts = |descendant| {
            [
                (
                    1,
                    entrance(
                        2,
                        descendant,
                        TemplateSwitchAncestor::Reference,
                        TemplateSwitchDirection::Reverse,
                    ),
                ),
                (2, AlignmentType::SecondaryMatch),
                (1, exit(2)),
            ]
        };
        let mut alignment = vec![(3, AlignmentType::PrimaryMatch)];
        alignment.extend(ts(TemplateSwitchDescendant::Query));
        alignment.push((1, AlignmentType::PrimaryInsertion));
        alignment.extend(ts(TemplateSwitchDescendant::Reference));
        alignment.push((3, AlignmentType::PrimaryMatch));
        let template_switches = parse(&Alignment::from(alignment), 0, 0).unwrap();

        assert_eq!(template_switches.len(), 2);
        assert_eq!(
            template_switches[0].descendant,
            TemplateSwitchDescendant::Query
        );
        assert_eq!(
            template_switches[1].descendant,
            TemplateSwitchDescendant::Reference
        );
        assert_eq!(
            (
                template_switches[1].sp1_offset.reference(),
                template_switches[1].sp1_offset.query()
            ),
            (5, 6)
        );
        // The downstream of the first TS ends at the next TS.
        assert_eq!(template_switches[0].downstream.cigar(), "1I");
    }

    #[test]
    fn no_template_switches() {
        let alignment = Alignment::from(vec![(5, AlignmentType::PrimaryMatch)]);
        assert!(parse(&alignment, 0, 0).unwrap().is_empty());
    }
}
