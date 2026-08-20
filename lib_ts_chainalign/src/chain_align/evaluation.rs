use std::iter;

use generic_a_star::cost::AStarCost;
use log::trace;
use num_traits::Zero;

use crate::{
    alignment::{
        Alignment, AlignmentType,
        coordinates::{AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates},
        sequences::AlignmentSequences,
    },
    anchors::Anchors,
    chain_align::chainer::Identifier,
    chaining_cost_function::ChainingCostFunction,
    costs::AlignmentCosts,
    exact_chaining::{
        gap_affine::GapAffineAligner, ts_12_jump::Ts12JumpAligner, ts_34_jump::Ts34JumpAligner,
    },
    panic_on_extend::PanicOnExtend,
};

pub struct ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    sequences: &'sequences AlignmentSequences,
    primary_aligner: GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,
    secondary_aligner: GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,
    ts_12_jump_aligner: Ts12JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,
    ts_34_jump_aligner: Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,

    additional_primary_targets_buffer: Vec<(PrimaryAlignmentCoordinates, Cost)>,
    additional_secondary_targets_buffer: Vec<(AnySecondaryAlignmentCoordinates, Cost)>,

    total_gaps: u64,
    total_gap_fillings: u64,
    total_redundant_gap_fillings: u64,
    gap_fill_alignments_per_chain: Vec<u32>,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost>
    ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
    ) -> Self {
        Self {
            sequences,
            primary_aligner: GapAffineAligner::new(
                sequences,
                &alignment_costs.primary_costs,
                rc_fn,
                max_match_run,
            ),
            secondary_aligner: GapAffineAligner::new(
                sequences,
                &alignment_costs.secondary_costs,
                rc_fn,
                max_match_run,
            ),
            ts_12_jump_aligner: Ts12JumpAligner::new(
                sequences,
                alignment_costs,
                rc_fn,
                max_match_run,
            ),
            ts_34_jump_aligner: Ts34JumpAligner::new(
                sequences,
                alignment_costs,
                rc_fn,
                max_match_run,
            ),

            additional_primary_targets_buffer: Default::default(),
            additional_secondary_targets_buffer: Default::default(),

            total_gaps: 0,
            total_gap_fillings: 0,
            total_redundant_gap_fillings: 0,
            gap_fill_alignments_per_chain: Default::default(),
        }
    }

    pub fn evaluate_chain(
        &mut self,
        anchors: &Anchors<Cost>,
        chain: &[Identifier],
        max_match_run: u32,
        chaining_cost_function: &mut ChainingCostFunction<Cost>,
        final_evaluation: bool,
    ) -> (Cost, Vec<Alignment>) {
        let k = usize::try_from(max_match_run + 1).unwrap();

        let mut current_upper_bound = Cost::zero();
        let mut alignments = Vec::new();
        let mut current_from_index = 0;
        let mut gap_fill_alignment_count = 0;

        loop {
            let from_anchor = chain[current_from_index];
            let Some((to_anchor_index, to_anchor)) = chain
                .iter()
                .copied()
                .enumerate()
                .skip(current_from_index + 1)
                .find(|(_, identifier)| match identifier {
                    Identifier::Start => true,
                    Identifier::StartToPrimary { .. } => false,
                    Identifier::StartToSecondary { .. } => false,
                    Identifier::PrimaryToPrimary { offset, .. }
                    | Identifier::PrimaryToSecondary { offset, .. }
                    | Identifier::SecondaryToSecondary { offset, .. }
                    | Identifier::SecondaryToPrimary { offset, .. } => offset.is_zero(),
                    Identifier::End => true,
                })
            else {
                break;
            };
            current_from_index = to_anchor_index;
            self.total_gaps += 1;

            match (from_anchor, to_anchor) {
                (Identifier::Start, Identifier::End) => {
                    if final_evaluation || !chaining_cost_function.is_start_to_end_exact() {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.primary_aligner.align(
                            self.sequences.primary_start(),
                            self.sequences.primary_end(),
                            &mut self.additional_primary_targets_buffer,
                            &mut PanicOnExtend,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!("Aligning from start to end costs {}", cost);
                        if !final_evaluation {
                            chaining_cost_function.update_start_to_end(cost, true);
                            chaining_cost_function.update_additional_primary_targets_from_start(
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(alignment);
                    }
                    current_upper_bound =
                        current_upper_bound.saturating_add(&chaining_cost_function.start_to_end());
                }
                (
                    Identifier::Start,
                    Identifier::PrimaryToPrimary { index, .. }
                    | Identifier::PrimaryToSecondary { index, .. },
                ) => {
                    let end = anchors.primary(index).start();
                    if final_evaluation
                        || !chaining_cost_function.is_primary_from_start_exact(index)
                    {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.primary_aligner.align(
                            self.sequences.primary_start(),
                            end,
                            &mut self.additional_primary_targets_buffer,
                            &mut PanicOnExtend,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from start to P{index}{} costs {}",
                            anchors.primary(index),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function.update_primary_from_start(index, cost, true);
                            chaining_cost_function.update_additional_primary_targets_from_start(
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.primary_from_start(index));
                }
                (
                    Identifier::Start,
                    Identifier::SecondaryToPrimary { index, ts_kind, .. }
                    | Identifier::SecondaryToSecondary { index, ts_kind, .. },
                ) => {
                    let end = anchors
                        .secondary(index, ts_kind)
                        .start()
                        .into_specific(ts_kind);
                    if final_evaluation
                        || !chaining_cost_function.is_jump_12_from_start_exact(index, ts_kind)
                    {
                        self.additional_secondary_targets_buffer.clear();
                        let (cost, alignment) = self.ts_12_jump_aligner.align(
                            self.sequences.primary_start(),
                            end,
                            &mut self.additional_secondary_targets_buffer,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_secondary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from start to S{}[{index}]{} costs {}",
                            ts_kind.digits(),
                            anchors.secondary(index, ts_kind),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function
                                .update_jump_12_from_start(index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_12_jump_targets_from_start(
                                &mut self.additional_secondary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.jump_12_from_start(index, ts_kind));
                }
                (Identifier::PrimaryToPrimary { index, .. }, Identifier::End) => {
                    let start = anchors.primary(index).end(k);
                    if final_evaluation || !chaining_cost_function.is_primary_to_end_exact(index) {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.primary_aligner.align(
                            start,
                            self.sequences.primary_end(),
                            &mut self.additional_primary_targets_buffer,
                            &mut PanicOnExtend,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from P{index}{} to end costs {}",
                            anchors.primary(index),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function.update_primary_to_end(index, cost, true);
                            chaining_cost_function.update_additional_primary_targets(
                                index,
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.primary_to_end(index));
                }
                (Identifier::SecondaryToPrimary { index, ts_kind, .. }, Identifier::End) => {
                    let start = anchors
                        .secondary(index, ts_kind)
                        .end(k)
                        .into_specific(ts_kind);
                    if final_evaluation
                        || !chaining_cost_function.is_jump_34_to_end_exact(index, ts_kind)
                    {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.ts_34_jump_aligner.align(
                            start,
                            self.sequences.primary_end(),
                            &mut self.additional_primary_targets_buffer,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from S{}[{index}]{} to end costs {}",
                            ts_kind.digits(),
                            anchors.secondary(index, ts_kind),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function
                                .update_jump_34_to_end(index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_34_jump_targets(
                                index,
                                &mut self.additional_primary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.jump_34_to_end(index, ts_kind));
                }
                (
                    Identifier::PrimaryToPrimary {
                        index: from_index, ..
                    },
                    Identifier::PrimaryToPrimary {
                        index: to_index, ..
                    }
                    | Identifier::PrimaryToSecondary {
                        index: to_index, ..
                    },
                ) => {
                    if anchors
                        .primary(from_index)
                        .is_direct_predecessor_of(anchors.primary(to_index))
                    {
                        alignments.push(Alignment::from(vec![AlignmentType::Match]));
                        continue;
                    }

                    let start = anchors.primary(from_index).end(k);
                    let end = anchors.primary(to_index).start();
                    if final_evaluation
                        || !chaining_cost_function.is_primary_exact(from_index, to_index)
                    {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.primary_aligner.align(
                            start,
                            end,
                            &mut self.additional_primary_targets_buffer,
                            &mut PanicOnExtend,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from P{from_index}{} to P{to_index}{} (from {start} to {end}) costs {}",
                            anchors.primary(from_index),
                            anchors.primary(to_index),
                            cost
                        );

                        if end.a() - start.a() > usize::try_from(max_match_run).unwrap() {
                            assert!(
                                !cost.is_zero(),
                                "Alignment is longer than max_match_run, but has zero cost: {}",
                                alignment
                            );
                        }

                        if !final_evaluation {
                            chaining_cost_function.update_primary(from_index, to_index, cost, true);
                            chaining_cost_function.update_additional_primary_targets(
                                from_index,
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.primary(from_index, to_index));
                }
                (
                    Identifier::PrimaryToSecondary {
                        index: from_index, ..
                    },
                    Identifier::SecondaryToSecondary {
                        index: to_index,
                        ts_kind,
                        ..
                    }
                    | Identifier::SecondaryToPrimary {
                        index: to_index,
                        ts_kind,
                        ..
                    },
                ) => {
                    let start = anchors.primary(from_index).end(k);
                    let end = anchors
                        .secondary(to_index, ts_kind)
                        .start()
                        .into_specific(ts_kind);
                    if final_evaluation
                        || !chaining_cost_function.is_jump_12_exact(from_index, to_index, ts_kind)
                    {
                        self.additional_secondary_targets_buffer.clear();
                        let (cost, alignment) = self.ts_12_jump_aligner.align(
                            start,
                            end,
                            &mut self.additional_secondary_targets_buffer,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_secondary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from P{from_index}{} to S{}[{to_index}]{} costs {}",
                            anchors.primary(from_index),
                            ts_kind.digits(),
                            anchors.secondary(to_index, ts_kind),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function
                                .update_jump_12(from_index, to_index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_12_jump_targets(
                                from_index,
                                &mut self.additional_secondary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound.saturating_add(
                        &chaining_cost_function.jump_12(from_index, to_index, ts_kind),
                    );
                }
                (
                    Identifier::SecondaryToSecondary {
                        index: from_index,
                        ts_kind,
                        ..
                    },
                    Identifier::SecondaryToSecondary {
                        index: to_index,
                        ts_kind: to_ts_kind,
                        ..
                    }
                    | Identifier::SecondaryToPrimary {
                        index: to_index,
                        ts_kind: to_ts_kind,
                        ..
                    },
                ) => {
                    assert_eq!(ts_kind, to_ts_kind);
                    if anchors
                        .secondary(from_index, ts_kind)
                        .is_direct_predecessor_of(anchors.secondary(to_index, ts_kind))
                    {
                        alignments.push(Alignment::from(vec![AlignmentType::Match]));
                        continue;
                    }
                    let start = anchors
                        .secondary(from_index, ts_kind)
                        .end(k)
                        .into_specific(ts_kind);
                    let end = anchors
                        .secondary(to_index, ts_kind)
                        .start()
                        .into_specific(ts_kind);
                    if final_evaluation
                        || !chaining_cost_function.is_secondary_exact(from_index, to_index, ts_kind)
                    {
                        self.additional_secondary_targets_buffer.clear();
                        let (cost, alignment) = self.secondary_aligner.align(
                            start,
                            end,
                            &mut PanicOnExtend,
                            &mut self.additional_secondary_targets_buffer,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_secondary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from S{}[{from_index}]{} to S{}[{to_index}]{} costs {}",
                            ts_kind.digits(),
                            anchors.secondary(from_index, ts_kind),
                            ts_kind.digits(),
                            anchors.secondary(to_index, ts_kind),
                            cost
                        );
                        if !final_evaluation {
                            chaining_cost_function
                                .update_secondary(from_index, to_index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_secondary_targets(
                                from_index,
                                &mut self.additional_secondary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound.saturating_add(
                        &chaining_cost_function.secondary(from_index, to_index, ts_kind),
                    );
                }
                (
                    Identifier::SecondaryToPrimary {
                        index: from_index,
                        ts_kind,
                        ..
                    },
                    Identifier::PrimaryToPrimary {
                        index: to_index, ..
                    }
                    | Identifier::PrimaryToSecondary {
                        index: to_index, ..
                    },
                ) => {
                    let start = anchors
                        .secondary(from_index, ts_kind)
                        .end(k)
                        .into_specific(ts_kind);
                    let end = anchors.primary(to_index).start();
                    if final_evaluation
                        || !chaining_cost_function.is_jump_34_exact(from_index, to_index, ts_kind)
                    {
                        self.additional_primary_targets_buffer.clear();
                        let (cost, alignment) = self.ts_34_jump_aligner.align(
                            start,
                            end,
                            &mut self.additional_primary_targets_buffer,
                        );
                        self.total_gap_fillings +=
                            u64::try_from(self.additional_primary_targets_buffer.len()).unwrap();
                        gap_fill_alignment_count += 1;

                        trace!(
                            "Aligning from S{}[{from_index}]{} to P{to_index}{} (S{} to P{}) costs {}",
                            ts_kind.digits(),
                            anchors.secondary(from_index, ts_kind),
                            anchors.primary(to_index),
                            start,
                            end,
                            cost,
                        );
                        if !final_evaluation {
                            chaining_cost_function
                                .update_jump_34(from_index, to_index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_34_jump_targets(
                                from_index,
                                &mut self.additional_primary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                        alignments.push(iter::repeat_n(AlignmentType::Match, k).collect());
                        alignments.push(alignment);
                    }
                    current_upper_bound = current_upper_bound.saturating_add(
                        &chaining_cost_function.jump_34(from_index, to_index, ts_kind),
                    );
                }
                (Identifier::End, _)
                | (_, Identifier::Start)
                | (
                    Identifier::SecondaryToPrimary { .. } | Identifier::PrimaryToPrimary { .. },
                    Identifier::SecondaryToPrimary { .. } | Identifier::SecondaryToSecondary { .. },
                )
                | (
                    Identifier::PrimaryToSecondary { .. } | Identifier::SecondaryToSecondary { .. },
                    Identifier::PrimaryToPrimary { .. }
                    | Identifier::PrimaryToSecondary { .. }
                    | Identifier::End,
                )
                | (Identifier::StartToPrimary { .. } | Identifier::StartToSecondary { .. }, _)
                | (_, Identifier::StartToPrimary { .. } | Identifier::StartToSecondary { .. }) => {
                    unreachable!()
                }
            }
        }

        self.gap_fill_alignments_per_chain
            .push(gap_fill_alignment_count);

        (current_upper_bound, alignments)
    }

    pub fn total_gaps(&self) -> u64 {
        self.total_gaps
    }

    pub fn total_gap_fillings(&self) -> u64 {
        self.total_gap_fillings
    }

    pub fn total_redundant_gap_fillings(&self) -> u64 {
        self.total_redundant_gap_fillings
    }

    pub fn gap_fill_alignments_per_chain(&self) -> &[u32] {
        &self.gap_fill_alignments_per_chain
    }
}
