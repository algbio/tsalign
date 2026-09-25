use generic_a_star::cost::AStarCost;
use log::trace;
use num_traits::Zero;

use crate::{
    alignment::{
        Alignment, AlignmentType,
        coordinates::{AnySecondaryAlignmentCoordinates, PrimaryAlignmentCoordinates},
        sequences::AlignmentSequences,
    },
    alignment_history::history_vec::AlignmentHistory,
    anchors::Anchors,
    chain_align::{chainer::Identifier, evaluation::ChainEvaluator},
    chaining_cost_function::ChainingCostFunction,
    costs::AlignmentCosts,
    exact_chaining,
    inexact_chaining::{
        gap_affine::GapAffineAligner, ts_12_jump::Ts12JumpAligner, ts_34_jump::Ts34JumpAligner,
    },
    panic_on_extend::PanicOnExtend,
};

pub struct InexactChainEvaluator<
    'sequences,
    'alignment_costs,
    'rc_fn,
    Cost: AStarCost,
    AlignmentHistoryVec,
> {
    sequences: &'sequences AlignmentSequences,
    primary_aligner:
        GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>,
    secondary_aligner:
        GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>,
    ts_12_jump_aligner:
        Ts12JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>,
    ts_34_jump_aligner:
        Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>,
    primary_anchor_aligner:
        exact_chaining::gap_affine::GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,
    secondary_anchor_aligner:
        exact_chaining::gap_affine::GapAffineAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>,

    additional_primary_targets_buffer: Vec<(PrimaryAlignmentCoordinates, Cost)>,
    additional_secondary_targets_buffer: Vec<(AnySecondaryAlignmentCoordinates, Cost)>,

    total_gaps: u64,
    total_gap_fillings: u64,
    total_redundant_gap_fillings: u64,
    gap_fill_alignments_per_chain: Vec<u32>,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory>
    InexactChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Self {
        Self {
            sequences,
            primary_aligner: GapAffineAligner::new(
                sequences,
                &alignment_costs.primary_costs,
                rc_fn,
                anchor_k,
                max_anchor_mutations,
            ),
            secondary_aligner: GapAffineAligner::new(
                sequences,
                &alignment_costs.secondary_costs,
                rc_fn,
                anchor_k,
                max_anchor_mutations,
            ),
            ts_12_jump_aligner: Ts12JumpAligner::new(
                sequences,
                alignment_costs,
                rc_fn,
                anchor_k,
                max_anchor_mutations,
            ),
            ts_34_jump_aligner: Ts34JumpAligner::new(
                sequences,
                alignment_costs,
                rc_fn,
                anchor_k,
                max_anchor_mutations,
            ),
            primary_anchor_aligner: exact_chaining::gap_affine::GapAffineAligner::new(
                sequences,
                &alignment_costs.primary_costs,
                rc_fn,
                (anchor_k - 1).into(),
            ),
            secondary_anchor_aligner: exact_chaining::gap_affine::GapAffineAligner::new(
                sequences,
                &alignment_costs.secondary_costs,
                rc_fn,
                (anchor_k - 1).into(),
            ),

            additional_primary_targets_buffer: Default::default(),
            additional_secondary_targets_buffer: Default::default(),

            total_gaps: 0,
            total_gap_fillings: 0,
            total_redundant_gap_fillings: 0,
            gap_fill_alignments_per_chain: Default::default(),
        }
    }
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost, AlignmentHistoryVec: AlignmentHistory>
    ChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost>
    for InexactChainEvaluator<'sequences, 'alignment_costs, 'rc_fn, Cost, AlignmentHistoryVec>
{
    fn evaluate_chain(
        &mut self,
        anchors: &Anchors<Cost>,
        chain: &[Identifier],
        anchor_k: u32,
        max_anchor_mutations: u8,
        chaining_cost_function: &mut ChainingCostFunction<Cost>,
        final_evaluation: bool,
    ) -> (Cost, Vec<Alignment>) {
        debug_assert_eq!(max_anchor_mutations, 0);

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
                        if final_evaluation {
                            alignments.push(alignment);
                        } else {
                            chaining_cost_function.update_start_to_end(cost, true);
                            chaining_cost_function.update_additional_primary_targets_from_start(
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
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
                        if final_evaluation {
                            alignments.push(alignment);
                        } else {
                            chaining_cost_function.update_primary_from_start(index, cost, true);
                            chaining_cost_function.update_additional_primary_targets_from_start(
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
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
                        if final_evaluation {
                            alignments.push(alignment);
                        } else {
                            chaining_cost_function
                                .update_jump_12_from_start(index, ts_kind, cost, true);
                            chaining_cost_function.update_additional_12_jump_targets_from_start(
                                &mut self.additional_secondary_targets_buffer,
                                ts_kind,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.jump_12_from_start(index, ts_kind));
                }
                (Identifier::PrimaryToPrimary { index, .. }, Identifier::End) => {
                    let start = anchors.primary(index).end();
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
                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.primary_anchor_aligner.align_anchor(
                                    anchors.primary(index).start(),
                                    anchors.primary(index).end(),
                                );

                            assert_eq!(anchor_alignment_cost, anchors.primary(index).cost());
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
                            chaining_cost_function.update_primary_to_end(index, cost, true);
                            chaining_cost_function.update_additional_primary_targets(
                                index,
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.primary_to_end(index))
                        .saturating_add(&anchors.primary(index).cost());
                }
                (Identifier::SecondaryToPrimary { index, ts_kind, .. }, Identifier::End) => {
                    let start = anchors
                        .secondary(index, ts_kind)
                        .end()
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
                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.secondary_anchor_aligner.align_anchor(
                                    anchors
                                        .secondary(index, ts_kind)
                                        .start()
                                        .into_specific(ts_kind),
                                    anchors
                                        .secondary(index, ts_kind)
                                        .end()
                                        .into_specific(ts_kind),
                                );

                            assert_eq!(
                                anchor_alignment_cost,
                                anchors.secondary(index, ts_kind).cost()
                            );
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
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
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.jump_34_to_end(index, ts_kind))
                        .saturating_add(&anchors.secondary(index, ts_kind).cost());
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
                        .is_direct_free_predecessor_of(anchors.primary(to_index))
                    {
                        if final_evaluation {
                            alignments.push(Alignment::from(vec![AlignmentType::Match]));
                        }
                        continue;
                    }

                    let start = anchors.primary(from_index).end();
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

                        if end.a() - start.a() > usize::try_from(anchor_k - 1).unwrap() {
                            assert!(
                                !cost.is_zero(),
                                "Alignment is longer than max_match_run, but has zero cost: {}",
                                alignment
                            );
                        }

                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.primary_anchor_aligner.align_anchor(
                                    anchors.primary(from_index).start(),
                                    anchors.primary(from_index).end(),
                                );

                            assert_eq!(anchor_alignment_cost, anchors.primary(from_index).cost());
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
                            chaining_cost_function.update_primary(from_index, to_index, cost, true);
                            chaining_cost_function.update_additional_primary_targets(
                                from_index,
                                &mut self.additional_primary_targets_buffer,
                                anchors,
                                &mut self.total_redundant_gap_fillings,
                            );
                        }
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(&chaining_cost_function.primary(from_index, to_index))
                        .saturating_add(&anchors.primary(from_index).cost());
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
                    let start = anchors.primary(from_index).end();
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
                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.primary_anchor_aligner.align_anchor(
                                    anchors.primary(from_index).start(),
                                    anchors.primary(from_index).end(),
                                );

                            assert_eq!(anchor_alignment_cost, anchors.primary(from_index).cost());
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
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
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(
                            &chaining_cost_function.jump_12(from_index, to_index, ts_kind),
                        )
                        .saturating_add(&anchors.primary(from_index).cost());
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
                        .is_direct_free_predecessor_of(anchors.secondary(to_index, ts_kind))
                    {
                        if final_evaluation {
                            alignments.push(Alignment::from(vec![AlignmentType::Match]));
                        }
                        continue;
                    }
                    let start = anchors
                        .secondary(from_index, ts_kind)
                        .end()
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
                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.secondary_anchor_aligner.align_anchor(
                                    anchors
                                        .secondary(from_index, ts_kind)
                                        .start()
                                        .into_specific(ts_kind),
                                    anchors
                                        .secondary(from_index, ts_kind)
                                        .end()
                                        .into_specific(ts_kind),
                                );

                            assert_eq!(
                                anchor_alignment_cost,
                                anchors.secondary(from_index, ts_kind).cost()
                            );
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
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
                    }
                    current_upper_bound = current_upper_bound
                        .saturating_add(
                            &chaining_cost_function.secondary(from_index, to_index, ts_kind),
                        )
                        .saturating_add(&anchors.secondary(from_index, ts_kind).cost());
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
                        .end()
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
                        if final_evaluation {
                            let (anchor_alignment_cost, anchor_alignment) =
                                self.secondary_anchor_aligner.align_anchor(
                                    anchors
                                        .secondary(from_index, ts_kind)
                                        .start()
                                        .into_specific(ts_kind),
                                    anchors
                                        .secondary(from_index, ts_kind)
                                        .end()
                                        .into_specific(ts_kind),
                                );

                            assert_eq!(
                                anchor_alignment_cost,
                                anchors.secondary(from_index, ts_kind).cost()
                            );
                            alignments.push(anchor_alignment);
                            alignments.push(alignment);
                        } else {
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
                    }
                    current_upper_bound = current_upper_bound.saturating_add(
                        &chaining_cost_function
                            .jump_34(from_index, to_index, ts_kind)
                            .saturating_add(&anchors.secondary(from_index, ts_kind).cost()),
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

    fn total_gaps(&self) -> u64 {
        self.total_gaps
    }

    fn total_gap_fillings(&self) -> u64 {
        self.total_gap_fillings
    }

    fn total_redundant_gap_fillings(&self) -> u64 {
        self.total_redundant_gap_fillings
    }

    fn gap_fill_alignments_per_chain(&self) -> &[u32] {
        &self.gap_fill_alignments_per_chain
    }
}
