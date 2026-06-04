use binary_heap_plus::BinaryHeap;
use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{
        alphabet::Alphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
};
use generic_a_star::{
    AStar, AStarResult,
    closed_lists::AStarClosedList,
    comparator::AStarNodeComparator,
    cost::AStarCost,
    open_lists::{AStarOpenList, linear_heap::LinearHeap},
};
use indicatif::ProgressBar;
use itertools::Itertools;
use lib_tsalign::a_star_aligner::{
    alignment_geometry::AlignmentRange,
    alignment_result::AlignmentResult,
    template_switch_distance::{EqualCostRange, TemplateSwitchDirection},
};
use log::{debug, info, trace};
use rustc_hash::FxHashMapSeed;
use std::{
    fmt::Write,
    iter,
    time::{Duration, Instant},
};

use crate::{
    alignment::{AlignmentType, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::{
        chainer::{Context, Identifier, Node, closed_list::ChainerClosedList},
        evaluation::ChainEvaluator,
        performance_parameters::{
            AlignmentPerformanceParameters, ChainingClosedList, ChainingOpenList,
        },
    },
    chaining_cost_function::ChainingCostFunction,
    costs::AlignmentCosts,
};

mod chainer;
mod evaluation;
pub mod performance_parameters;

pub fn align<AlphabetType: Alphabet, Cost: AStarCost>(
    sequences: &AlignmentSequences,
    performance_parameters: &AlignmentPerformanceParameters<Cost>,
    alignment_costs: &AlignmentCosts<Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    max_match_run: u32,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    match performance_parameters.open_list {
        ChainingOpenList::StdHeap => {
            choose_closed_list::<AlphabetType, _, BinaryHeap<Node<Cost>, AStarNodeComparator>>(
                sequences,
                performance_parameters,
                alignment_costs,
                rc_fn,
                max_match_run,
                anchors,
                chaining_cost_function,
            )
        }
        ChainingOpenList::LinearHeap => choose_closed_list::<AlphabetType, _, LinearHeap<_>>(
            sequences,
            performance_parameters,
            alignment_costs,
            rc_fn,
            max_match_run,
            anchors,
            chaining_cost_function,
        ),
    }
}

pub fn choose_closed_list<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    OpenList: AStarOpenList<Node<Cost>>,
>(
    sequences: &AlignmentSequences,
    performance_parameters: &AlignmentPerformanceParameters<Cost>,
    alignment_costs: &AlignmentCosts<Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    max_match_run: u32,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    match performance_parameters.closed_list {
        ChainingClosedList::FxHashMap => {
            actually_align::<AlphabetType, _, FxHashMapSeed<_, _>, OpenList>(
                sequences,
                performance_parameters,
                alignment_costs,
                rc_fn,
                max_match_run,
                anchors,
                chaining_cost_function,
            )
        }
        ChainingClosedList::Special => {
            actually_align::<AlphabetType, _, ChainerClosedList<_>, OpenList>(
                sequences,
                performance_parameters,
                alignment_costs,
                rc_fn,
                max_match_run,
                anchors,
                chaining_cost_function,
            )
        }
    }
}

fn actually_align<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    ClosedList: AStarClosedList<Node<Cost>>,
    OpenList: AStarOpenList<Node<Cost>>,
>(
    sequences: &AlignmentSequences,
    performance_parameters: &AlignmentPerformanceParameters<Cost>,
    alignment_costs: &AlignmentCosts<Cost>,
    rc_fn: &dyn Fn(u8) -> u8,
    max_match_run: u32,
    anchors: &Anchors,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    info!("Aligning...");
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.enable_steady_tick(Duration::from_millis(200));

    let start_time = Instant::now();
    let mut chaining_duration = Duration::default();
    let mut evaluation_duration = Duration::default();

    let k = usize::try_from(max_match_run + 1).unwrap();
    let context = Context::new(
        anchors,
        chaining_cost_function,
        &alignment_costs.ts_limits,
        k,
        performance_parameters.max_successors,
    );
    let mut astar = AStar::<_, ClosedList, OpenList>::new(context);

    let mut chain_evaluator = ChainEvaluator::new(sequences, alignment_costs, rc_fn, max_match_run);

    let mut chaining_execution_count = 0;
    let mut current_lower_bound = Cost::zero();
    let mut current_upper_bound = Cost::max_value();
    let mut total_chaining_opened_nodes = 0;
    let mut total_chaining_suboptimal_opened_nodes = 0;
    let mut total_chaining_closed_nodes = 0;

    let (chain, result) = loop {
        progress_bar.inc(1);
        progress_bar.set_message(format!(
            "Computing chain number {} (cost {}->{})",
            chaining_execution_count + 1,
            current_lower_bound,
            current_upper_bound,
        ));
        let chaining_start_time = Instant::now();

        astar.reset();
        astar.initialise();
        let result = astar.search();
        total_chaining_opened_nodes += astar.performance_counters().opened_nodes;
        total_chaining_suboptimal_opened_nodes +=
            astar.performance_counters().suboptimal_opened_nodes;
        total_chaining_closed_nodes += astar.performance_counters().closed_nodes;

        chaining_execution_count += 1;
        let chain = match result {
            AStarResult::FoundTarget { cost, .. } => {
                trace!("Found chain with cost {cost}");
                current_lower_bound = cost;

                let mut chain = astar.reconstruct_path();
                chain.push(Identifier::End);
                trace!("Chain (len: {}):\n{}", chain.len(), {
                    let mut s = String::new();
                    let mut once = true;
                    for identifier in &chain {
                        if once {
                            once = false;
                        } else {
                            writeln!(s).unwrap();
                        }
                        match identifier {
                            Identifier::Start => write!(s, "start").unwrap(),
                            Identifier::StartToPrimary { offset } => {
                                write!(s, "start-to-primary-{offset}").unwrap()
                            }
                            Identifier::StartToSecondary { ts_kind, offset } => {
                                write!(s, "start-to-secondary{}-{offset}", ts_kind.digits())
                                    .unwrap()
                            }
                            Identifier::PrimaryToPrimary { index, offset } => {
                                write!(s, "P{}-to-primary-{offset}", anchors.primary(*index))
                                    .unwrap()
                            }
                            Identifier::PrimaryToSecondary {
                                index,
                                ts_kind,
                                offset,
                            } => write!(
                                s,
                                "P{}-to-secondary{}-{offset}",
                                anchors.primary(*index),
                                ts_kind.digits()
                            )
                            .unwrap(),
                            Identifier::SecondaryToSecondary {
                                index,
                                ts_kind,
                                first_secondary_index,
                                offset,
                            } => write!(
                                s,
                                "S{}{}->{}-to-secondary-{offset}",
                                ts_kind.digits(),
                                anchors.secondary(*first_secondary_index, *ts_kind),
                                anchors.secondary(*index, *ts_kind),
                            )
                            .unwrap(),
                            Identifier::SecondaryToPrimary {
                                index,
                                ts_kind,
                                first_secondary_index,
                                offset,
                            } => write!(
                                s,
                                "S{}{}->{}-to-primary-{offset}",
                                ts_kind.digits(),
                                anchors.secondary(*first_secondary_index, *ts_kind),
                                anchors.secondary(*index, *ts_kind),
                            )
                            .unwrap(),

                            Identifier::End => write!(s, "end").unwrap(),
                        }
                    }
                    s
                });
                chain
            }
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Memory limit is None"),
            AStarResult::NoTarget => panic!("No chain found"),
        };

        let chaining_end_time = Instant::now();
        chaining_duration += chaining_end_time - chaining_start_time;

        let evaluation_start_time = Instant::now();

        let (evaluated_cost, _) = chain_evaluator.evaluate_chain(
            anchors,
            &chain,
            max_match_run,
            astar.context_mut().chaining_cost_function,
            false,
        );
        let cost_increased = evaluated_cost > current_lower_bound;
        current_upper_bound = evaluated_cost;

        let evaluation_end_time = Instant::now();
        evaluation_duration += evaluation_end_time - evaluation_start_time;

        if !cost_increased {
            break (chain, result);
        }
    };

    progress_bar.finish_and_clear();
    let total_gap_fill_alignments: u32 =
        chain_evaluator.gap_fill_alignments_per_chain().iter().sum();
    let chains_with_at_most_exp2x_gap_alignments_amount: Vec<_> = iter::once(
        chain_evaluator
            .gap_fill_alignments_per_chain()
            .iter()
            .filter(|a| **a == 0)
            .count(),
    )
    .chain((0..u32::BITS).map(|e| {
        chain_evaluator
            .gap_fill_alignments_per_chain()
            .iter()
            .filter(|a| **a <= 1 << e)
            .count()
    }))
    .take_while_inclusive(|amount| *amount < chaining_execution_count)
    .collect();

    debug!("Computed {chaining_execution_count} chains");
    debug!("Chaining took {:.1}s", chaining_duration.as_secs_f64());
    debug!("Evaluation took {:.1}s", evaluation_duration.as_secs_f64());
    debug!("Chaining opened nodes: {total_chaining_opened_nodes}");
    debug!(
        "Chaining suboptimal openend nodes: {} ({:.0}% of opened nodes)",
        total_chaining_suboptimal_opened_nodes,
        total_chaining_suboptimal_opened_nodes as f64 / total_chaining_opened_nodes as f64 * 1e2,
    );
    debug!("Chaining closed nodes: {total_chaining_closed_nodes}");
    debug!("Total chain gaps: {}", chain_evaluator.total_gaps());
    debug!(
        "Total chain gap alignments: {} ({:.0}% of total gaps)",
        total_gap_fill_alignments,
        total_gap_fill_alignments as f64 / chain_evaluator.total_gaps() as f64 * 1e2
    );
    debug!(
        "Fraction of chains with [0, 0] gap alignments: {:.0}%",
        *chains_with_at_most_exp2x_gap_alignments_amount
            .first()
            .unwrap() as f64
            / chaining_execution_count as f64
            * 1e2,
    );
    for (e, amount) in chains_with_at_most_exp2x_gap_alignments_amount
        .windows(2)
        .enumerate()
    {
        debug!(
            "Fraction of chains with [{}, {}] gap alignments: {:.0}%",
            if e > 0 { (1u32 << (e - 1)) + 1 } else { 1 },
            1 << e,
            (amount[1] - amount[0]) as f64 / chaining_execution_count as f64 * 1e2,
        );
    }
    debug!(
        "Total chain gap fillings: {} ({:.2}x total gaps, {:.0}% redundant)",
        chain_evaluator.total_gap_fillings(),
        chain_evaluator.total_gap_fillings() as f64 / chain_evaluator.total_gaps() as f64,
        chain_evaluator.total_redundant_gap_fillings() as f64
            / chain_evaluator.total_gap_fillings() as f64
            * 1e2,
    );

    let mut tsalign_alignment =
        lib_tsalign::a_star_aligner::alignment_result::alignment::Alignment::new();
    let mut is_primary = true;
    let mut anti_primary_gap = 0isize;

    debug!("Evaluating final chain");
    let (evaluated_cost, alignments) = chain_evaluator.evaluate_chain(
        anchors,
        &chain,
        max_match_run,
        astar.context_mut().chaining_cost_function,
        true,
    );
    assert_eq!(evaluated_cost, result.cost());

    for alignment in alignments {
        use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType as TsAlignAlignmentType;

        for (multiplicity, alignment_type) in alignment.alignment {
            match alignment_type {
                AlignmentType::Match => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryMatch
                    } else {
                        anti_primary_gap -= multiplicity as isize;
                        TsAlignAlignmentType::SecondaryMatch
                    },
                ),
                AlignmentType::Substitution => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimarySubstitution
                    } else {
                        anti_primary_gap -= multiplicity as isize;
                        TsAlignAlignmentType::SecondarySubstitution
                    },
                ),
                AlignmentType::GapA => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryInsertion
                    } else {
                        TsAlignAlignmentType::SecondaryInsertion
                    },
                ),
                AlignmentType::GapB => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        TsAlignAlignmentType::PrimaryDeletion
                    } else {
                        anti_primary_gap -= multiplicity as isize;
                        TsAlignAlignmentType::SecondaryDeletion
                    },
                ),
                AlignmentType::TsStart { jump, ts_kind } => {
                    assert!(is_primary);
                    assert_eq!(multiplicity, 1);
                    is_primary = false;
                    anti_primary_gap = jump;

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchEntrance {
                            first_offset: jump,
                            equal_cost_range: EqualCostRange::new_invalid(),
                            descendant: ts_kind.descendant.into_tsalign_primary(),
                            ancestor: ts_kind.ancestor.into_tsalign_secondary(),
                            direction: TemplateSwitchDirection::Reverse,
                        },
                    );
                }
                AlignmentType::TsEnd { jump } => {
                    assert!(!is_primary);
                    assert_eq!(multiplicity, 1);
                    is_primary = true;

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchExit {
                            anti_descendant_gap: anti_primary_gap + jump,
                        },
                    )
                }
            }
        }
    }

    let end_time = Instant::now();
    let duration_seconds = (end_time - start_time).as_secs_f64();

    let alignment_range = AlignmentRange::new(
        sequences.primary_start().primary_ordinate_a().unwrap(),
        sequences.primary_start().primary_ordinate_b().unwrap(),
        sequences.primary_end().primary_ordinate_a().unwrap(),
        sequences.primary_end().primary_ordinate_b().unwrap(),
    );

    AlignmentResult::new_with_target::<AlphabetType, _>(
        tsalign_alignment.into_inner(),
        VectorGenome::from_slice_u8(sequences.seq1())
            .unwrap()
            .as_genome_subsequence(),
        VectorGenome::from_slice_u8(sequences.seq2())
            .unwrap()
            .as_genome_subsequence(),
        sequences.seq1_name(),
        sequences.seq2_name(),
        alignment_range,
        result.without_node_identifier(),
        duration_seconds,
        0,
        0,
        0,
        sequences.seq1().len(),
        sequences.seq2().len(),
    )
}
