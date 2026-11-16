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
    alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    template_switch_distance::TemplateSwitchDirection,
};
use log::{debug, info, trace};
use rustc_hash::FxHashMapSeed;
use std::{
    fmt::Write,
    iter,
    time::{Duration, Instant},
};

use crate::{
    alignment::{
        Alignment, AlignmentType,
        sequences::AlignmentSequences,
        ts_kind::{TsAncestor, TsDescendant},
    },
    anchors::Anchors,
    chain_align::{
        chainer::{Context, Identifier, Node, closed_list::ChainerClosedList},
        evaluation::build_chain_evaluator,
        performance_parameters::{
            AlignmentPerformanceParameters, ChainingClosedList, ChainingOpenList,
        },
    },
    chaining_cost_function::ChainingCostFunction,
    config::ChainingLowerBoundConfig,
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
    config: &ChainingLowerBoundConfig,
    anchors: &Anchors<Cost>,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    match performance_parameters.open_list {
        ChainingOpenList::StdHeap => {
            choose_closed_list::<AlphabetType, _, BinaryHeap<Node<Cost>, AStarNodeComparator>>(
                sequences,
                performance_parameters,
                alignment_costs,
                rc_fn,
                config,
                anchors,
                chaining_cost_function,
            )
        }
        ChainingOpenList::LinearHeap => choose_closed_list::<AlphabetType, _, LinearHeap<_>>(
            sequences,
            performance_parameters,
            alignment_costs,
            rc_fn,
            config,
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
    config: &ChainingLowerBoundConfig,
    anchors: &Anchors<Cost>,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    match performance_parameters.closed_list {
        ChainingClosedList::FxHashMap => {
            actually_align::<AlphabetType, _, FxHashMapSeed<_, _>, OpenList>(
                sequences,
                performance_parameters,
                alignment_costs,
                rc_fn,
                config,
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
                config,
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
    config: &ChainingLowerBoundConfig,
    anchors: &Anchors<Cost>,
    chaining_cost_function: &mut ChainingCostFunction<Cost>,
) -> AlignmentResult<lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType, Cost> {
    info!("Aligning...");
    let progress_bar = ProgressBar::new_spinner();
    progress_bar.enable_steady_tick(Duration::from_millis(200));

    let start_time = Instant::now();
    let mut chaining_duration = Duration::default();
    let mut evaluation_duration = Duration::default();

    let context = Context::new(
        anchors,
        chaining_cost_function,
        &alignment_costs.ts_limits,
        performance_parameters.max_successors,
    );
    let mut astar = AStar::<_, ClosedList, OpenList>::new(context);

    let mut chain_evaluator = build_chain_evaluator(sequences, alignment_costs, rc_fn, config);

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
            config,
            astar.context_mut().chaining_cost_function,
            false,
        );
        let cost_increased = evaluated_cost > current_lower_bound;
        // PROTOTYPE: per-iteration statistics.
        debug!(
            "CHAIN_ITER {} lower_bound={} evaluated={} opened={} chain_len={}",
            chaining_execution_count,
            current_lower_bound,
            evaluated_cost,
            astar.performance_counters().opened_nodes,
            chain.len()
        );
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
    // The current primary coordinates, needed to convert the jumps of intra template switches.
    let mut primary_a = sequences.primary_start().a();
    let mut primary_b = sequences.primary_start().b();
    // The descendant and anti-descendant indices at the entrance of the current template switch,
    // and the amount of descendant characters consumed by it so far.
    let mut ts_descendant_is_seq1 = false;
    let mut entrance_descendant = 0usize;
    let mut entrance_anti_descendant = 0usize;
    let mut descendant_consumed = 0usize;

    debug!("Evaluating final chain");
    let (chain_cost, alignments) = chain_evaluator.evaluate_chain(
        anchors,
        &chain,
        config,
        astar.context_mut().chaining_cost_function,
        true,
    );
    // Both are chain costs: the sum of the costs of the gaps and anchors of the chain.
    assert_eq!(chain_cost, result.cost());

    let cost = alignment_costs.alignment_cost(alignments.iter().flat_map(Alignment::iter_unpacked));
    debug!("Chain cost: {chain_cost}");
    debug!("Cost of the final alignment: {cost}");
    check_final_alignment_cost(cost, chain_cost, alignment_costs, config);

    for alignment in alignments {
        use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType as TsAlignAlignmentType;

        for (multiplicity, alignment_type) in alignment.alignment {
            match alignment_type {
                AlignmentType::Match => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        primary_a += multiplicity;
                        primary_b += multiplicity;
                        TsAlignAlignmentType::PrimaryMatch
                    } else {
                        anti_primary_gap -= multiplicity as isize;
                        descendant_consumed += multiplicity;
                        TsAlignAlignmentType::SecondaryMatch
                    },
                ),
                AlignmentType::Substitution => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        primary_a += multiplicity;
                        primary_b += multiplicity;
                        TsAlignAlignmentType::PrimarySubstitution
                    } else {
                        anti_primary_gap -= multiplicity as isize;
                        descendant_consumed += multiplicity;
                        TsAlignAlignmentType::SecondarySubstitution
                    },
                ),
                AlignmentType::GapA => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        primary_b += multiplicity;
                        TsAlignAlignmentType::PrimaryInsertion
                    } else {
                        descendant_consumed += multiplicity;
                        TsAlignAlignmentType::SecondaryInsertion
                    },
                ),
                AlignmentType::GapB => tsalign_alignment.push_n(
                    multiplicity,
                    if is_primary {
                        primary_a += multiplicity;
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

                    ts_descendant_is_seq1 = ts_kind.descendant == TsDescendant::Seq1;
                    (entrance_descendant, entrance_anti_descendant) = if ts_descendant_is_seq1 {
                        (primary_a, primary_b)
                    } else {
                        (primary_b, primary_a)
                    };
                    descendant_consumed = 0;

                    // The jumps are relative to the ancestor, and the 34-jump ends in the anti-descendant.
                    // For inter template switches, the ancestor is the anti-descendant,
                    // so the sum of the jumps and the ancestor movement is the anti-descendant gap.
                    // For intra template switches, the ancestor is the descendant,
                    // so the 12-jump needs to be shifted from the descendant to the anti-descendant.
                    let is_intra = matches!(
                        (ts_kind.ancestor, ts_kind.descendant),
                        (TsAncestor::Seq1, TsDescendant::Seq1)
                            | (TsAncestor::Seq2, TsDescendant::Seq2)
                    );
                    anti_primary_gap = if is_intra {
                        jump + entrance_descendant as isize - entrance_anti_descendant as isize
                    } else {
                        jump
                    };

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchEntrance {
                            first_offset: jump,
                            uncertainty_range: None,
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
                    let anti_descendant_gap = anti_primary_gap + jump;

                    // Continue the primary alignment after the template switch.
                    let descendant = entrance_descendant + descendant_consumed;
                    let anti_descendant =
                        usize::try_from(entrance_anti_descendant as isize + anti_descendant_gap)
                            .unwrap();
                    (primary_a, primary_b) = if ts_descendant_is_seq1 {
                        (descendant, anti_descendant)
                    } else {
                        (anti_descendant, descendant)
                    };

                    tsalign_alignment.push_n(
                        multiplicity,
                        TsAlignAlignmentType::TemplateSwitchExit {
                            anti_descendant_gap,
                        },
                    )
                }
            }
        }
    }

    let end_time = Instant::now();
    let duration_seconds = (end_time - start_time).as_secs_f64();

    let alignment_range = AlignmentRange::new(
        sequences.primary_start().a(),
        sequences.primary_start().b(),
        sequences.primary_end().a(),
        sequences.primary_end().b(),
    );

    let mut alignment_result = AlignmentResult::new_with_target::<AlphabetType, _>(
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
        AStarResult::FoundTarget {
            identifier: (),
            cost,
        },
        duration_seconds,
        0,
        0,
        0,
        sequences.seq1().len(),
        sequences.seq2().len(),
    );
    alignment_result.statistics_mut().chain_cost = Some(chain_cost.as_f64().try_into().unwrap());
    alignment_result
}

/// Checks the cost of the final alignment against the cost of the chain it was built from.
///
/// tschainalign finds alignments that are optimal with regard to the anchors, not globally optimal alignments.
/// The chain cost is the sum of the costs of the gaps and anchors of the chain, which are aligned independently.
/// For inexact anchors, the cost of the final alignment may be lower than the chain cost
/// (see [issue #221](https://github.com/algbio/tsalign/issues/221)):
///
/// * If a gap (or anchor) alignment ends with a gap and the following anchor (or gap) alignment starts with a gap
///   of the same kind, then the final alignment contains a single gap, which pays the gap open cost only once.
///   This can happen between a gap and an anchor, an anchor and a gap, two adjacent anchors (with an empty gap),
///   and between the primary (secondary) parts of 12- and 34-jumps and the adjacent primary (secondary) anchors.
///   The gap state is reset at the start and the end of a template switch, so it is not carried over the jumps.
/// * The final evaluation re-aligns each anchor without the limit on the number of mutations,
///   which may be cheaper than the anchor cost, which is limited to `max_anchor_mutations` mutations.
///
/// Merging gaps never increases the cost if the gap extend cost is at most the gap open cost,
/// so then the cost of the final alignment is at most the chain cost.
/// Exact anchors consist of matches only, so they neither merge gaps nor have a cheaper re-alignment,
/// and the costs are equal.
fn check_final_alignment_cost<Cost: AStarCost>(
    cost: Cost,
    chain_cost: Cost,
    alignment_costs: &AlignmentCosts<Cost>,
    config: &ChainingLowerBoundConfig,
) {
    if alignment_costs.gap_merging_never_increases_cost() {
        assert!(
            cost <= chain_cost,
            "The cost {cost} of the final alignment is larger than the chain cost {chain_cost}"
        );
    }
    if config.max_anchor_mutations == 0 {
        assert_eq!(
            cost, chain_cost,
            "With exact anchors, the cost of the final alignment must equal the chain cost"
        );
    }
}
