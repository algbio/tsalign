use std::{
    collections::BTreeMap,
    fmt::{Display, Formatter, Result, Write},
};

use a_star_sequences::SequencePair;
use alignment::Alignment;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStarResult, cost::AStarCost};
use log::{debug, trace, warn};
use noisy_float::types::{R64, r64};
use num_traits::{Float, Zero};

use crate::{alignment_geometry::AlignmentRange, config::TemplateSwitchConfig};

pub mod a_star_sequences;
pub mod alignment;
pub mod alignment_type;

pub use alignment_type::{
    AlignmentType,
    tsm_uncertainty_range::{TSMUncertaintyRange, TSMUncertaintyRangeExtensionMode},
};

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
pub enum AlignmentResult<Cost> {
    WithTarget {
        #[cfg_attr(feature = "serde", serde(flatten))]
        alignment: Alignment,

        #[cfg_attr(feature = "serde", serde(flatten))]
        statistics: AlignmentStatistics<Cost>,
    },

    WithoutTarget {
        #[cfg_attr(feature = "serde", serde(flatten))]
        statistics: AlignmentStatistics<Cost>,
    },
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[must_use]
pub struct AlignmentStatistics<Cost> {
    pub result: AStarResult<(), Cost>,

    // Input data used by show.
    pub sequences: SequencePair,
    pub reference_offset: usize,
    pub query_offset: usize,
    pub reference_limit: usize,
    pub query_limit: usize,

    pub cost: R64,
    pub cost_per_base: R64,
    pub duration_seconds: R64,
    pub opened_nodes: R64,
    pub closed_nodes: R64,
    pub suboptimal_opened_nodes: R64,
    pub suboptimal_opened_nodes_ratio: R64,
    pub template_switch_amount: R64,

    /// The cost of the chain that the alignment was built from, if it was computed by tschainalign.
    ///
    /// The chain cost is the sum of the costs of the gaps and anchors of the chain, which are aligned separately.
    /// It may be larger than the cost of the alignment, see [issue #221](https://github.com/algbio/tsalign/issues/221).
    /// This is for debugging only, and is ignored by the piecewise operations on statistics.
    #[cfg_attr(
        feature = "serde",
        serde(default, skip_serializing_if = "Option::is_none")
    )]
    pub chain_cost: Option<R64>,

    /// The iteration limit of tschainalign, if tschainalign stopped because it reached it.
    ///
    /// The iteration limit is the maximum number of chaining rounds.
    /// If it is set, then there is no alignment, `result` is [`AStarResult::NoTarget`],
    /// and `cost` is the lower bound of the chain cost reached before stopping.
    #[cfg_attr(
        feature = "serde",
        serde(default, skip_serializing_if = "Option::is_none")
    )]
    pub exceeded_iteration_limit: Option<u64>,

    /// Statistics about the anchors and the chaining of tschainalign, if the alignment was computed by tschainalign.
    ///
    /// This is ignored by the piecewise operations on statistics.
    #[cfg_attr(
        feature = "serde",
        serde(default, skip_serializing_if = "Option::is_none")
    )]
    pub chain_statistics: Option<ChainStatistics>,

    /// Runtime in seconds.
    ///
    /// To be filled by some other tool, not collected by tsalign.
    #[cfg_attr(feature = "serde", serde(default))]
    pub runtime: R64,

    /// Memory in bytes.
    ///
    /// To be filled by some other tool, not collected by tsalign.
    #[cfg_attr(feature = "serde", serde(default))]
    pub memory: R64,
}

/// Statistics about the anchors and the chaining of tschainalign.
///
/// The anchor counts are summed over all template switch kinds for secondary anchors.
/// Pruning statistics are keyed by the name of the pruning rule, e.g. `primary-prefix-domination`
/// (see `lib_ts_chainalign::chain_align::performance_parameters::InexactPruningRule`).
/// The exact pruning mode reports its removed anchors as `primary-exact-domination`.
#[derive(Debug, Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(default))]
pub struct ChainStatistics {
    /// The anchor pruning mode (`none`, `exact` or `inexact`).
    pub anchor_pruning: String,
    /// The enabled inexact pruning rules as a comma-separated list, or `none`.
    /// Empty if the pruning mode is not `inexact`.
    pub inexact_pruning_rules: String,
    /// The number of primary anchors before pruning.
    pub primary_anchors: u64,
    /// The number of primary anchors after pruning.
    pub primary_anchors_after_pruning: u64,
    /// The number of secondary anchors before pruning.
    pub secondary_anchors: u64,
    /// The number of secondary anchors after pruning.
    pub secondary_anchors_after_pruning: u64,
    /// The time for generating the anchors, in seconds.
    pub anchor_seconds: R64,
    /// The time for pruning the anchors, in seconds.
    pub pruning_seconds: R64,
    /// The time for initialising the chaining cost function, in seconds.
    pub cost_function_seconds: R64,
    /// The time spent in the chaining searches (all rounds), in seconds.
    pub chaining_seconds: R64,
    /// The time spent evaluating chains (all rounds), in seconds.
    pub evaluation_seconds: R64,
    /// The number of chaining rounds (chains computed).
    pub chaining_rounds: u64,
    /// The number of nodes opened by the chaining searches (all rounds).
    pub chaining_opened_nodes: u64,
    /// The number of nodes closed by the chaining searches (all rounds).
    pub chaining_closed_nodes: u64,
    /// The number of anchors removed by each pruning rule.
    /// Each removed anchor is attributed to exactly one rule.
    pub pruned_anchors: BTreeMap<String, u64>,
    /// The number of anchors matched by each pruning rule, independently of the other rules.
    pub matched_anchors: BTreeMap<String, u64>,
}

macro_rules! each_statistic {
    ($action:path) => {{
        $action!(cost);
        $action!(cost_per_base);
        $action!(duration_seconds);
        $action!(opened_nodes);
        $action!(closed_nodes);
        $action!(suboptimal_opened_nodes);
        $action!(suboptimal_opened_nodes_ratio);
        $action!(template_switch_amount);
        $action!(runtime);
        $action!(memory);
    }};
}

impl<Cost: AStarCost> AlignmentResult<Cost> {
    #[expect(clippy::too_many_arguments)]
    pub fn new_with_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        alignment: Vec<(usize, AlignmentType)>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        alignment_range: AlignmentRange,
        result: AStarResult<(), Cost>,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self::new(
            Some(alignment),
            reference,
            query,
            reference_name,
            query_name,
            alignment_range,
            result,
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            reference_length,
            query_length,
        )
    }

    #[expect(clippy::too_many_arguments)]
    pub fn new_without_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        result: AStarResult<(), Cost>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        alignment_range: AlignmentRange,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self::new(
            None,
            reference,
            query,
            reference_name,
            query_name,
            alignment_range,
            result,
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            reference_length,
            query_length,
        )
    }

    /// Create a result without alignment for an aligner that stopped because it reached its iteration limit.
    ///
    /// `cost_lower_bound` is a lower bound of the optimal cost known when the aligner stopped.
    /// See [`AlignmentStatistics::exceeded_iteration_limit`].
    #[expect(clippy::too_many_arguments)]
    pub fn new_exceeded_iteration_limit<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        iteration_limit: u64,
        cost_lower_bound: Cost,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        alignment_range: AlignmentRange,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        // `new` takes the cost from the result, which has no cost if there is no target.
        let mut result = Self::new(
            None,
            reference,
            query,
            reference_name,
            query_name,
            alignment_range,
            AStarResult::ExceededCostLimit {
                cost_limit: cost_lower_bound,
            },
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            reference_length,
            query_length,
        );
        let statistics = result.statistics_mut();
        statistics.result = AStarResult::NoTarget;
        statistics.exceeded_iteration_limit = Some(iteration_limit);
        result
    }

    #[expect(clippy::too_many_arguments)]
    fn new<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        alignment: Option<Vec<(usize, AlignmentType)>>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        alignment_range: AlignmentRange,
        result: AStarResult<(), Cost>,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        let cost = result.cost();
        let statistics = AlignmentStatistics {
            result,
            sequences: SequencePair::new(reference, query, reference_name, query_name),
            reference_offset: alignment_range.reference_offset(),
            query_offset: alignment_range.query_offset(),
            reference_limit: alignment_range.reference_limit(),
            query_limit: alignment_range.query_limit(),

            cost: (cost.as_f64()).try_into().unwrap(),
            cost_per_base: if reference_length + query_length == 0 {
                // Avoid NaN when aligning two empty sequences.
                r64(0.0)
            } else {
                ((cost.as_f64() * 2.0) / (reference_length + query_length) as f64)
                    .try_into()
                    .unwrap()
            },
            duration_seconds: duration_seconds.try_into().unwrap(),
            opened_nodes: (opened_nodes as f64).try_into().unwrap(),
            closed_nodes: (closed_nodes as f64).try_into().unwrap(),
            suboptimal_opened_nodes: (suboptimal_opened_nodes as f64).try_into().unwrap(),
            suboptimal_opened_nodes_ratio: if opened_nodes == 0 && suboptimal_opened_nodes == 0 {
                r64(0.0)
            } else {
                (suboptimal_opened_nodes as f64 / (opened_nodes - suboptimal_opened_nodes) as f64)
                    .try_into()
                    .unwrap()
            },
            template_switch_amount: r64(alignment
                .as_ref()
                .map(|alignment| {
                    alignment
                        .iter()
                        .filter(|(_, alignment_type)| alignment_type.is_template_switch_exit())
                        .count() as f64
                })
                .unwrap_or(0.0)),
            chain_cost: None,
            exceeded_iteration_limit: None,
            chain_statistics: None,
            runtime: r64(0.0),
            memory: r64(0.0),
        };

        if let Some(alignment) = alignment {
            Self::WithTarget {
                alignment: alignment.into(),
                statistics,
            }
        } else {
            Self::WithoutTarget { statistics }
        }
    }

    pub fn alignment_range(&self) -> AlignmentRange {
        let (reference_offset, query_offset, reference_limit, query_limit) = match self {
            Self::WithTarget { statistics, .. } | Self::WithoutTarget { statistics } => (
                statistics.reference_offset,
                statistics.query_offset,
                statistics.reference_limit,
                statistics.query_limit,
            ),
        };
        AlignmentRange::new(reference_offset, query_offset, reference_limit, query_limit)
    }
}

impl<Cost: AStarCost + From<u64>> AlignmentResult<Cost> {
    /// Heuristically extend the alignment beyond the alignment range as long as this does not increase the alignment cost.
    ///
    /// Returns the amount of steps the alignment was extended.
    #[allow(clippy::too_many_lines)]
    pub fn extend_beyond_range_without_increasing_cost<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        range: &mut AlignmentRange,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) -> usize {
        let Self::WithTarget {
            alignment,
            statistics,
        } = self
        else {
            trace!("There is no alignment, therefore we cannot postprocess it.");
            return 0;
        };
        if config.left_flank_length > 0 || config.right_flank_length > 0 {
            warn!("Alignment extension does not support flanks");
            return 0;
        }

        // Compute cost before extending.
        let mut current_cost = alignment.compute_complete_cost(
            reference,
            query,
            range.reference_offset(),
            range.query_offset(),
            config,
        );
        let alignment_cost = (statistics.cost.round().raw() as u64).into();
        debug_assert_eq!(
            current_cost,
            alignment_cost,
            "computed cost {current_cost} != alignment cost {alignment_cost}; {}",
            alignment.cigar()
        );
        let mut extension_steps = 0;

        // Extend left with equal cost.
        while let Some(new_range) = range.move_offsets_left() {
            // Determine extension alignment type.
            let reference_character = reference[new_range.reference_offset()].clone();
            let query_character = query[new_range.query_offset()].clone();
            let alignment_type = if reference_character == query_character {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert alignment type at beginning of alignment.
            if alignment
                .inner_mut()
                .first()
                .is_some_and(|(_, ty)| *ty == alignment_type)
            {
                alignment.inner_mut()[0].0 += 1;
            } else {
                alignment.inner_mut().insert(0, (1, alignment_type));
            }

            // Compute cost.
            let new_cost = alignment.compute_complete_cost(
                reference,
                query,
                new_range.reference_offset(),
                new_range.query_offset(),
                config,
            );

            trace!("Extending left to {new_range}. Cost {current_cost} -> {new_cost}");

            if new_cost > current_cost {
                // If cost is increasing, revert alignment and break.
                alignment.inner_mut()[0].0 -= 1;
                if alignment.inner_mut()[0].0 == 0 {
                    alignment.inner_mut().remove(0);
                }
                break;
            }
            current_cost = new_cost;
            *range = new_range;
            extension_steps += 1;
        }

        // Extend right with equal cost.
        while range.reference_limit() < reference.len() && range.query_limit() < query.len() {
            // Update ranges.
            let Some(new_range) = range.move_limits_right() else {
                // Never occurs due to loop invariant
                break;
            };

            // Determine extension alignment type.
            let reference_character = reference[new_range.reference_limit() - 1].clone();
            let query_character = query[new_range.query_limit() - 1].clone();
            let alignment_type = if reference_character == query_character {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert alignment type at end of alignment.
            if alignment
                .inner_mut()
                .last()
                .is_some_and(|(_, ty)| *ty == alignment_type)
            {
                let Some((n, _)) = alignment.inner_mut().last_mut() else {
                    unreachable!("Due to if condition")
                };
                *n += 1;
            } else {
                alignment.inner_mut().push((1, alignment_type));
            }

            // Compute cost.
            let new_cost = alignment.compute_complete_cost(
                reference,
                query,
                new_range.reference_offset(),
                new_range.query_offset(),
                config,
            );

            trace!("Extending right to {new_range}. Cost {current_cost} -> {new_cost}");

            if new_cost > current_cost {
                // If cost is increasing, revert alignment and break.
                let Some((n, _)) = alignment.inner_mut().last_mut() else {
                    unreachable!("As we have just added at least one element")
                };
                *n -= 1;
                if *n == 0 {
                    alignment.inner_mut().pop();
                }
                break;
            }
            current_cost = new_cost;
            *range = new_range;
            extension_steps += 1;
        }

        // Update statistics from updated range
        statistics.reference_offset = range.reference_offset();
        statistics.query_offset = range.query_offset();
        statistics.reference_limit = range.reference_limit();
        statistics.query_limit = range.query_limit();
        extension_steps
    }

    #[allow(clippy::too_many_lines)]
    pub fn compute_ts_uncertainty_ranges<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        range: &AlignmentRange,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
        mode: TSMUncertaintyRangeExtensionMode,
    ) {
        let Self::WithTarget {
            alignment,
            statistics,
        } = self
        else {
            trace!("There is no alignment, therefore we cannot postprocess it.");
            return;
        };
        if mode == TSMUncertaintyRangeExtensionMode::None {
            debug!("TSM uncertainty range extension is disabled");
            return;
        }
        if config.left_flank_length > 0 || config.right_flank_length > 0 {
            warn!("TS extension does not support flanks");
            return;
        }

        let reference_offset = range.reference_offset();
        let query_offset = range.query_offset();
        let ignore_geometry_cost = match mode {
            TSMUncertaintyRangeExtensionMode::None => unreachable!(),
            TSMUncertaintyRangeExtensionMode::EqualCost => false,
            TSMUncertaintyRangeExtensionMode::EqualCostIgnoreGeometry => true,
        };

        for i in 0..alignment.inner_mut().len() {
            let alignment_type = alignment.inner_mut()[i].1;

            if let AlignmentType::TemplateSwitchEntrance {
                uncertainty_range: None,
                ..
            } = alignment_type
            {
                let mut uncertainty_range = TSMUncertaintyRange::new();

                // Keep track of the current cost, which might get less while extending but never increase
                let mut current_cost = alignment.compute_cost(
                    reference,
                    query,
                    reference_offset,
                    query_offset,
                    config,
                    ignore_geometry_cost,
                );
                if !ignore_geometry_cost {
                    debug_assert_eq!(current_cost, (statistics.cost.round().raw() as u64).into());
                }

                // Move start to the left
                {
                    let mut min_start_alignment = alignment.clone();
                    let mut i = i;

                    while min_start_alignment.move_template_switch_start_backwards(
                        reference,
                        query,
                        reference_offset,
                        query_offset,
                        &mut i,
                    ) {
                        let Some(new_start_left_shift) =
                            uncertainty_range.start_left_shift.checked_add(1)
                        else {
                            trace!(
                                "Stopping moving TS start backwards because min_start index is exhausted."
                            );
                            break;
                        };
                        let new_cost = min_start_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                            ignore_geometry_cost,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS start backwards because cost {new_cost} > {current_cost} with alignment {min_start_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        uncertainty_range.start_left_shift = new_start_left_shift;
                    }
                }

                // Move start to the right
                {
                    let mut max_start_alignment = alignment.clone();
                    let mut i = i;

                    while max_start_alignment.move_template_switch_start_forwards(
                        reference,
                        query,
                        reference_offset,
                        query_offset,
                        &mut i,
                    ) {
                        let Some(new_start_right_shift) =
                            uncertainty_range.start_right_shift.checked_add(1)
                        else {
                            trace!(
                                "Stopping moving TS start forwards because max_start index is exhausted."
                            );
                            break;
                        };
                        let new_cost = max_start_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                            ignore_geometry_cost,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS start forwards because cost {new_cost} > {current_cost} with alignment {max_start_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        uncertainty_range.start_right_shift = new_start_right_shift;
                    }
                }

                // Move end to the left
                {
                    let mut min_end_alignment = alignment.clone();
                    while min_end_alignment.move_template_switch_end_backwards(
                        reference,
                        query,
                        reference_offset,
                        query_offset,
                        i,
                    ) {
                        let Some(new_end_left_shift) =
                            uncertainty_range.end_left_shift.checked_add(1)
                        else {
                            trace!(
                                "Stopping moving TS end backwards because min_end index is exhausted."
                            );
                            break;
                        };
                        let new_cost = min_end_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                            ignore_geometry_cost,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS end backwards because cost {new_cost} > {current_cost} with alignment {min_end_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        uncertainty_range.end_left_shift = new_end_left_shift;
                    }
                }

                // Move end to the right
                {
                    let mut max_end_alignment = alignment.clone();
                    while max_end_alignment.move_template_switch_end_forwards(
                        reference,
                        query,
                        reference_offset,
                        query_offset,
                        i,
                    ) {
                        let Some(new_end_right_shift) =
                            uncertainty_range.end_right_shift.checked_add(1)
                        else {
                            trace!(
                                "Stopping moving TS end forwards because max_end index is exhausted."
                            );
                            break;
                        };
                        let new_cost = max_end_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                            ignore_geometry_cost,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS end forwards because cost {new_cost} > {current_cost} with alignment {max_end_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        uncertainty_range.end_right_shift = new_end_right_shift;
                    }
                }

                let AlignmentType::TemplateSwitchEntrance {
                    uncertainty_range: alignment_uncertainty_range,
                    ..
                } = &mut alignment.inner_mut()[i].1
                else {
                    unreachable!()
                };
                *alignment_uncertainty_range = Some(uncertainty_range);
            }
        }
    }
}

impl<Cost> AlignmentResult<Cost> {
    pub fn statistics(&self) -> &AlignmentStatistics<Cost> {
        match self {
            AlignmentResult::WithTarget { statistics, .. } => statistics,
            AlignmentResult::WithoutTarget { statistics } => statistics,
        }
    }

    pub fn statistics_mut(&mut self) -> &mut AlignmentStatistics<Cost> {
        match self {
            AlignmentResult::WithTarget { statistics, .. } => statistics,
            AlignmentResult::WithoutTarget { statistics } => statistics,
        }
    }
}

impl<Cost> AlignmentResult<Cost> {
    pub fn cigar(&self) -> String {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl Write) -> Result {
        let Self::WithTarget { alignment, .. } = self else {
            return Ok(());
        };

        alignment.write_cigar(writer)
    }
}

impl<Cost: Clone> AlignmentStatistics<Cost> {
    pub fn min_value() -> Self {
        let mut result = Self::default();

        macro_rules! create_min_value {
            ($field:ident) => {
                result.$field = R64::min_value();
            };
        }
        each_statistic!(create_min_value);

        result
    }

    pub fn max_value() -> Self {
        let mut result = Self::default();

        macro_rules! create_max_value {
            ($field:ident) => {
                result.$field = R64::max_value();
            };
        }
        each_statistic!(create_max_value);

        result
    }

    pub fn zero() -> Self {
        let mut result = Self::default();

        macro_rules! create_zero {
            ($field:ident) => {
                result.$field = R64::zero();
            };
        }
        each_statistic!(create_zero);

        result
    }

    pub fn piecewise_min(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! piecewise_min {
            ($field:ident) => {
                result.$field = result.$field.min(other.$field);
            };
        }
        each_statistic!(piecewise_min);

        result
    }

    pub fn piecewise_max(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! piecewise_max {
            ($field:ident) => {
                result.$field = result.$field.max(other.$field);
            };
        }
        each_statistic!(piecewise_max);

        result
    }

    pub fn piecewise_add(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! add {
            ($field:ident) => {
                result.$field += other.$field;
            };
        }
        each_statistic!(add);

        result
    }

    pub fn piecewise_div(&self, divisor: R64) -> Self {
        let mut result = self.clone();

        macro_rules! div {
            ($field:ident) => {
                result.$field /= divisor;
            };
        }
        each_statistic!(div);

        result
    }

    pub fn piecewise_percentile(statistics: &[Self], percentile: R64) -> Self {
        assert!(percentile >= 0.0);
        assert!(percentile <= 1.0);
        assert!(!statistics.is_empty());
        let mut result = Self::zero();
        let mut buffer = Vec::new();

        macro_rules! percentile {
            ($field:ident) => {
                buffer.clear();
                buffer.extend(statistics.iter().map(|s| s.$field));
                buffer.sort();

                // Scale percentile.
                let index = (percentile * buffer.len() as f64).floor().raw() as usize;
                let index = if index == buffer.len() {
                    // Edge case if percentile == 1.0
                    index - 1
                } else {
                    index
                };
                result.$field = buffer[index];
            };
        }
        each_statistic!(percentile);

        result
    }
}

impl<Cost: Display> Display for AlignmentResult<Cost> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if let AlignmentResult::WithTarget { .. } = self {
            write!(f, "CIGAR: ")?;
            self.write_cigar(f)?;
            writeln!(f)?;
        } else {
            writeln!(f, "No alignment found")?;
        }

        let (AlignmentResult::WithTarget { statistics, .. }
        | AlignmentResult::WithoutTarget { statistics }) = self;

        statistics.fmt(f)
    }
}

impl<Cost: Display> Display for AlignmentStatistics<Cost> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if let Some(iteration_limit) = self.exceeded_iteration_limit {
            writeln!(
                f,
                "Exceeded iteration limit of {iteration_limit}, the chain cost is at least {}",
                self.cost
            )?;
        } else {
            writeln!(f, "{}", self.result)?;
        }
        if let Some(chain_cost) = self.chain_cost {
            writeln!(f, "Chain cost: {chain_cost}")?;
        }
        writeln!(f, "Reference offset: {}", self.reference_offset)?;
        writeln!(f, "Query offset: {}", self.query_offset)?;
        writeln!(f, "Cost per base: {:.2}", self.cost_per_base)?;
        writeln!(f, "Opened nodes: {}", self.opened_nodes)?;
        writeln!(f, "Closed nodes: {}", self.closed_nodes)?;
        writeln!(
            f,
            "Suboptimal openend nodes: {}",
            self.suboptimal_opened_nodes
        )?;
        writeln!(
            f,
            "Suboptimal openend nodes per optimal opened node: {:.2}",
            self.suboptimal_opened_nodes_ratio
        )?;
        write!(f, "Duration: {:.2}s", self.duration_seconds)?;

        Ok(())
    }
}

impl<Cost> Default for AlignmentStatistics<Cost> {
    fn default() -> Self {
        Self {
            result: Default::default(),
            sequences: Default::default(),
            reference_offset: Default::default(),
            query_offset: Default::default(),
            reference_limit: usize::MAX,
            query_limit: usize::MAX,
            cost: Default::default(),
            cost_per_base: Default::default(),
            duration_seconds: Default::default(),
            opened_nodes: Default::default(),
            closed_nodes: Default::default(),
            suboptimal_opened_nodes: Default::default(),
            suboptimal_opened_nodes_ratio: Default::default(),
            template_switch_amount: Default::default(),
            chain_cost: None,
            exceeded_iteration_limit: None,
            chain_statistics: None,
            runtime: Default::default(),
            memory: Default::default(),
        }
    }
}
