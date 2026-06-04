use std::fmt::{Display, Formatter, Result, Write};

use a_star_sequences::SequencePair;
use alignment::Alignment;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStarResult, cost::AStarCost};
use log::{trace, warn};
use noisy_float::types::{R64, r64};
use num_traits::{Float, Zero};

use crate::{
    a_star_aligner::{
        alignment_geometry::AlignmentCoordinates, template_switch_distance::AlignmentType,
    },
    config::TemplateSwitchConfig,
};

use super::alignment_geometry::AlignmentRange;

pub mod a_star_sequences;
pub mod alignment;

pub trait IAlignmentType {
    fn is_repeatable(&self) -> bool;

    fn is_repeated(&self, previous: &Self) -> bool;

    fn is_internal(&self) -> bool;

    fn is_template_switch_entrance(&self) -> bool;

    fn is_template_switch_exit(&self) -> bool;

    fn alternative_start(&self) -> Option<AlignmentCoordinates>;

    fn alternative_end(&self) -> Option<AlignmentCoordinates>;
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
pub enum AlignmentResult<AlignmentType, Cost> {
    WithTarget {
        #[cfg_attr(feature = "serde", serde(flatten))]
        alignment: Alignment<AlignmentType>,

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

impl<AlignmentType: IAlignmentType, Cost: AStarCost> AlignmentResult<AlignmentType, Cost> {
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
            cost_per_base: ((cost.as_f64() * 2.0) / (reference_length + query_length) as f64)
                .try_into()
                .unwrap(),
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

impl<Cost: AStarCost + From<u64>>
    AlignmentResult<super::template_switch_distance::AlignmentType, Cost>
{
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
        let mut current_cost = alignment.compute_cost(
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
            let new_cost = alignment.compute_cost(
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
            let new_cost = alignment.compute_cost(
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
        extension_steps
    }

    #[allow(clippy::too_many_lines)]
    pub fn compute_ts_equal_cost_ranges<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        range: &AlignmentRange,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) {
        let Self::WithTarget {
            alignment,
            statistics,
        } = self
        else {
            trace!("There is no alignment, therefore we cannot postprocess it.");
            return;
        };
        if config.left_flank_length > 0 || config.right_flank_length > 0 {
            warn!("TS extension does not support flanks");
            return;
        }

        let reference_offset = range.reference_offset();
        let query_offset = range.query_offset();

        for i in 0..alignment.inner_mut().len() {
            let alignment_type = alignment.inner_mut()[i].1;

            if let super::template_switch_distance::AlignmentType::TemplateSwitchEntrance {
                mut equal_cost_range,
                ..
            } = alignment_type
            {
                equal_cost_range.min_start = 0;
                equal_cost_range.max_start = 0;
                equal_cost_range.min_end = 0;
                equal_cost_range.max_end = 0;

                // Keep track of the current cost, which might get less while extending but never increase
                let mut current_cost = alignment.compute_cost(
                    reference,
                    query,
                    reference_offset,
                    query_offset,
                    config,
                );
                debug_assert_eq!(current_cost, (statistics.cost.round().raw() as u64).into());

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
                        let new_cost = min_start_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS start backwards because cost {new_cost} > {current_cost} with alignment {min_start_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        equal_cost_range.min_start -= 1;
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
                        let new_cost = max_start_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS start forwards because cost {new_cost} > {current_cost} with alignment {max_start_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        equal_cost_range.max_start += 1;
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
                        let new_cost = min_end_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS end backwards because cost {new_cost} > {current_cost} with alignment {min_end_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        equal_cost_range.min_end -= 1;
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
                        let new_cost = max_end_alignment.compute_cost(
                            reference,
                            query,
                            reference_offset,
                            query_offset,
                            config,
                        );
                        if new_cost > current_cost {
                            trace!(
                                "Stopping moving TS end forwards because cost {new_cost} > {current_cost} with alignment {max_end_alignment}"
                            );
                            break;
                        }
                        current_cost = new_cost;
                        equal_cost_range.max_end += 1;
                    }
                }

                let super::template_switch_distance::AlignmentType::TemplateSwitchEntrance {
                    equal_cost_range: alignment_equal_cost_range,
                    ..
                } = &mut alignment.inner_mut()[i].1
                else {
                    unreachable!()
                };
                *alignment_equal_cost_range = equal_cost_range;
            }
        }
    }
}

impl<AlignmentType, Cost> AlignmentResult<AlignmentType, Cost> {
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

impl<AlignmentType: IAlignmentType, Cost> AlignmentResult<AlignmentType, Cost> {
    pub fn cigar(&self) -> String
    where
        AlignmentType: Display,
    {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl Write) -> Result
    where
        AlignmentType: Display,
    {
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

impl<AlignmentType: Display + IAlignmentType, Cost: Display> Display
    for AlignmentResult<AlignmentType, Cost>
{
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
        writeln!(f, "{}", self.result)?;
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
            runtime: Default::default(),
            memory: Default::default(),
        }
    }
}
