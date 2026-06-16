use std::{fmt::Debug, str::FromStr};

use compact_genome::{
    implementation::{alphabets::dna_alphabet_or_n::DnaAlphabetOrN, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U64Cost;
use traitsequence::interface::Sequence;

use crate::{
    a_star_aligner::{
        template_switch_distance::{
            context::DynamicStrategies,
            strategies::{
                AlignmentStrategySelection,
                allow_ts_14_out_of_range::{
                    AdditionalExplicitTSMStartsAndEnds, Ts14OutOfRangeStrategy,
                },
                descendant::{
                    AnyTemplateSwitchDescendantStrategy, OnlyEqualTemplateSwitchDescendantStrategy,
                    TemplateSwitchDescendantStrategy,
                },
                primary_range::RangePrunePrimaryRangeStrategy,
                secondary_deletion::AllowSecondaryDeletionStrategy,
                shortcut::NoShortcutStrategy,
                template_switch_min_length::{
                    PreprocessedLookaheadTemplateSwitchMinLengthStrategy,
                    PreprocessedTemplateSwitchMinLengthStrategy,
                },
            },
        },
        template_switch_distance_a_star_align,
    },
    config::TemplateSwitchConfig,
};

use super::{
    alignment_geometry::AlignmentRange,
    alignment_result::AlignmentResult,
    template_switch_distance::{
        AlignmentType,
        strategies::{
            chaining::{ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy},
            node_ord::AntiDiagonalNodeOrdStrategy,
            primary_match::AllowPrimaryMatchStrategy,
            template_switch_count::{
                MaxTemplateSwitchCountStrategy, NoTemplateSwitchCountStrategy,
                TemplateSwitchCountStrategy,
            },
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                TemplateSwitchMinLengthStrategy,
            },
            template_switch_total_length::{
                MaxTemplateSwitchTotalLengthStrategy, NoTemplateSwitchTotalLengthStrategy,
                TemplateSwitchTotalLengthStrategy,
            },
        },
    },
};

pub use compact_genome::implementation::alphabets;
pub use compact_genome::interface::alphabet::Alphabet;

#[cfg(test)]
mod tests;

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum MinLengthStrategySelector {
    None,
    #[default]
    Lookahead,
    /// Check for perfect matches of minimum length and compute a lower bound for imperfect matches.
    PreprocessPrice,
    /// Check for perfect matches of minimum length and filter out imperfect matches.
    PreprocessFilter,
    /// Check for perfect matches of minimum length and cache the costs of imperfect matches.
    PreprocessLookahead,
}

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum ChainingStrategySelector {
    #[default]
    None,
    // PrecomputeOnly,
    LowerBound,
}

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum TotalLengthStrategySelector {
    None,
    #[default]
    Maximise,
}

/// Select the descendant strategy.
#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum DescendantStrategySelector {
    /// The TSM descendant can be any of the two input sequences for every TSM.
    #[default]
    AllowAny,
    /// All TSMs must have the same descendant, but it can be either of the two input sequences.
    AllowOnlyAllEqual,
}

/// Just used in this file to bundle the query parameters to make the code more readable
struct QueryData<'a> {
    reference_name: &'a str,
    reference: &'a [u8],
    query_name: &'a str,
    query: &'a [u8],
    ranges: Option<AlignmentRange>,
    gap_characters: &'a [u8],
    cost_limit: Option<u64>,
    memory_limit: Option<usize>,
    extend_beyond_range: bool,
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default, bound = ""))]
pub struct Aligner<AlphabetType: Alphabet = DnaAlphabetOrN> {
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,

    // ↓ Settings for how the alignment (definition)
    min_length_strategy: MinLengthStrategySelector,
    chaining_strategy: ChainingStrategySelector,
    total_length_strategy: TotalLengthStrategySelector,
    descendant_strategy: DescendantStrategySelector,
    ts_14_out_of_range_strategy: Ts14OutOfRangeStrategy,
    no_ts: bool,
}

impl Aligner<DnaAlphabetOrN> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// We implement this ourselves to avoid the overly restrictive trait bound
/// `AlphabetType: Default` that would be generated with the derive macro
impl<AlphabetType: Alphabet> Default for Aligner<AlphabetType> {
    fn default() -> Self {
        Self {
            costs: TemplateSwitchConfig::<AlphabetType, U64Cost>::default(),
            min_length_strategy: Default::default(),
            chaining_strategy: Default::default(),
            total_length_strategy: Default::default(),
            descendant_strategy: Default::default(),
            ts_14_out_of_range_strategy: Default::default(),
            no_ts: false,
        }
    }
}

impl<AlphabetType: Alphabet> Aligner<AlphabetType> {
    /// Parse a cost string and update the aligner cost
    ///
    /// # Errors
    /// Returns an error if the string cannot be parsed as `TemplateSwitchConfig<AlphabetType>`
    pub fn set_costs_parse(&mut self, costs: &str) -> Result<&mut Self, crate::error::Error> {
        self.costs = TemplateSwitchConfig::from_str(costs)?;
        Ok(self)
    }

    pub fn set_costs(&mut self, costs: TemplateSwitchConfig<AlphabetType, U64Cost>) -> &mut Self {
        self.costs = costs;
        self
    }

    pub fn set_min_length_strategy(
        &mut self,
        min_length_strategy: MinLengthStrategySelector,
    ) -> &mut Self {
        self.min_length_strategy = min_length_strategy;
        self
    }

    pub fn set_chaining_strategy(
        &mut self,
        chaining_strategy: ChainingStrategySelector,
    ) -> &mut Self {
        self.chaining_strategy = chaining_strategy;
        self
    }

    pub fn set_total_length_strategy(
        &mut self,
        total_length_strategy: TotalLengthStrategySelector,
    ) -> &mut Self {
        self.total_length_strategy = total_length_strategy;
        self
    }

    pub fn set_descendant_strategy(
        &mut self,
        descendant_strategy: DescendantStrategySelector,
    ) -> &mut Self {
        self.descendant_strategy = descendant_strategy;
        self
    }

    pub fn set_no_ts(&mut self, no_ts: bool) -> &mut Self {
        self.no_ts = no_ts;
        self
    }

    pub fn set_ts_14_out_of_range_strategy(
        &mut self,
        ts_14_out_of_range_strategy: Ts14OutOfRangeStrategy,
    ) -> &mut Self {
        self.ts_14_out_of_range_strategy = ts_14_out_of_range_strategy;
        self
    }

    /// Perform the actual alignment.
    ///
    /// Special gap characters can be specified in order to interpret the given sequences as an alignment with gaps instead of raw sequences.
    /// In this case, the range coordinates must refer to the ungapped sequences.
    ///
    /// Depending on the inputs, this takes quite long.
    /// Consider spawning tasks to not block the user application.
    #[allow(
        clippy::too_many_arguments,
        reason = "For now, this is the intended API"
    )]
    #[must_use]
    pub fn align(
        &self,
        reference_name: &str,
        reference: &[u8],
        query_name: &str,
        query: &[u8],
        ranges: Option<AlignmentRange>,
        gap_characters: &[u8],
        cost_limit: Option<u64>,
        memory_limit: Option<usize>,
        extend_beyond_range: bool,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        let data = QueryData {
            reference_name,
            reference,
            query_name,
            query,
            ranges,
            gap_characters,
            cost_limit,
            memory_limit,
            extend_beyond_range,
        };
        self.align_select_min_length_strategy(data)
    }

    fn align_select_min_length_strategy(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.min_length_strategy {
            MinLengthStrategySelector::None => self
                .align_select_chaining_strategy::<NoTemplateSwitchMinLengthStrategy<U64Cost>>(data),
            MinLengthStrategySelector::Lookahead => self.align_select_chaining_strategy::<LookaheadTemplateSwitchMinLengthStrategy<U64Cost>>(data),
            MinLengthStrategySelector::PreprocessPrice => self.align_select_chaining_strategy::<PreprocessedTemplateSwitchMinLengthStrategy<false, U64Cost>>(data),
            MinLengthStrategySelector::PreprocessFilter => self.align_select_chaining_strategy::<PreprocessedTemplateSwitchMinLengthStrategy<true, U64Cost>>(data),
            MinLengthStrategySelector::PreprocessLookahead => self.align_select_chaining_strategy::<PreprocessedLookaheadTemplateSwitchMinLengthStrategy<U64Cost>>(data),
        }
    }

    fn align_select_chaining_strategy<ML: TemplateSwitchMinLengthStrategy<U64Cost>>(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.chaining_strategy {
            ChainingStrategySelector::None => {
                self.align_select_total_length_strategy::<ML, NoChainingStrategy<U64Cost>>(data)
            }
            ChainingStrategySelector::LowerBound => self
                .align_select_total_length_strategy::<ML, LowerBoundChainingStrategy<U64Cost>>(
                    data,
                ),
        }
    }

    fn align_select_total_length_strategy<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
    >(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.total_length_strategy {
            TotalLengthStrategySelector::None => self
                .align_select_no_ts_strategy::<ML, CH, NoTemplateSwitchTotalLengthStrategy>(data),
            TotalLengthStrategySelector::Maximise => self
                .align_select_no_ts_strategy::<ML, CH, MaxTemplateSwitchTotalLengthStrategy>(data),
        }
    }

    fn align_select_no_ts_strategy<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
    >(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        if self.no_ts {
            self.align_select_descendant_strategy::<ML, CH, TL, MaxTemplateSwitchCountStrategy>(
                data, 0,
            )
        } else {
            self.align_select_descendant_strategy::<ML, CH, TL, NoTemplateSwitchCountStrategy>(
                data,
                (),
            )
        }
    }

    fn align_select_descendant_strategy<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
    >(
        &self,
        data: QueryData,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.descendant_strategy {
            DescendantStrategySelector::AllowAny => self
                .align_call::<ML, CH, TL, TC, AnyTemplateSwitchDescendantStrategy>(
                    data,
                    count_strategy_memory,
                ),
            DescendantStrategySelector::AllowOnlyAllEqual => {
                self.align_call::<ML, CH, TL, TC, OnlyEqualTemplateSwitchDescendantStrategy>(
                    data,
                    count_strategy_memory,
                )
            }
        }
    }

    fn align_call<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
        DS: TemplateSwitchDescendantStrategy,
    >(
        &self,
        data: QueryData,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        // TODO error handling
        let reference = VectorGenome::<AlphabetType>::from_iter_u8(
            data.reference
                .iter()
                .copied()
                .filter(|c| !data.gap_characters.contains(c)),
        )
        .unwrap();
        let query = VectorGenome::from_iter_u8(
            data.query
                .iter()
                .copied()
                .filter(|c| !data.gap_characters.contains(c)),
        )
        .unwrap();

        let range = data
            .ranges
            .unwrap_or_else(|| AlignmentRange::new_complete(reference.len(), query.len()));
        let additional_tsm_starts_and_ends =
            if self.ts_14_out_of_range_strategy == Ts14OutOfRangeStrategy::Allow {
                AdditionalExplicitTSMStartsAndEnds::new(
                    &data
                        .reference
                        .iter()
                        .map(|c| *c as char)
                        .collect::<String>(),
                    &data.query.iter().map(|c| *c as char).collect::<String>(),
                    &range,
                    &data
                        .gap_characters
                        .iter()
                        .map(|c| *c as char)
                        .collect::<Vec<_>>(),
                    false,
                )
                .unwrap()
            } else {
                AdditionalExplicitTSMStartsAndEnds::default()
            };

        let cost_limit = data.cost_limit.map(U64Cost::from);

        template_switch_distance_a_star_align::<
            AlignmentStrategySelection<
                AlphabetType,
                U64Cost,
                AntiDiagonalNodeOrdStrategy,
                ML,
                CH,
                TC,
                AllowSecondaryDeletionStrategy,
                NoShortcutStrategy<U64Cost>,
                AllowPrimaryMatchStrategy,
                RangePrunePrimaryRangeStrategy,
                TL,
                DS,
            >,
            _,
        >(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            data.reference_name,
            data.query_name,
            range,
            additional_tsm_starts_and_ends,
            &self.costs,
            DynamicStrategies {
                ts_14_out_of_range: self.ts_14_out_of_range_strategy,
            },
            cost_limit,
            data.memory_limit,
            false,
            data.extend_beyond_range,
            count_strategy_memory,
        )
    }
}
