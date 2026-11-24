use std::{fmt::Debug, str::FromStr};

use compact_genome::{
    implementation::{alphabets::dna_alphabet_or_n::DnaAlphabetOrN, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U64Cost;
use traitsequence::interface::Sequence;

use crate::{
    aligner::{
        context::DynamicStrategies,
        strategies::{
            AlignmentStrategySelection,
            descendant::{
                AnyTemplateSwitchDescendantStrategy, OnlyEqualTemplateSwitchDescendantStrategy,
                TemplateSwitchDescendantStrategy,
            },
            primary_range::{
                NoPrunePrimaryRangeStrategy, PrimaryRangeStrategy, RangePrunePrimaryRangeStrategy,
            },
            secondary_deletion::AllowSecondaryDeletionStrategy,
            shortcut::NoShortcutStrategy,
            template_switch_min_length::{
                PreprocessedLookaheadTemplateSwitchMinLengthStrategy,
                PreprocessedTemplateSwitchMinLengthStrategy,
            },
        },
        template_switch_distance_a_star_align,
        ts_14_out_of_range::{AdditionalExplicitTSMStartsAndEnds, Ts14OutOfRangeStrategy},
    },
    alignment_result::TSMUncertaintyRangeExtensionMode,
    config::TemplateSwitchConfig,
};

use crate::{
    aligner::strategies::{
        chaining::{ChainingStrategy, NoChainingStrategy},
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
    alignment_geometry::AlignmentRange,
    alignment_result::AlignmentResult,
};

pub use compact_genome::implementation::alphabets;
pub use compact_genome::interface::alphabet::Alphabet;

#[cfg(test)]
use crate::aligner::strategies::chaining::LowerBoundChainingStrategy;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod validation_tests;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
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

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum ChainingStrategySelector {
    #[default]
    None,
    // PrecomputeOnly,
    /// Use [`LowerBoundChainingStrategy`].
    ///
    /// **Disabled** (issue #248): this variant exists only in tests of this crate,
    /// so it cannot be selected through the public API or deserialised (e.g. from Python).
    /// The memory and time of the precomputation grow with about n^4.8, and the chain lower bound was 0 in practice.
    /// To re-enable it, remove the `#[cfg(test)]` attributes here, on its match arm in `align_select_chaining_strategy`
    /// and on the import of `LowerBoundChainingStrategy`, and drop the `ValueError` check in the Python wrapper
    /// (`python_bindings/python/tsalign/__init__.py`).
    /// The CLI (`--ts-chaining-strategy lower-bound`) and serde then pick up the variant automatically.
    #[cfg(test)]
    LowerBound,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum TotalLengthStrategySelector {
    None,
    #[default]
    Maximise,
}

/// Select the descendant strategy.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum DescendantStrategySelector {
    /// The TSM descendant can be any of the two input sequences for every TSM.
    #[default]
    AllowAny,
    /// All TSMs must have the same descendant, but it can be either of the two input sequences.
    AllowOnlyAllEqual,
}

/// Select whether the primary alignment is restricted to the alignment range.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum PrimaryRangeStrategySelector {
    /// The primary alignment is restricted to the alignment range.
    #[default]
    Prune,
    /// The primary alignment is restricted to the input sequences only, and may hence leave the
    /// alignment range.
    NoPrune,
}

/// Just used in this file to bundle the query parameters to make the code more readable
struct QueryData<'a, SubsequenceType: ?Sized> {
    reference_name: &'a str,
    reference: &'a SubsequenceType,
    query_name: &'a str,
    query: &'a SubsequenceType,
    range: AlignmentRange,
    additional_tsm_starts_and_ends: AdditionalExplicitTSMStartsAndEnds,
    cost_limit: Option<U64Cost>,
    memory_limit: Option<usize>,
    extend_beyond_range: bool,
    tsm_uncertainty_range_extension_mode: TSMUncertaintyRangeExtensionMode,
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
    primary_range_strategy: PrimaryRangeStrategySelector,
    no_ts: bool,

    /// Force the A* search to run in label-correcting mode. This is for debug purposes only.
    force_label_correcting: bool,
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
            primary_range_strategy: Default::default(),
            no_ts: false,
            force_label_correcting: false,
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

    pub fn set_primary_range_strategy(
        &mut self,
        primary_range_strategy: PrimaryRangeStrategySelector,
    ) -> &mut Self {
        self.primary_range_strategy = primary_range_strategy;
        self
    }

    /// Force the A* search to run in label-correcting mode.
    ///
    /// This is for debug purposes only.
    pub fn set_force_label_correcting(&mut self, force_label_correcting: bool) -> &mut Self {
        self.force_label_correcting = force_label_correcting;
        self
    }

    /// Perform the actual alignment.
    ///
    /// Special gap characters can be specified in order to interpret the given sequences as an alignment with gaps instead of raw sequences.
    /// In this case, the range coordinates must refer to the ungapped sequences.
    ///
    /// The `cost_limit` is inclusive: if there is no alignment with cost at most `cost_limit`, then no alignment is returned.
    /// The `memory_limit` is an approximate memory limit in bytes; if it is exceeded, then no alignment is returned.
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
        tsm_uncertainty_range_extension_mode: TSMUncertaintyRangeExtensionMode,
    ) -> AlignmentResult<U64Cost> {
        // TODO error handling
        let reference_genome = VectorGenome::<AlphabetType>::from_iter_u8(
            reference
                .iter()
                .copied()
                .filter(|c| !gap_characters.contains(c)),
        )
        .unwrap();
        let query_genome = VectorGenome::from_iter_u8(
            query
                .iter()
                .copied()
                .filter(|c| !gap_characters.contains(c)),
        )
        .unwrap();

        let range = ranges.unwrap_or_else(|| {
            AlignmentRange::new_complete(reference_genome.len(), query_genome.len())
        });
        let additional_tsm_starts_and_ends =
            if self.ts_14_out_of_range_strategy == Ts14OutOfRangeStrategy::Allow {
                AdditionalExplicitTSMStartsAndEnds::new(
                    &reference.iter().map(|c| *c as char).collect::<String>(),
                    &query.iter().map(|c| *c as char).collect::<String>(),
                    &range,
                    &gap_characters
                        .iter()
                        .map(|c| *c as char)
                        .collect::<Vec<_>>(),
                    false,
                )
                .unwrap()
            } else {
                AdditionalExplicitTSMStartsAndEnds::default()
            };

        self.align_sequences(
            reference_name,
            reference_genome.as_genome_subsequence(),
            query_name,
            query_genome.as_genome_subsequence(),
            range,
            additional_tsm_starts_and_ends,
            cost_limit.map(U64Cost::from),
            memory_limit,
            extend_beyond_range,
            tsm_uncertainty_range_extension_mode,
        )
    }

    /// Align two sequences that are already free of gap characters.
    ///
    /// This is the lower-level variant of [`Self::align`] for callers that parse and filter the input themselves, like the CLI.
    ///
    /// * `range` is the alignment range on the given sequences.
    /// * `additional_tsm_starts_and_ends` are the TSM starts and ends outside of `range` that are allowed if the
    ///   [`Ts14OutOfRangeStrategy`] is [`Allow`](Ts14OutOfRangeStrategy::Allow).
    ///   Pass [`AdditionalExplicitTSMStartsAndEnds::default()`] otherwise.
    ///
    /// The other parameters are as in [`Self::align`].
    #[allow(
        clippy::too_many_arguments,
        reason = "For now, this is the intended API"
    )]
    #[must_use]
    pub fn align_sequences<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference_name: &str,
        reference: &SubsequenceType,
        query_name: &str,
        query: &SubsequenceType,
        range: AlignmentRange,
        additional_tsm_starts_and_ends: AdditionalExplicitTSMStartsAndEnds,
        cost_limit: Option<U64Cost>,
        memory_limit: Option<usize>,
        extend_beyond_range: bool,
        tsm_uncertainty_range_extension_mode: TSMUncertaintyRangeExtensionMode,
    ) -> AlignmentResult<U64Cost> {
        let data = QueryData {
            reference_name,
            reference,
            query_name,
            query,
            range,
            additional_tsm_starts_and_ends,
            cost_limit,
            memory_limit,
            extend_beyond_range,
            tsm_uncertainty_range_extension_mode,
        };
        self.align_select_min_length_strategy(data)
    }

    fn align_select_min_length_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        data: QueryData<SubsequenceType>,
    ) -> AlignmentResult<U64Cost> {
        match self.min_length_strategy {
            MinLengthStrategySelector::None => self
                .align_select_chaining_strategy::<_, NoTemplateSwitchMinLengthStrategy<U64Cost>>(data),
            MinLengthStrategySelector::Lookahead => self.align_select_chaining_strategy::<_, LookaheadTemplateSwitchMinLengthStrategy<U64Cost>>(data),
            MinLengthStrategySelector::PreprocessPrice => self.align_select_chaining_strategy::<_, PreprocessedTemplateSwitchMinLengthStrategy<false, U64Cost>>(data),
            MinLengthStrategySelector::PreprocessFilter => self.align_select_chaining_strategy::<_, PreprocessedTemplateSwitchMinLengthStrategy<true, U64Cost>>(data),
            MinLengthStrategySelector::PreprocessLookahead => self.align_select_chaining_strategy::<_, PreprocessedLookaheadTemplateSwitchMinLengthStrategy<U64Cost>>(data),
        }
    }

    fn align_select_chaining_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
    >(
        &self,
        data: QueryData<SubsequenceType>,
    ) -> AlignmentResult<U64Cost> {
        match self.chaining_strategy {
            ChainingStrategySelector::None => {
                self.align_select_total_length_strategy::<_, ML, NoChainingStrategy<U64Cost>>(data)
            }
            #[cfg(test)]
            ChainingStrategySelector::LowerBound => self
                .align_select_total_length_strategy::<_, ML, LowerBoundChainingStrategy<U64Cost>>(
                    data,
                ),
        }
    }

    fn align_select_total_length_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
    >(
        &self,
        data: QueryData<SubsequenceType>,
    ) -> AlignmentResult<U64Cost> {
        match self.total_length_strategy {
            TotalLengthStrategySelector::None => self
                .align_select_no_ts_strategy::<_, ML, CH, NoTemplateSwitchTotalLengthStrategy>(
                    data,
                ),
            TotalLengthStrategySelector::Maximise => self
                .align_select_no_ts_strategy::<_, ML, CH, MaxTemplateSwitchTotalLengthStrategy>(
                    data,
                ),
        }
    }

    fn align_select_no_ts_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
    >(
        &self,
        data: QueryData<SubsequenceType>,
    ) -> AlignmentResult<U64Cost> {
        if self.no_ts {
            self.align_select_descendant_strategy::<_, ML, CH, TL, MaxTemplateSwitchCountStrategy>(
                data, 0,
            )
        } else {
            self.align_select_descendant_strategy::<_, ML, CH, TL, NoTemplateSwitchCountStrategy>(
                data,
                (),
            )
        }
    }

    fn align_select_descendant_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
    >(
        &self,
        data: QueryData<SubsequenceType>,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<U64Cost> {
        match self.descendant_strategy {
            DescendantStrategySelector::AllowAny => self
                .align_select_primary_range_strategy::<_, ML, CH, TL, TC, AnyTemplateSwitchDescendantStrategy>(
                    data,
                    count_strategy_memory,
                ),
            DescendantStrategySelector::AllowOnlyAllEqual => self
                .align_select_primary_range_strategy::<_, ML, CH, TL, TC, OnlyEqualTemplateSwitchDescendantStrategy>(
                    data,
                    count_strategy_memory,
                ),
        }
    }

    fn align_select_primary_range_strategy<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
        DS: TemplateSwitchDescendantStrategy,
    >(
        &self,
        data: QueryData<SubsequenceType>,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<U64Cost> {
        match self.primary_range_strategy {
            PrimaryRangeStrategySelector::Prune => self
                .align_call::<_, ML, CH, TL, TC, DS, RangePrunePrimaryRangeStrategy>(
                    data,
                    count_strategy_memory,
                ),
            PrimaryRangeStrategySelector::NoPrune => self
                .align_call::<_, ML, CH, TL, TC, DS, NoPrunePrimaryRangeStrategy>(
                    data,
                    count_strategy_memory,
                ),
        }
    }

    fn align_call<
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
        DS: TemplateSwitchDescendantStrategy,
        PR: PrimaryRangeStrategy,
    >(
        &self,
        data: QueryData<SubsequenceType>,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<U64Cost> {
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
                PR,
                TL,
                DS,
            >,
            _,
        >(
            data.reference,
            data.query,
            data.reference_name,
            data.query_name,
            data.range,
            data.additional_tsm_starts_and_ends,
            &self.costs,
            DynamicStrategies {
                ts_14_out_of_range: self.ts_14_out_of_range_strategy,
            },
            data.cost_limit,
            data.memory_limit,
            self.force_label_correcting,
            data.extend_beyond_range,
            data.tsm_uncertainty_range_extension_mode,
            count_strategy_memory,
        )
    }
}
