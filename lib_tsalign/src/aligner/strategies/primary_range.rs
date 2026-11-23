use std::ops::Range;

use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait PrimaryRangeStrategy: AlignmentStrategy {
    fn reference_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize>;

    fn query_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize>;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoPrunePrimaryRangeStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct RangePrunePrimaryRangeStrategy;

impl PrimaryRangeStrategy for NoPrunePrimaryRangeStrategy {
    fn reference_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize> {
        0..context.reference.len()
    }

    fn query_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize> {
        0..context.query.len()
    }
}

impl PrimaryRangeStrategy for RangePrunePrimaryRangeStrategy {
    fn reference_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize> {
        context.range.reference_range()
    }

    fn query_range<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Range<usize> {
        context.range.query_range()
    }
}

impl AlignmentStrategy for NoPrunePrimaryRangeStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as super::primary_match::PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }
}

impl AlignmentStrategy for RangePrunePrimaryRangeStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as super::primary_match::PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }
}
