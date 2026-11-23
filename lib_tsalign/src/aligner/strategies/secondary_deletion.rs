use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait SecondaryDeletionStrategy: AlignmentStrategy {
    fn allow_secondary_deletions() -> bool;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AllowSecondaryDeletionStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct ForbidSecondaryDeletionStrategy;

impl SecondaryDeletionStrategy for AllowSecondaryDeletionStrategy {
    fn allow_secondary_deletions() -> bool {
        true
    }
}

impl SecondaryDeletionStrategy for ForbidSecondaryDeletionStrategy {
    fn allow_secondary_deletions() -> bool {
        false
    }
}

impl AlignmentStrategy for AllowSecondaryDeletionStrategy {
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
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl AlignmentStrategy for ForbidSecondaryDeletionStrategy {
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
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
