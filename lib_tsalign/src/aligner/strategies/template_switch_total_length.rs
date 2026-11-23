use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait TemplateSwitchTotalLengthStrategy: AlignmentStrategy {
    fn template_switch_total_length(&self) -> usize;

    /// True if this choice of strategy makes the search label correcting.
    fn makes_label_correcting() -> bool;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoTemplateSwitchTotalLengthStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxTemplateSwitchTotalLengthStrategy {
    template_switch_total_length: usize,
}

impl TemplateSwitchTotalLengthStrategy for NoTemplateSwitchTotalLengthStrategy {
    fn template_switch_total_length(&self) -> usize {
        0
    }

    fn makes_label_correcting() -> bool {
        false
    }
}

impl TemplateSwitchTotalLengthStrategy for MaxTemplateSwitchTotalLengthStrategy {
    fn template_switch_total_length(&self) -> usize {
        self.template_switch_total_length
    }

    fn makes_label_correcting() -> bool {
        true
    }
}

impl AlignmentStrategy for NoTemplateSwitchTotalLengthStrategy {
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

impl AlignmentStrategy for MaxTemplateSwitchTotalLengthStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            template_switch_total_length: 0,
        }
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
        alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        let increment = if matches!(
            alignment_type,
            AlignmentType::SecondaryMatch
                | AlignmentType::SecondarySubstitution
                | AlignmentType::SecondaryInsertion
        ) {
            1
        } else {
            0
        };

        Self {
            template_switch_total_length: self.template_switch_total_length + increment,
        }
    }
}
