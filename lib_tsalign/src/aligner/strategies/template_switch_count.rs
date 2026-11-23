use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait TemplateSwitchCountStrategy: AlignmentStrategy {
    type Memory;

    /// Called when a template switch has ended.
    fn increment_count(&mut self);

    /// If false, then no further template switch can be started from this node.
    fn can_start_another_template_switch<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<TemplateSwitchCount = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoTemplateSwitchCountStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxTemplateSwitchCountStrategy {
    template_switch_count: usize,
}

impl TemplateSwitchCountStrategy for NoTemplateSwitchCountStrategy {
    type Memory = ();

    fn increment_count(&mut self) {}

    fn can_start_another_template_switch<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<TemplateSwitchCount = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        true
    }
}

impl TemplateSwitchCountStrategy for MaxTemplateSwitchCountStrategy {
    type Memory = usize;

    fn increment_count(&mut self) {
        self.template_switch_count += 1;
    }

    fn can_start_another_template_switch<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<TemplateSwitchCount = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        self.template_switch_count < context.memory.template_switch_count
    }
}

impl AlignmentStrategy for NoTemplateSwitchCountStrategy {
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

impl AlignmentStrategy for MaxTemplateSwitchCountStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            template_switch_count: 0,
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
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
