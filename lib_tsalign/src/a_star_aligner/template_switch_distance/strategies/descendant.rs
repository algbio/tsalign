use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{
    AlignmentType, Context, Identifier, TemplateSwitchPrimary,
};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait TemplateSwitchDescendantStrategy: AlignmentStrategy {
    /// True if a TSM with this descendant can be started.
    fn is_descendant_allowed(&self, descendant: TemplateSwitchPrimary) -> bool;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AnyTemplateSwitchDescendantStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct OnlyEqualTemplateSwitchDescendantStrategy {
    allowed_descendant: Option<TemplateSwitchPrimary>,
}

impl TemplateSwitchDescendantStrategy for AnyTemplateSwitchDescendantStrategy {
    fn is_descendant_allowed(&self, _descendant: TemplateSwitchPrimary) -> bool {
        true
    }
}

impl TemplateSwitchDescendantStrategy for OnlyEqualTemplateSwitchDescendantStrategy {
    fn is_descendant_allowed(&self, descendant: TemplateSwitchPrimary) -> bool {
        if let Some(allowed_descendant) = self.allowed_descendant {
            allowed_descendant == descendant
        } else {
            true
        }
    }
}

impl AlignmentStrategy for AnyTemplateSwitchDescendantStrategy {
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

impl AlignmentStrategy for OnlyEqualTemplateSwitchDescendantStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            allowed_descendant: None,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        let mut successor = *self;
        if let Identifier::TemplateSwitchEntrance {
            template_switch_primary: descendant,
            ..
        } = identifier
        {
            if let Some(allowed_descendant) = successor.allowed_descendant {
                assert_eq!(allowed_descendant, descendant);
            }

            successor.allowed_descendant = Some(descendant);
        }

        successor
    }
}
