use std::marker::PhantomData;

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStarContext, AStarNode, cost::AStarCost};

use crate::{
    a_star_aligner::template_switch_distance::{
        AlignmentType, Context, Identifier,
        lower_bounds::template_switch::TemplateSwitchLowerBoundMatrix,
    },
    config::TemplateSwitchConfig,
};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait ShortcutStrategy<Cost>: AlignmentStrategy {
    type Memory;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) -> Self::Memory;

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, Shortcut = Self>,
    >(
        node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        context: &Context<SubsequenceType, Strategies>,
        opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    );
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoShortcutStrategy<Cost> {
    phantom_data: PhantomData<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct TemplateSwitchLowerBoundShortcutStrategy<Cost> {
    phantom_data: PhantomData<Cost>,
}

impl<Cost: AStarCost> ShortcutStrategy<Cost> for NoShortcutStrategy<Cost> {
    type Memory = ();

    fn initialise_memory<AlphabetType: Alphabet>(
        _config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) -> Self::Memory {
        // Do nothing.
    }

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, Shortcut = Self>,
    >(
        _node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        _context: &Context<SubsequenceType, Strategies>,
        _opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    ) {
        // Do nothing.
    }
}

impl<Cost: AStarCost> ShortcutStrategy<Cost> for TemplateSwitchLowerBoundShortcutStrategy<Cost> {
    type Memory = TemplateSwitchLowerBoundMatrix<Cost>;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) -> Self::Memory {
        TemplateSwitchLowerBoundMatrix::new(config)
    }

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, Shortcut = Self>,
    >(
        node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        context: &Context<SubsequenceType, Strategies>,
        opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    ) {
        match *node.identifier() {
            Identifier::Primary { flank_index, .. }
            | Identifier::PrimaryReentry { flank_index, .. }
                if flank_index == context.config.left_flank_length =>
            {
                opened_nodes_output.extend(
                    context
                        .memory
                        .shortcut
                        .iter()
                        .flat_map(|entry| {
                            node.generate_template_switch_shortcut_successor(
                                entry.x(),
                                entry.y(),
                                entry.cost(),
                                context,
                            )
                        })
                        .map(Into::into),
                );
            }

            _ => { /* Do nothing. */ }
        }
    }
}

impl<Cost: AStarCost> AlignmentStrategy for NoShortcutStrategy<Cost> {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            phantom_data: PhantomData,
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

impl<Cost: AStarCost> AlignmentStrategy for TemplateSwitchLowerBoundShortcutStrategy<Cost> {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            phantom_data: PhantomData,
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
