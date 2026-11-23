use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::{AStarNode, cost::AStarCost};

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier, Node};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait NodeOrdStrategy<Cost, PrimaryMatch: PrimaryMatchStrategy<Cost>>:
    AlignmentStrategy
{
    fn cmp<Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct CostOnlyNodeOrdStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AntiDiagonalNodeOrdStrategy;

impl<Cost: AStarCost, PrimaryMatch: PrimaryMatchStrategy<Cost>> NodeOrdStrategy<Cost, PrimaryMatch>
    for CostOnlyNodeOrdStrategy
{
    fn cmp<Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering {
        // Nodes with the lower-bound costs must be ordered by their actual cost.
        // This is because nodes with the same identifier but higher cost and lower lower bound result in worse solutions.
        n1.node_data
            .lower_bound_cost()
            .cmp(&n2.node_data.lower_bound_cost())
            .then_with(|| n1.cost().cmp(&n2.cost()))
    }
}

impl<Cost: AStarCost, PrimaryMatch: PrimaryMatchStrategy<Cost>> NodeOrdStrategy<Cost, PrimaryMatch>
    for AntiDiagonalNodeOrdStrategy
{
    fn cmp<Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering {
        // Nodes with the lower-bound costs must be ordered by their actual cost.
        // This is because nodes with the same identifier but higher cost and lower lower bound result in worse solutions.
        match n1
            .node_data
            .lower_bound_cost()
            .cmp(&n2.node_data.lower_bound_cost())
            .then_with(|| n1.cost().cmp(&n2.cost()))
        {
            // This secondary ordering may make things actually slower.
            // While it does reduce the number of visited nodes a little bit,
            // it also makes heap operations more expensive.
            // Preliminary testing showed that this would be a slowdown.
            std::cmp::Ordering::Equal => n2
                .node_data
                .identifier
                .anti_diagonal()
                .cmp(&n1.node_data.identifier.anti_diagonal()),
            ordering => ordering,
        }
    }
}

impl AlignmentStrategy for CostOnlyNodeOrdStrategy {
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

impl AlignmentStrategy for AntiDiagonalNodeOrdStrategy {
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
