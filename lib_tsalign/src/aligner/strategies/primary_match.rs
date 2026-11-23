use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::cost::AStarCost;
use num_traits::Bounded;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::AlignmentStrategySelector;

pub trait PrimaryMatchStrategy<Cost>: Eq + Clone + Debug + Display {
    type Memory;
    type IdentifierPrimaryExtraData: Eq + Copy + Debug + Ord + Hash;

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn create_root_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData;

    fn generate_successor_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData;

    /// If false, then this node cannot have a primary match outgoing edge outside of a flank.
    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    /// If false, then this node cannot have a primary match outgoing edge within a flank.
    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost;

    fn always_generate_substitution() -> bool;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AllowPrimaryMatchStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxConsecutivePrimaryMatchStrategy;

pub struct MaxConsecutivePrimaryMatchMemory<Cost> {
    pub max_consecutive_primary_matches: usize,
    pub root_available_primary_matches: usize,
    pub fake_substitution_cost: Cost,
}

impl<Cost: Bounded> PrimaryMatchStrategy<Cost> for AllowPrimaryMatchStrategy {
    type Memory = ();
    type IdentifierPrimaryExtraData = ();

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }

    fn create_root_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData {
    }

    fn generate_successor_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData {
    }

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a match.
        true
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a match.
        true
    }

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost {
        Cost::max_value()
    }

    fn always_generate_substitution() -> bool {
        false
    }
}

impl<Cost: AStarCost> PrimaryMatchStrategy<Cost> for MaxConsecutivePrimaryMatchStrategy {
    type Memory = MaxConsecutivePrimaryMatchMemory<Cost>;
    type IdentifierPrimaryExtraData = usize;

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }

    fn create_root_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData {
        context.memory.primary_match.root_available_primary_matches
    }

    fn generate_successor_identifier_primary_extra_data<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::IdentifierPrimaryExtraData {
        match (identifier, alignment_type) {
            (
                Identifier::Primary { data, .. } | Identifier::PrimaryReentry { data, .. },
                AlignmentType::PrimaryMatch,
            ) => {
                debug_assert!(data > 0);
                data - 1
            }
            _ => context.memory.primary_match.max_consecutive_primary_matches,
        }
    }

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        let (Identifier::Primary { data, .. } | Identifier::PrimaryReentry { data, .. }) =
            identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        debug_assert!(data < isize::MAX as usize);
        data > 0
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Identifier<Self::IdentifierPrimaryExtraData>,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a flank match.
        true
    }

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Cost = Cost, PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost {
        context.memory.primary_match.fake_substitution_cost
    }

    fn always_generate_substitution() -> bool {
        true
    }
}

impl Display for AllowPrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrimaryMatch allowed")
    }
}

impl Display for MaxConsecutivePrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrimaryMatch")
    }
}
