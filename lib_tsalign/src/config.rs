use compact_genome::{
    implementation::alphabets::dna_alphabet_or_n::DnaAlphabetOrN, interface::alphabet::Alphabet,
};
use generic_a_star::cost::AStarCost;
use num_traits::{Bounded, bounds::UpperBounded};

use crate::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    error::{Error, Result},
};

pub mod io;

#[derive(Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        deserialize = "Cost: serde::Deserialize<'de>", // omit AlphabetType
        serialize = "Cost: serde::Serialize" // omit AlphabetType
    ))
)]
pub struct TemplateSwitchConfig<AlphabetType, Cost> {
    /// The minimum length of a template switch is the first length for which
    /// the length costs are finite.
    pub template_switch_min_length: usize,

    // Base cost
    pub base_cost: BaseCost<Cost>,

    // Edit costs
    pub primary_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_forward_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_reverse_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,

    // Jump costs
    pub rq_qr_offset_costs: CostFunction<isize, Cost>,
    pub rr_qq_offset_costs: CostFunction<isize, Cost>,
    pub length_costs: CostFunction<usize, Cost>,
    pub length_difference_costs: CostFunction<isize, Cost>,
    pub forward_anti_descendant_gap_costs: CostFunction<isize, Cost>,
    pub reverse_anti_descendant_gap_costs: CostFunction<isize, Cost>,
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCost<Cost> {
    /// Descendant: reference; ancestor: reference; direction: forward.
    pub rrf: Cost,
    /// Descendant: reference; ancestor: query; direction: forward.
    pub rqf: Cost,
    /// Descendant: query; ancestor: reference; direction: forward.
    pub qrf: Cost,
    /// Descendant: query; ancestor: query; direction: forward.
    pub qqf: Cost,
    /// Descendant: reference; ancestor: reference; direction: reverse.
    pub rrr: Cost,
    /// Descendant: reference; ancestor: query; direction: reverse.
    pub rqr: Cost,
    /// Descendant: query; ancestor: reference; direction: reverse.
    pub qrr: Cost,
    /// Descendant: query; ancestor: query; direction: reverse.
    pub qqr: Cost,
}

impl<AlphabetType, Cost: Bounded + Ord> TemplateSwitchConfig<AlphabetType, Cost> {
    /// Returns an error if any cost function is malformed.
    pub fn verify(&self) -> Result<()> {
        if !self.rq_qr_offset_costs.is_v_shaped() {
            Err(Error::RQQROffsetCostsNotVShaped)
        } else if !self.rr_qq_offset_costs.is_v_shaped() {
            Err(Error::RRQQOffsetCostsNotVShaped)
        } else if !self.length_difference_costs.is_v_shaped() {
            Err(Error::LengthDifferenceCostsNotVShaped)
        } else {
            Ok(())
        }
    }
}

impl<AlphabetType, Cost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn secondary_edit_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &GapAffineAlignmentCostTable<AlphabetType, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.secondary_forward_edit_costs,
            TemplateSwitchDirection::Reverse => &self.secondary_reverse_edit_costs,
        }
    }

    pub fn anti_descendant_gap_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &CostFunction<isize, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.forward_anti_descendant_gap_costs,
            TemplateSwitchDirection::Reverse => &self.reverse_anti_descendant_gap_costs,
        }
    }

    pub fn offset_costs(
        &self,
        descendant: TemplateSwitchDescendant,
        ancestor: TemplateSwitchAncestor,
    ) -> &CostFunction<isize, Cost> {
        match (descendant, ancestor) {
            (TemplateSwitchDescendant::Reference, TemplateSwitchAncestor::Reference) => {
                &self.rr_qq_offset_costs
            }
            (TemplateSwitchDescendant::Reference, TemplateSwitchAncestor::Query) => {
                &self.rq_qr_offset_costs
            }
            (TemplateSwitchDescendant::Query, TemplateSwitchAncestor::Reference) => {
                &self.rq_qr_offset_costs
            }
            (TemplateSwitchDescendant::Query, TemplateSwitchAncestor::Query) => {
                &self.rr_qq_offset_costs
            }
        }
    }
}

impl<Cost: UpperBounded> BaseCost<Cost> {
    pub fn new_max() -> Self {
        Self {
            rrf: Cost::max_value(),
            rqf: Cost::max_value(),
            qrf: Cost::max_value(),
            qqf: Cost::max_value(),
            rrr: Cost::max_value(),
            rqr: Cost::max_value(),
            qrr: Cost::max_value(),
            qqr: Cost::max_value(),
        }
    }
}

impl<Cost: Clone> BaseCost<Cost> {
    pub fn get(
        &self,
        descendant: TemplateSwitchDescendant,
        ancestor: TemplateSwitchAncestor,
        direction: TemplateSwitchDirection,
    ) -> Cost {
        match (descendant, ancestor, direction) {
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Forward,
            ) => self.rrf.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Reverse,
            ) => self.rrr.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Forward,
            ) => self.rqf.clone(),
            (
                TemplateSwitchDescendant::Reference,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Reverse,
            ) => self.rqr.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Forward,
            ) => self.qrf.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Reference,
                TemplateSwitchDirection::Reverse,
            ) => self.qrr.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Forward,
            ) => self.qqf.clone(),
            (
                TemplateSwitchDescendant::Query,
                TemplateSwitchAncestor::Query,
                TemplateSwitchDirection::Reverse,
            ) => self.qqr.clone(),
        }
    }
}

impl<AlphabetType: Alphabet, Cost: Clone> Clone for TemplateSwitchConfig<AlphabetType, Cost> {
    fn clone(&self) -> Self {
        Self {
            template_switch_min_length: self.template_switch_min_length,
            base_cost: self.base_cost.clone(),
            primary_edit_costs: self.primary_edit_costs.clone(),
            secondary_forward_edit_costs: self.secondary_forward_edit_costs.clone(),
            secondary_reverse_edit_costs: self.secondary_reverse_edit_costs.clone(),
            rq_qr_offset_costs: self.rq_qr_offset_costs.clone(),
            rr_qq_offset_costs: self.rr_qq_offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
            forward_anti_descendant_gap_costs: self.forward_anti_descendant_gap_costs.clone(),
            reverse_anti_descendant_gap_costs: self.reverse_anti_descendant_gap_costs.clone(),
        }
    }
}

/// The default config in the plain `.tsa` format of [`TemplateSwitchConfig::read_plain`],
/// over the alphabet [`DnaAlphabetOrN`].
///
/// It is `lib_tsalign/default_config.tsa` of the repository, which equals `sample_tsa_config/config.tsa`.
/// [`TemplateSwitchConfig::default`] parses it.
pub const DEFAULT_CONFIG_TSA: &str = include_str!("../default_config.tsa");

/// Parses [`DEFAULT_CONFIG_TSA`] (a few microseconds).
///
/// The default config is given over [`DnaAlphabetOrN`].
/// For other alphabets, each character takes the costs of the character with the same ASCII representation,
/// except that `U` takes the costs of `T`, and characters that do not exist in [`DnaAlphabetOrN`]
/// (for example IUPAC ambiguity codes) take the costs of `N`.
impl<AlphabetType: Alphabet, Cost: AStarCost> Default for TemplateSwitchConfig<AlphabetType, Cost> {
    fn default() -> Self {
        let config: TemplateSwitchConfig<DnaAlphabetOrN, Cost> = DEFAULT_CONFIG_TSA
            .parse()
            .expect("the built-in default config is valid");
        config.translate_alphabet(|character| {
            let character = character.to_ascii_uppercase();
            if character == b'U' {
                b'T'
            } else if DnaAlphabetOrN::ascii_to_character(character).is_ok() {
                character
            } else {
                b'N'
            }
        })
    }
}

impl<SourceAlphabet: Alphabet, Cost: AStarCost> TemplateSwitchConfig<SourceAlphabet, Cost> {
    /// Convert this config into a config over `TargetAlphabet`.
    ///
    /// Each character of `TargetAlphabet` takes the costs of the character of `SourceAlphabet`
    /// whose ASCII representation is `translate` applied to its own ASCII representation.
    ///
    /// # Panics
    ///
    /// If `translate` returns an ASCII character that is not in `SourceAlphabet`.
    fn translate_alphabet<TargetAlphabet: Alphabet>(
        &self,
        translate: impl Fn(u8) -> u8,
    ) -> TemplateSwitchConfig<TargetAlphabet, Cost> {
        let translate_table = |table: &GapAffineAlignmentCostTable<SourceAlphabet, Cost>| {
            table
                .translate_alphabet(&translate)
                .expect("translate returns characters of the source alphabet")
        };

        TemplateSwitchConfig {
            template_switch_min_length: self.template_switch_min_length,
            base_cost: self.base_cost.clone(),
            primary_edit_costs: translate_table(&self.primary_edit_costs),
            secondary_forward_edit_costs: translate_table(&self.secondary_forward_edit_costs),
            secondary_reverse_edit_costs: translate_table(&self.secondary_reverse_edit_costs),
            rq_qr_offset_costs: self.rq_qr_offset_costs.clone(),
            rr_qq_offset_costs: self.rr_qq_offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
            forward_anti_descendant_gap_costs: self.forward_anti_descendant_gap_costs.clone(),
            reverse_anti_descendant_gap_costs: self.reverse_anti_descendant_gap_costs.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::alphabets::{
            dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
            dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet,
            rna_alphabet::RnaAlphabet,
        },
        interface::alphabet::Alphabet,
    };
    use generic_a_star::cost::{U32Cost, U64Cost};

    use super::{DEFAULT_CONFIG_TSA, TemplateSwitchConfig};

    #[test]
    fn default_config_equals_sample_config() {
        let sample = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../sample_tsa_config/config.tsa");
        assert_eq!(
            DEFAULT_CONFIG_TSA,
            std::fs::read_to_string(&sample).unwrap(),
            "the built-in default config differs from {}",
            sample.display()
        );
    }

    #[test]
    fn default_is_the_parsed_default_config() {
        let parsed: TemplateSwitchConfig<DnaAlphabetOrN, U64Cost> =
            DEFAULT_CONFIG_TSA.parse().unwrap();
        assert_eq!(TemplateSwitchConfig::default(), parsed);
        TemplateSwitchConfig::<DnaAlphabetOrN, U32Cost>::default()
            .verify()
            .unwrap();
    }

    /// Assert that `character` of `AlphabetType` has the default costs of `expected` of [`DnaAlphabetOrN`].
    fn assert_same_costs<AlphabetType: Alphabet>(character: u8, expected: u8) {
        let config = TemplateSwitchConfig::<AlphabetType, U64Cost>::default();
        let reference = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default();
        let character = AlphabetType::ascii_to_character(character).unwrap();
        let expected = DnaAlphabetOrN::ascii_to_character(expected).unwrap();
        let other = AlphabetType::ascii_to_character(b'A').unwrap();
        let reference_other = DnaAlphabetOrN::ascii_to_character(b'A').unwrap();

        for (table, reference_table) in [
            (&config.primary_edit_costs, &reference.primary_edit_costs),
            (
                &config.secondary_forward_edit_costs,
                &reference.secondary_forward_edit_costs,
            ),
            (
                &config.secondary_reverse_edit_costs,
                &reference.secondary_reverse_edit_costs,
            ),
        ] {
            assert_eq!(
                table.gap_open_cost(character.clone()),
                reference_table.gap_open_cost(expected)
            );
            assert_eq!(
                table.gap_extend_cost(character.clone()),
                reference_table.gap_extend_cost(expected)
            );
            assert_eq!(
                table.match_or_substitution_cost(character.clone(), other.clone()),
                reference_table.match_or_substitution_cost(expected, reference_other)
            );
            assert_eq!(
                table.match_or_substitution_cost(character.clone(), character.clone()),
                reference_table.match_or_substitution_cost(expected, expected)
            );
        }
        assert_eq!(config.base_cost, reference.base_cost);
        assert_eq!(config.length_costs, reference.length_costs);
    }

    #[test]
    fn default_for_other_alphabets() {
        for character in *b"ACGT" {
            assert_same_costs::<DnaAlphabet>(character, character);
            assert_same_costs::<DnaIupacNucleicAcidAlphabet>(character, character);
        }
        assert_same_costs::<DnaIupacNucleicAcidAlphabet>(b'N', b'N');
        // Ambiguity codes take the costs of `N`.
        assert_same_costs::<DnaIupacNucleicAcidAlphabet>(b'R', b'N');
        assert_same_costs::<RnaAlphabet>(b'U', b'T');
        assert_same_costs::<RnaAlphabet>(b'G', b'G');
    }
}
