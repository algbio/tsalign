//! The input alphabets of the CLI and the dispatch from the runtime selection to the alphabet type.

use std::fmt::Debug;

use clap::ValueEnum;
use compact_genome::interface::alphabet::Alphabet;

#[cfg(feature = "alphabet_dna")]
use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
#[cfg(feature = "alphabet_dna_n")]
use compact_genome::implementation::alphabets::dna_alphabet_or_n::DnaAlphabetOrN;
#[cfg(feature = "alphabet_dna_iupac")]
use compact_genome::implementation::alphabets::dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet;
#[cfg(feature = "alphabet_rna")]
use compact_genome::implementation::alphabets::rna_alphabet::RnaAlphabet;
#[cfg(feature = "alphabet_rna_n")]
use compact_genome::implementation::alphabets::rna_alphabet_or_n::RnaAlphabetOrN;
#[cfg(feature = "alphabet_rna_iupac")]
use compact_genome::implementation::alphabets::rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet;

/// The alphabets that can be selected on the command line.
///
/// Each alphabet is only available if its cargo feature is enabled.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum InputAlphabet {
    #[cfg(feature = "alphabet_dna")]
    Dna,
    #[cfg(feature = "alphabet_dna_n")]
    DnaN,
    #[cfg(feature = "alphabet_rna")]
    Rna,
    #[cfg(feature = "alphabet_rna_n")]
    RnaN,
    #[cfg(feature = "alphabet_dna_iupac")]
    DnaIupac,
    #[cfg(feature = "alphabet_rna_iupac")]
    RnaIupac,
}

/// A computation that is generic over the alphabet type.
pub trait WithAlphabet {
    type Output;

    fn run<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(self) -> Self::Output;
}

impl InputAlphabet {
    /// Run `computation` with the alphabet type selected by `self`.
    pub fn dispatch<Computation: WithAlphabet>(
        self,
        computation: Computation,
    ) -> Computation::Output {
        match self {
            #[cfg(feature = "alphabet_dna")]
            Self::Dna => computation.run::<DnaAlphabet>(),
            #[cfg(feature = "alphabet_dna_n")]
            Self::DnaN => computation.run::<DnaAlphabetOrN>(),
            #[cfg(feature = "alphabet_rna")]
            Self::Rna => computation.run::<RnaAlphabet>(),
            #[cfg(feature = "alphabet_rna_n")]
            Self::RnaN => computation.run::<RnaAlphabetOrN>(),
            #[cfg(feature = "alphabet_dna_iupac")]
            Self::DnaIupac => computation.run::<DnaIupacNucleicAcidAlphabet>(),
            #[cfg(feature = "alphabet_rna_iupac")]
            Self::RnaIupac => computation.run::<RnaIupacNucleicAcidAlphabet>(),
        }
    }
}
