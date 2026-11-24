use std::{fmt::Debug, hash::Hash};

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

#[derive(Debug, Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SequencePair {
    pub reference_name: String,
    pub reference: String,
    pub reference_rc: String,
    pub query_name: String,
    pub query: String,
    pub query_rc: String,
}

impl SequencePair {
    pub fn new<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
    ) -> Self {
        Self {
            reference_name: reference_name.to_owned(),
            reference: reference.as_string(),
            reference_rc: reference
                .reverse_complement_iter()
                .map(Into::<char>::into)
                .collect(),
            query_name: query_name.to_owned(),
            query: query.as_string(),
            query_rc: query
                .reverse_complement_iter()
                .map(Into::<char>::into)
                .collect(),
        }
    }
}
