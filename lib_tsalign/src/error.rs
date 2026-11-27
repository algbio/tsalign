use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("An IO error occurred: {0}.")]
    Io(#[from] std::io::Error),

    #[error("A parsing error of kind '{kind:?}' occurred when the remaining input was '{input}'.")]
    Parser {
        input: String,
        kind: nom::error::ErrorKind,
    },

    #[error("Parsing was unsuccessful due to incomplete input: {0:?}.")]
    ParserIncomplete(nom::Needed),

    #[error("The cost table name {0} was encountered twice.")]
    DuplicateCostTableName(String),

    #[error(
        "The template switch cost file contained a wrong set of cost table names. Expected: {expected:?}. Actual: {actual:?}."
    )]
    WrongCostTableNames {
        actual: Vec<String>,
        expected: Vec<String>,
    },

    #[error(
        "A cost function was attempted to create from a sequence whose index does not strictly increase at {index}."
    )]
    CostFunctionIndexNotIncreasing { index: usize },

    #[error(
        "The rq and qr offset costs are not V-shaped. They must be non-decreasing when starting from zero in both positive and negative direction."
    )]
    RQQROffsetCostsNotVShaped,

    #[error(
        "The rr and qq offset costs are not V-shaped. They must be non-decreasing when starting from zero in both positive and negative direction."
    )]
    RRQQOffsetCostsNotVShaped,

    #[error(
        "The template switch cost file sets non-zero flank lengths (left_flank_length = {left}, right_flank_length = {right}), but flanks are not supported anymore. Remove the '# Limits' section and the flank edit cost tables from the file."
    )]
    NonZeroFlankLength { left: isize, right: isize },

    #[error(
        "The length difference costs are not V-shaped. They must be non-decreasing when starting from zero in both positive and negative direction."
    )]
    LengthDifferenceCostsNotVShaped,
}

/// One of the two input sequences of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceRole {
    Reference,
    Query,
}

impl SequenceRole {
    /// The name of the sequence with a capital first letter, for the start of a sentence.
    pub fn capitalised(&self) -> &'static str {
        match self {
            Self::Reference => "Reference",
            Self::Query => "Query",
        }
    }
}

impl std::fmt::Display for SequenceRole {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Reference => write!(f, "reference"),
            Self::Query => write!(f, "query"),
        }
    }
}

/// A reason why [`Aligner::align`](crate::Aligner::align) or
/// [`Aligner::align_sequences`](crate::Aligner::align_sequences) fails:
/// the inputs are invalid, or the A* search failed ([`Self::Search`]).
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum AlignerInputError {
    /// A sequence contains a character that is not part of the alphabet of the aligner.
    #[error("{} contains non-alphabet character: {error}", .sequence.capitalised())]
    NonAlphabetCharacter {
        sequence: SequenceRole,
        error: compact_genome::interface::alphabet::AlphabetError,
    },

    /// The alignment range does not fit the sequences.
    #[error("Invalid alignment range: {0}")]
    InvalidRange(#[from] crate::alignment_geometry::InvalidAlignmentRange),

    /// The additional TSM starts and ends cannot be computed.
    #[error("{0}")]
    AdditionalExplicitTSMStartsAndEnds(
        #[from] crate::aligner::ts_14_out_of_range::AdditionalExplicitTSMStartsAndEndsError,
    ),

    /// TSMs that start or end outside of the alignment range
    /// ([`Ts14OutOfRangeStrategy::Allow`](crate::aligner::ts_14_out_of_range::Ts14OutOfRangeStrategy::Allow))
    /// are not supported by the aligner yet.
    #[error(
        "Allowing template switches to start or end outside of the alignment range is not supported yet"
    )]
    Ts14OutOfRangeNotSupported,

    /// Additional TSM starts and ends were given,
    /// but the [`Ts14OutOfRangeStrategy`](crate::aligner::ts_14_out_of_range::Ts14OutOfRangeStrategy)
    /// does not allow TSMs to start or end outside of the alignment range, so they would be ignored.
    #[error(
        "Additional explicit TSM starts and ends were given, but they are ignored, \
         because TSMs are not allowed to start or end outside of the alignment range"
    )]
    IgnoredAdditionalExplicitTSMStartsAndEnds,

    /// The A* search failed.
    ///
    /// This happens if the search needs more nodes at the same time than the arena index type allows
    /// ([`AStarError::ArenaFull`](generic_a_star::AStarError::ArenaFull)).
    /// Then the alignment can be repeated with a larger index type, see
    /// [`Aligner::set_arena_index`](crate::Aligner::set_arena_index).
    #[error(transparent)]
    Search(#[from] generic_a_star::AStarError),
}

/// Convert an ASCII sequence into a sequence of the given alphabet.
///
/// Returns [`AlignerInputError::NonAlphabetCharacter`] if the sequence contains a character that is not part of the alphabet.
pub fn sequence_from_ascii<AlphabetType: compact_genome::interface::alphabet::Alphabet>(
    sequence: SequenceRole,
    ascii: &[u8],
) -> std::result::Result<
    compact_genome::implementation::vec_sequence::VectorGenome<AlphabetType>,
    AlignerInputError,
> {
    use compact_genome::interface::sequence::OwnedGenomeSequence;

    compact_genome::implementation::vec_sequence::VectorGenome::from_slice_u8(ascii)
        .map_err(|error| AlignerInputError::NonAlphabetCharacter { sequence, error })
}
