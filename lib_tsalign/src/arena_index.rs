//! The choice of the integer type of the node indices of the A* search.

use std::fmt::Display;

/// The integer type of the node indices of the A* search of an aligner.
///
/// The nodes of the search are stored in an arena and referred to by their index
/// (see [`generic_a_star::arena`]), and the index type limits how many nodes the search can store at the same time.
///
/// * [`U32`](Self::U32) (the default) allows at most `u32::MAX` = 4 294 967 295 nodes.
///   This is enough for searches that fit into a few hundred GB of memory.
/// * [`U64`](Self::U64) has no practical limit,
///   but the open and closed lists need more memory (about 3 % more in total for the classic aligner).
///
/// If the arena of a search is full, then the aligner returns an error
/// ([`AStarError::ArenaFull`](generic_a_star::AStarError::ArenaFull)) instead of an alignment,
/// and the alignment can be repeated with [`U64`](Self::U64).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum ArenaIndexType {
    /// 32-bit node indices: at most 4 294 967 295 nodes at the same time.
    #[default]
    U32,
    /// 64-bit node indices: no practical limit, but more memory per node.
    U64,
    /// 16-bit node indices: at most 65 535 nodes, to test the handling of a full arena.
    #[cfg(test)]
    U16,
}

impl ArenaIndexType {
    /// The name of the index type, as in the CLI and the Python bindings, e.g. `"u32"`.
    pub fn name(&self) -> &'static str {
        match self {
            Self::U32 => "u32",
            Self::U64 => "u64",
            #[cfg(test)]
            Self::U16 => "u16",
        }
    }
}

impl Display for ArenaIndexType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}
