use serde::{Deserialize, Serialize};

/// Holds the parameters for the chaining algorithm.
#[derive(Clone, Serialize, Deserialize)]
pub struct ChainingLowerBoundConfig {
    pub anchor_k: u32,
    pub max_anchor_mutations: u8,
    pub inexact_lower_bound: InexactLowerBoundKind,
}

/// Chooses the kind of lower bound to use for chaining with inexact anchors.
#[derive(Clone, Serialize, Deserialize)]
pub enum InexactLowerBoundKind {
    /// Compute a lower bound assuming that in every window of size `anchor_k` there have to be at least `max_anchor_mutations + 1` mutations.
    WindowedMinMutations,
    /// Compute a lower bound assuming that there is no match run longer than `anchor_k - 1 - max_anchor_mutations`,
    /// except for alignments that are shorter than `anchor_k` on both sequences.
    MaxMatchRun,
}
