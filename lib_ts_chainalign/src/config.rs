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
    /// Compute a tight lower bound.
    Tight,
}
