use clap::ValueEnum;

pub use crate::anchors::inexact_pruning::{
    InexactPruningParameters, InexactPruningRule, InexactPruningRules,
};

/// Performance parameters for chainalign.
///
/// Use the `default()` method for a reasonable choice of parameters that works well in many cases.
pub struct AlignmentPerformanceParameters<Cost> {
    /// The step width for generating successors during chaining.
    ///
    /// At most `max_successors` will be generated at a time, but at least all with minimum chaining cost.
    pub max_successors: usize,

    /// The cost until which the cost function is initialised exactly.
    ///
    /// Anchor chainings with a higher exact cost are initialised based on the lower bound.
    pub max_exact_cost_function_cost: Cost,

    /// The closed list type to use for chaining.
    pub closed_list: ChainingClosedList,

    /// The open list type to use for chaining.
    pub open_list: ChainingOpenList,

    /// Which inexact anchors to remove before chaining.
    ///
    /// See [`AnchorPruning`].
    pub anchor_pruning: AnchorPruning,

    /// The rules and parameters of [`AnchorPruning::Inexact`].
    ///
    /// Ignored for the other modes.
    pub inexact_pruning: InexactPruningParameters,
}

/// The mode of removing inexact anchors before chaining.
///
/// Anchor pruning only applies to inexact anchors (`max_anchor_mutations >= 1`).
/// Exact anchors are never pruned.
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq, ValueEnum)]
pub enum AnchorPruning {
    /// Keep all anchors.
    ///
    /// This is useful to compare against the chaining with all anchors.
    None,
    /// Remove dominated inexact primary anchors, but only where this provably does not change the cost
    /// of the optimal chain.
    ///
    /// See the module `anchors::domination` for the rule and its proof.
    /// It only takes effect with `max_anchor_mutations == 1`, the max-match-run lower bound and
    /// primary costs with `gap_extend <= gap_open`.
    #[default]
    Exact,
    /// Remove primary and secondary inexact anchors with heuristic rules, for any `max_anchor_mutations >= 1`.
    ///
    /// This may increase the cost of the optimal chain (it never decreases it), and may lose template switches.
    /// The rules can be switched on and off individually, see [`InexactPruningRule`].
    Inexact,
}

/// The closed list type to use for chaining.
#[derive(Debug, Clone, ValueEnum)]
pub enum ChainingClosedList {
    /// Use [`HashMap`](std::collections::HashMap) as closed list with [`FxHasher`](rustc_hash::FxHasher) as hasher.
    FxHashMap,
    /// Use a special-purpose closed list.
    Special,
}

/// The open list type to use for chaining.
#[derive(Debug, Clone, ValueEnum)]
pub enum ChainingOpenList {
    /// Use [`BinaryHeap`](std::collections::BinaryHeap) as open list.
    StdHeap,
    /// Use [`LinearHeap`](generic_a_star::open_lists::linear_heap::LinearHeap) as open list.
    LinearHeap,
}

impl<Cost: From<u8>> Default for AlignmentPerformanceParameters<Cost> {
    fn default() -> Self {
        Self {
            max_successors: 1,
            max_exact_cost_function_cost: 1u8.into(),
            closed_list: ChainingClosedList::Special,
            open_list: ChainingOpenList::LinearHeap,
            anchor_pruning: AnchorPruning::default(),
            inexact_pruning: InexactPruningParameters::default(),
        }
    }
}
