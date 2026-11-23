use clap::ValueEnum;

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

    /// Remove dominated inexact primary anchors before chaining.
    ///
    /// This does not change the cost of the optimal chain, see the module `anchors::domination`.
    /// It is only applied in the cases in which this is proven.
    /// Disabling it is useful to compare against the chaining with all anchors.
    pub prune_dominated_anchors: bool,
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
            prune_dominated_anchors: true,
        }
    }
}
