//! Pareto fronts for multi-criteria searches over alignments with alignment histories.
//!
//! A search label consists of the cost of an alignment, its gap type, and its alignment history.
//! Label `A` dominates label `B` at the same alignment coordinates if
//!
//! * the history of `A` gives at least as much freedom to future alignments as the history of `B` (see [`HistoryProfile`]), and
//! * for every continuation `C`, the cost of `AC` is at most the cost of `BC`.
//!
//! The gap type only influences the cost of the first operation of `C`, and only if this operation is a gap.
//! Hence, the second condition holds if and only if the cost of `A` is at most that of `B`,
//! and the cost of `A` plus the cost of a gap in a or b is at most the cost of `B` plus the cost of the same gap.
//!
//! If `A` dominates `B`, then `B` can be discarded without changing the cost of any optimal alignment.
//! Since the history of `A` is at least as free as that of `B`, the same holds recursively for all extensions of `A` and `B`.

use generic_a_star::cost::AStarCost;

use crate::{
    alignment::GapType,
    alignment_history::{
        extension_graph::{HistoryExtensionGraph, HistoryExtensionGraphNodeIndex},
        history_vec::AlignmentHistory,
        profile::HistoryProfile,
    },
    costs::GapAffineCosts,
};

#[cfg(test)]
mod tests;

/// A search label that can be compared by domination.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParetoLabel<Cost> {
    cost: Cost,
    gap_type: GapType,
    history: HistoryKey,
}

/// The representation of an alignment history used for domination checks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HistoryKey {
    /// The profile of the history.
    Profile(HistoryProfile),
    /// The history itself, used if the profile cannot be computed for the given parameters.
    /// Then, a history only dominates itself.
    Node(HistoryExtensionGraphNodeIndex),
}

/// A set of labels at the same alignment coordinates such that no label dominates another.
#[derive(Debug, Clone)]
pub struct ParetoFront<Cost> {
    /// Most fronts contain a single label, so the first label is stored inline to avoid an allocation.
    ///
    /// If this is `None`, then `other_labels` is empty.
    first_label: Option<ParetoLabel<Cost>>,
    other_labels: Vec<ParetoLabel<Cost>>,
}

impl<Cost: AStarCost> ParetoLabel<Cost> {
    pub fn new<AlignmentHistoryVec: AlignmentHistory>(
        cost: Cost,
        gap_type: GapType,
        history: HistoryExtensionGraphNodeIndex,
        history_graph: &HistoryExtensionGraph<AlignmentHistoryVec>,
    ) -> Self {
        Self {
            cost,
            gap_type,
            history: if let Some(profile) = history_graph.profile(history) {
                HistoryKey::Profile(profile)
            } else {
                HistoryKey::Node(history)
            },
        }
    }

    pub fn cost(&self) -> Cost {
        self.cost
    }

    pub fn gap_type(&self) -> GapType {
        self.gap_type
    }

    /// Returns true if this label dominates the `other` label.
    ///
    /// Both labels must be at the same alignment coordinates, and `costs` must be the costs of the operations that follow.
    pub fn dominates(&self, other: &Self, costs: &GapAffineCosts<Cost>) -> bool {
        self.history_dominates(other) && self.cost_dominates(other, costs)
    }

    fn history_dominates(&self, other: &Self) -> bool {
        match (&self.history, &other.history) {
            (HistoryKey::Profile(profile), HistoryKey::Profile(other_profile)) => {
                profile.is_at_most(other_profile)
            }
            (HistoryKey::Node(node), HistoryKey::Node(other_node)) => node == other_node,
            _ => false,
        }
    }

    fn cost_dominates(&self, other: &Self, costs: &GapAffineCosts<Cost>) -> bool {
        self.cost <= other.cost
            && [GapType::InA, GapType::InB].into_iter().all(|gap_type| {
                self.cost + gap_cost(costs, self.gap_type, gap_type)
                    <= other.cost + gap_cost(costs, other.gap_type, gap_type)
            })
    }
}

/// The cost of a gap character of type `gap_type` that follows an alignment ending with `previous_gap_type`.
fn gap_cost<Cost: Copy>(
    costs: &GapAffineCosts<Cost>,
    previous_gap_type: GapType,
    gap_type: GapType,
) -> Cost {
    if previous_gap_type == gap_type {
        costs.gap_extend
    } else {
        costs.gap_open
    }
}

impl<Cost: AStarCost> ParetoFront<Cost> {
    pub fn new() -> Self {
        Self {
            first_label: None,
            other_labels: Vec::new(),
        }
    }

    /// Returns true if a label in this front dominates the given label.
    pub fn is_dominated(&self, label: &ParetoLabel<Cost>, costs: &GapAffineCosts<Cost>) -> bool {
        self.iter().any(|other| other.dominates(label, costs))
    }

    /// Inserts the given label into the front, unless it is dominated by a label in the front.
    /// Labels in the front that are dominated by the given label are removed.
    ///
    /// Returns true if the label was inserted.
    pub fn insert(&mut self, label: ParetoLabel<Cost>, costs: &GapAffineCosts<Cost>) -> bool {
        if self.is_dominated(&label, costs) {
            return false;
        }

        self.other_labels
            .retain(|other| !label.dominates(other, costs));
        if self
            .first_label
            .is_some_and(|first_label| label.dominates(&first_label, costs))
        {
            self.first_label = self.other_labels.pop();
        }

        if self.first_label.is_none() {
            self.first_label = Some(label);
        } else {
            self.other_labels.push(label);
        }
        true
    }

    pub fn len(&self) -> usize {
        usize::from(self.first_label.is_some()) + self.other_labels.len()
    }

    pub fn is_empty(&self) -> bool {
        self.first_label.is_none()
    }

    pub fn iter(&self) -> impl Iterator<Item = &ParetoLabel<Cost>> {
        self.first_label.iter().chain(&self.other_labels)
    }
}

impl<Cost: AStarCost> Default for ParetoFront<Cost> {
    fn default() -> Self {
        Self::new()
    }
}
