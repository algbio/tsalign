use generic_a_star::cost::U32Cost;

use crate::{
    alignment::GapType,
    alignment_history::{
        extension_graph::{HistoryExtensionGraph, HistoryExtensionGraphNodeIndex},
        history_alignment_operations::AlignmentHistoryOperation::{self, *},
        history_vec::UnsignedIntAlignmentHistoryVec,
        pareto_front::{ParetoFront, ParetoLabel},
    },
    costs::GapAffineCosts,
};

type Graph = HistoryExtensionGraph<UnsignedIntAlignmentHistoryVec<u64>>;

const ANCHOR_K: u8 = 5;
const MAX_ANCHOR_MUTATIONS: u8 = 1;

fn costs() -> GapAffineCosts<U32Cost> {
    GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8))
}

/// Returns the history node reached by extending the empty history with the given operations (oldest first).
fn node(
    graph: &mut Graph,
    operations: &[AlignmentHistoryOperation],
    anchor_k: u8,
    max_anchor_mutations: u8,
) -> HistoryExtensionGraphNodeIndex {
    operations
        .iter()
        .fold(graph.empty_node_id(), |node, operation| {
            graph
                .try_extend(node, *operation, anchor_k, max_anchor_mutations)
                .unwrap()
        })
}

fn label(
    graph: &Graph,
    cost: u32,
    gap_type: GapType,
    history: HistoryExtensionGraphNodeIndex,
) -> ParetoLabel<U32Cost> {
    ParetoLabel::new(U32Cost::from(cost), gap_type, history, graph)
}

#[test]
fn identical_labels_dominate_each_other() {
    let mut graph = Graph::new();
    let history = node(
        &mut graph,
        &[Match, Substitution],
        ANCHOR_K,
        MAX_ANCHOR_MUTATIONS,
    );
    let a = label(&graph, 4, GapType::None, history);
    assert!(a.dominates(&a, &costs()));
}

#[test]
fn cheaper_label_dominates() {
    let mut graph = Graph::new();
    let history = node(&mut graph, &[Match, GapInA], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let cheap = label(&graph, 3, GapType::InA, history);
    let expensive = label(&graph, 4, GapType::InA, history);
    assert!(cheap.dominates(&expensive, &costs()));
    assert!(!expensive.dominates(&cheap, &costs()));
}

#[test]
fn freer_history_dominates() {
    let mut graph = Graph::new();
    // Profiles: free = [(0, 0), (1, 1)], restricted = [(1, 1), (2, 2)].
    let free = node(&mut graph, &[Substitution], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let restricted = node(
        &mut graph,
        &[Substitution, Match],
        ANCHOR_K,
        MAX_ANCHOR_MUTATIONS,
    );

    let free_label = label(&graph, 2, GapType::None, free);
    let restricted_label = label(&graph, 2, GapType::None, restricted);
    assert!(free_label.dominates(&restricted_label, &costs()));
    assert!(!restricted_label.dominates(&free_label, &costs()));

    // Freedom does not make up for higher cost.
    let expensive_free_label = label(&graph, 3, GapType::None, free);
    assert!(!expensive_free_label.dominates(&restricted_label, &costs()));
    assert!(!restricted_label.dominates(&expensive_free_label, &costs()));
}

#[test]
fn incomparable_histories() {
    let mut graph = Graph::new();
    // Profiles: gap_in_a = [(0, 0), (0, 1)], gap_in_b = [(0, 0), (1, 0)].
    let gap_in_a = node(&mut graph, &[GapInA], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let gap_in_b = node(&mut graph, &[GapInB], ANCHOR_K, MAX_ANCHOR_MUTATIONS);

    for (gap_type_a, gap_type_b) in [(GapType::InA, GapType::InB), (GapType::None, GapType::None)] {
        let a = label(&graph, 3, gap_type_a, gap_in_a);
        let b = label(&graph, 3, gap_type_b, gap_in_b);
        assert!(!a.dominates(&b, &costs()));
        assert!(!b.dominates(&a, &costs()));
    }
}

#[test]
fn gap_type_slack() {
    let mut graph = Graph::new();
    let history = node(&mut graph, &[Match], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let costs = costs();

    // Ending in a gap is never worse than ending without a gap, since gap extension is cheaper than gap opening.
    let in_b = label(&graph, 5, GapType::InB, history);
    let none = label(&graph, 5, GapType::None, history);
    assert!(in_b.dominates(&none, &costs));
    assert!(!none.dominates(&in_b, &costs));

    // Ending without a gap dominates if it is cheaper by at least `gap_open - gap_extend`.
    assert!(label(&graph, 3, GapType::None, history).dominates(&in_b, &costs));
    assert!(!label(&graph, 4, GapType::None, history).dominates(&in_b, &costs));

    // Different gap types are incomparable unless the cost difference is at least `gap_open - gap_extend`.
    let in_a = label(&graph, 5, GapType::InA, history);
    assert!(!in_a.dominates(&in_b, &costs));
    assert!(!in_b.dominates(&in_a, &costs));
    assert!(label(&graph, 3, GapType::InB, history).dominates(&in_a, &costs));
    assert!(!label(&graph, 4, GapType::InB, history).dominates(&in_a, &costs));
}

#[test]
fn gap_type_slack_with_expensive_gap_extension() {
    let mut graph = Graph::new();
    let history = node(&mut graph, &[Match], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    // Unusual costs where extending a gap is more expensive than opening one.
    let costs = GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(1u8), U32Cost::from(3u8));

    let in_b = label(&graph, 5, GapType::InB, history);
    let none = label(&graph, 5, GapType::None, history);
    assert!(!in_b.dominates(&none, &costs));
    assert!(none.dominates(&in_b, &costs));
    assert!(label(&graph, 3, GapType::InB, history).dominates(&none, &costs));
    assert!(!label(&graph, 4, GapType::InB, history).dominates(&none, &costs));
}

#[test]
fn front_insertion() {
    let mut graph = Graph::new();
    let free = node(&mut graph, &[Substitution], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let restricted = node(
        &mut graph,
        &[Substitution, Match],
        ANCHOR_K,
        MAX_ANCHOR_MUTATIONS,
    );
    let gap_in_a = node(&mut graph, &[GapInA], ANCHOR_K, MAX_ANCHOR_MUTATIONS);
    let costs = costs();

    let mut front = ParetoFront::new();
    assert!(front.is_empty());

    let restricted_label = label(&graph, 2, GapType::None, restricted);
    assert!(front.insert(restricted_label, &costs));
    assert_eq!(front.len(), 1);

    // An identical label is dominated.
    assert!(front.is_dominated(&restricted_label, &costs));
    assert!(!front.insert(restricted_label, &costs));
    assert_eq!(front.len(), 1);

    // A dominating label replaces the dominated label.
    let free_label = label(&graph, 2, GapType::None, free);
    assert!(!front.is_dominated(&free_label, &costs));
    assert!(front.insert(free_label, &costs));
    assert_eq!(front.iter().copied().collect::<Vec<_>>(), vec![free_label]);
    assert!(!front.insert(restricted_label, &costs));

    // An incomparable label is added.
    let gap_label = label(&graph, 3, GapType::InA, gap_in_a);
    assert!(front.insert(gap_label, &costs));
    assert_eq!(front.len(), 2);

    // A label that dominates both replaces both.
    let empty_label = label(&graph, 1, GapType::None, graph.empty_node_id());
    assert!(front.insert(empty_label, &costs));
    assert_eq!(front.iter().copied().collect::<Vec<_>>(), vec![empty_label]);
}

/// If profiles are not supported for the parameters, then only identical histories dominate each other.
#[test]
fn unsupported_profiles_fall_back_to_identity() {
    let (anchor_k, max_anchor_mutations) = (20, 8);
    let mut graph = Graph::new();
    let free = node(&mut graph, &[Substitution], anchor_k, max_anchor_mutations);
    let restricted = node(
        &mut graph,
        &[Substitution, Match],
        anchor_k,
        max_anchor_mutations,
    );
    assert!(graph.profile(free).is_none());
    assert!(graph.profile(restricted).is_none());
    let costs = costs();

    let free_label = label(&graph, 2, GapType::None, free);
    let restricted_label = label(&graph, 2, GapType::None, restricted);
    assert!(free_label.dominates(&free_label, &costs));
    assert!(!free_label.dominates(&restricted_label, &costs));
    assert!(!restricted_label.dominates(&free_label, &costs));
    assert!(free_label.dominates(&label(&graph, 3, GapType::None, free), &costs));

    // The empty history still has a profile, but it is incomparable to histories without profile.
    let empty_label = label(&graph, 0, GapType::None, graph.empty_node_id());
    assert!(empty_label.dominates(&empty_label, &costs));
    assert!(!empty_label.dominates(&free_label, &costs));
    assert!(!free_label.dominates(&empty_label, &costs));
}
