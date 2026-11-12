use std::collections::HashMap;

use multicriteria_dijkstra::{
    MulticriteriaDijkstra, MulticriteriaDijkstraContext,
    node::{MulticriteriaDijkstraCost, MulticriteriaDijkstraIdentifier, MulticriteriaDijkstraNode},
    prune_list::{
        IndexedNode, MulticriteriaDijkstraPruneList, NodeWithIsTarget, NodeWithOptionalPredecessor,
        NodeWithOptionalPredecessorAndIsTarget,
    },
};

use super::{PruneList, TwoCriteriaCost};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Node {
    identifier: Identifier,
    cost: Cost,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Identifier {
    position_a: usize,
    position_b: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Cost {
    hops: usize,
    cost: usize,
}

impl MulticriteriaDijkstraNode for Node {
    type Identifier = Identifier;

    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> &Self::Cost {
        &self.cost
    }
}

impl MulticriteriaDijkstraIdentifier for Identifier {}

impl MulticriteriaDijkstraCost for Cost {
    fn dominates(&self, other: &Self) -> bool {
        self.hops <= other.hops && self.cost <= other.cost
    }
}

impl TwoCriteriaCost for Cost {
    type First = usize;
    type Second = usize;

    fn first(&self) -> Self::First {
        self.cost
    }

    fn second(&self) -> Self::Second {
        self.hops
    }
}

/// Edit distance alignment minimising both the number of operations and the cost.
struct Context {
    sequence_a: Vec<u8>,
    sequence_b: Vec<u8>,
}

impl Context {
    fn is_target_at(&self, position_a: usize, position_b: usize) -> bool {
        position_a == self.sequence_a.len() && position_b == self.sequence_b.len()
    }
}

impl MulticriteriaDijkstraContext for Context {
    type Node = Node;

    fn generate_root_nodes(&self, output: &mut impl Extend<NodeWithIsTarget<Self::Node>>) {
        output.extend(Some(NodeWithIsTarget {
            node: Node {
                identifier: Identifier {
                    position_a: 0,
                    position_b: 0,
                },
                cost: Cost { hops: 0, cost: 0 },
            },
            is_target: self.is_target_at(0, 0),
        }));
    }

    fn generate_successors<Index: Copy>(
        &self,
        IndexedNode { node, index }: IndexedNode<Self::Node, Index>,
        output: &mut impl Extend<NodeWithOptionalPredecessorAndIsTarget<Self::Node, Index>>,
    ) {
        let na = node.identifier.position_a;
        let nb = node.identifier.position_b;
        let mut successor = |position_a, position_b, cost| {
            output.extend(Some(NodeWithOptionalPredecessorAndIsTarget {
                node: Node {
                    identifier: Identifier {
                        position_a,
                        position_b,
                    },
                    cost: Cost {
                        hops: node.cost.hops + 1,
                        cost: node.cost.cost + cost,
                    },
                },
                predecessor_index: Some(index),
                is_target: self.is_target_at(position_a, position_b),
            }));
        };

        if na < self.sequence_a.len() && nb < self.sequence_b.len() {
            let cost = if self.sequence_a[na] == self.sequence_b[nb] {
                0
            } else {
                3
            };
            successor(na + 1, nb + 1, cost);
        }
        if na < self.sequence_a.len() {
            successor(na + 1, nb, 1);
        }
        if nb < self.sequence_b.len() {
            successor(na, nb + 1, 1);
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        self.is_target_at(node.identifier.position_a, node.identifier.position_b)
    }
}

fn count_path_operations(path: &[&Node]) -> (usize, usize, usize) {
    let mut masub = 0;
    let mut insertion = 0;
    let mut deletion = 0;
    for window in path.windows(2) {
        let from = &window[0].identifier;
        let to = &window[1].identifier;
        if from.position_a + 1 == to.position_a && from.position_b + 1 == to.position_b {
            masub += 1;
        } else if from.position_a + 1 == to.position_a && from.position_b == to.position_b {
            deletion += 1;
        } else if from.position_a == to.position_a && from.position_b + 1 == to.position_b {
            insertion += 1;
        } else {
            panic!("Invalid path");
        }
    }
    (masub, insertion, deletion)
}

fn check_alignment(
    sequence_a: &[u8],
    sequence_b: &[u8],
    mut expected_costs: Vec<Cost>,
    mut expected_paths: Vec<(usize, usize, usize)>,
) {
    let context = Context {
        sequence_a: sequence_a.to_vec(),
        sequence_b: sequence_b.to_vec(),
    };
    let mut dijkstra = MulticriteriaDijkstra::<PruneList<Node>, _>::new(context);
    dijkstra.search();

    let target_identifier = Identifier {
        position_a: sequence_a.len(),
        position_b: sequence_b.len(),
    };

    expected_costs.sort_unstable();
    let mut actual_costs = dijkstra
        .iter_costs(&target_identifier)
        .cloned()
        .collect::<Vec<_>>();
    actual_costs.sort_unstable();
    assert_eq!(actual_costs, expected_costs);

    expected_paths.sort_unstable();
    let mut actual_paths: Vec<_> = dijkstra
        .backtrack(&target_identifier)
        .unwrap()
        .into_iter()
        .map(|path| count_path_operations(&path))
        .collect();
    actual_paths.sort_unstable();
    assert_eq!(actual_paths, expected_paths);
}

#[test]
fn minimum_hop_alignment_1() {
    check_alignment(
        b"G",
        b"T",
        vec![Cost { hops: 1, cost: 3 }, Cost { hops: 2, cost: 2 }],
        vec![(1, 0, 0), (0, 1, 1)],
    );
}

#[test]
fn minimum_hop_alignment_4() {
    check_alignment(
        b"GGGG",
        b"TTTT",
        vec![
            Cost { hops: 4, cost: 12 },
            Cost { hops: 5, cost: 11 },
            Cost { hops: 6, cost: 10 },
            Cost { hops: 7, cost: 9 },
            Cost { hops: 8, cost: 8 },
        ],
        vec![(4, 0, 0), (3, 1, 1), (2, 2, 2), (1, 3, 3), (0, 4, 4)],
    );
}

#[test]
fn empty_alignment() {
    check_alignment(b"", b"", vec![Cost { hops: 0, cost: 0 }], vec![(0, 0, 0)]);
}

/// The naive prune list from the tests of the `multicriteria-dijkstra` crate.
#[derive(Debug, Default)]
struct NaivePruneList {
    active_nodes: HashMap<Identifier, Vec<usize>>,
    nodes: Vec<(Node, Option<usize>, bool, bool)>,
    target_nodes: Vec<usize>,
}

impl NaivePruneList {
    fn insert(&mut self, node: Node, is_target: bool) -> Option<usize> {
        for &index in self
            .active_nodes
            .get(&node.identifier)
            .unwrap_or(&Vec::new())
            .iter()
            .chain(self.target_nodes.iter())
        {
            if self.nodes[index].0.dominates(&node) {
                return None;
            }
        }

        let nodes = &mut self.nodes;
        let mut prune = |vec: &mut Vec<usize>| {
            vec.retain(|&index| {
                let result = !node.dominates(&nodes[index].0);
                nodes[index].3 |= !result;
                result
            })
        };
        if is_target {
            self.active_nodes.values_mut().for_each(&mut prune);
        } else if let Some(vec) = self.active_nodes.get_mut(&node.identifier) {
            prune(vec);
        }
        self.target_nodes.retain(|&index| !self.nodes[index].3);

        let index = self.nodes.len();
        self.active_nodes
            .entry(node.identifier.clone())
            .or_default()
            .push(index);
        self.nodes.push((node, None, true, false));
        if is_target {
            self.target_nodes.push(index);
        }
        Some(index)
    }

    fn costs_at(&self, identifier: &Identifier) -> Vec<Cost> {
        let mut result: Vec<_> = self
            .active_nodes
            .get(identifier)
            .into_iter()
            .flatten()
            .filter(|&&index| !self.nodes[index].3)
            .map(|&index| self.nodes[index].0.cost)
            .collect();
        result.sort_unstable();
        result
    }
}

/// A simple deterministic pseudo-random number generator.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self, bound: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) % bound as u64) as usize
    }
}

/// Insert random nodes into both prune lists and compare the results.
#[test]
fn random_insertions_match_naive() {
    for seed in 0..200 {
        let mut random = Lcg(seed);
        let mut prune_list = PruneList::<Node>::default();
        let mut naive = NaivePruneList::default();
        let identifier_count = 1 + random.next(7);
        let insertion_count = 1 + random.next(60);
        let cost_bound = 1 + random.next(12);

        for _ in 0..insertion_count {
            let node = Node {
                identifier: Identifier {
                    position_a: random.next(identifier_count),
                    position_b: 0,
                },
                cost: Cost {
                    hops: random.next(cost_bound),
                    cost: random.next(cost_bound),
                },
            };
            // Whether a node is a target must only depend on its identifier.
            let is_target = node.identifier.position_a % 3 == 0;

            let expected = naive.insert(node.clone(), is_target);
            let actual = prune_list.insert(NodeWithOptionalPredecessorAndIsTarget {
                node,
                predecessor_index: None,
                is_target,
            });
            assert_eq!(actual.is_some(), expected.is_some(), "seed {seed}");

            for position_a in 0..identifier_count {
                let identifier = Identifier {
                    position_a,
                    position_b: 0,
                };
                let mut actual_costs: Vec<_> = prune_list
                    .iter_costs_at_identifier(&identifier)
                    .copied()
                    .collect();
                actual_costs.sort_unstable();
                assert_eq!(actual_costs, naive.costs_at(&identifier), "seed {seed}");
            }
        }

        // All open nodes returned by `next` must be non-dominated.
        while let Some(IndexedNode { node, index }) = prune_list.next().map(|node| node.cloned()) {
            let NodeWithOptionalPredecessor { node: stored, .. } = prune_list.get(index).unwrap();
            assert_eq!(&node, stored);
            assert!(
                naive.costs_at(&node.identifier).contains(&node.cost),
                "seed {seed}"
            );
        }

        prune_list.clear();
        assert!(prune_list.next().is_none());
    }
}
