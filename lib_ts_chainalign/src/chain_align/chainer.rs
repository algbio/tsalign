use std::{fmt::Display, iter};

use generic_a_star::{AStarContext, AStarIdentifier, AStarNode, cost::AStarCost, reset::Reset};
use num_traits::Zero;

use crate::{
    alignment::ts_kind::TsKind,
    anchors::{Anchors, index::AnchorIndex},
    chain_align::chainer::max_successors_iterator::MaxSuccessorsIterator,
    chaining_cost_function::ChainingCostFunction,
    costs::TsLimits,
};

pub mod closed_list;
mod max_successors_iterator;

const DEBUG_CHAINER: bool = false;

pub struct Context<'anchors, 'chaining_cost_function, 'ts_limits, Cost> {
    pub anchors: &'anchors Anchors<Cost>,
    pub chaining_cost_function: &'chaining_cost_function mut ChainingCostFunction<Cost>,
    pub ts_limits: &'ts_limits TsLimits,
    pub k: usize,
    pub max_successors: usize,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
    offset_zero_cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Identifier {
    Start,
    StartToPrimary {
        offset: AnchorIndex,
    },
    StartToSecondary {
        ts_kind: TsKind,
        offset: AnchorIndex,
    },
    PrimaryToPrimary {
        index: AnchorIndex,
        offset: AnchorIndex,
    },
    PrimaryToSecondary {
        index: AnchorIndex,
        ts_kind: TsKind,
        offset: AnchorIndex,
    },
    SecondaryToSecondary {
        index: AnchorIndex,
        ts_kind: TsKind,
        /// The first secondary anchor that is part of the current template switch.
        ///
        /// Used to estimate the length of the resulting template switch.
        first_secondary_index: AnchorIndex,
        offset: AnchorIndex,
    },
    SecondaryToPrimary {
        index: AnchorIndex,
        ts_kind: TsKind,
        /// The first secondary anchor that is part of the current template switch.
        ///
        /// Used to estimate the length of the resulting template switch.
        first_secondary_index: AnchorIndex,
        offset: AnchorIndex,
    },
    End,
}

impl<'anchors, 'chaining_cost_function, 'ts_limits, Cost>
    Context<'anchors, 'chaining_cost_function, 'ts_limits, Cost>
{
    pub fn new(
        anchors: &'anchors Anchors<Cost>,
        chaining_cost_function: &'chaining_cost_function mut ChainingCostFunction<Cost>,
        ts_limits: &'ts_limits TsLimits,
        k: usize,
        max_successors: usize,
    ) -> Self {
        Self {
            anchors,
            chaining_cost_function,
            ts_limits,
            k,
            max_successors,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier::Start,
            predecessor: None,
            cost: Cost::zero(),
            offset_zero_cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let predecessor = Some(node.identifier.with_offset_zero());
        let predecessor_cost = node.cost;
        let offset_zero_cost = node.offset_zero_cost;
        let primary_end_anchor_index = self.anchors.primary_len();

        if DEBUG_CHAINER {
            println!("Generating successors of {node}");
        }

        match node.identifier {
            Identifier::Start => output.extend(
                iter::once(Identifier::StartToPrimary {
                    offset: AnchorIndex::zero(),
                })
                .chain(TsKind::iter().map(|ts_kind| Identifier::StartToSecondary {
                    ts_kind,
                    offset: AnchorIndex::zero(),
                }))
                .map(|identifier| Node {
                    identifier,
                    predecessor,
                    cost: predecessor_cost,
                    offset_zero_cost: predecessor_cost,
                })
                .chain(iter::once(Node {
                    identifier: Identifier::End,
                    predecessor,
                    cost: self.chaining_cost_function.start_to_end(),
                    offset_zero_cost: Cost::zero(),
                })),
            ),

            Identifier::StartToPrimary { offset } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_primary_from_start_in_cost_order(offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                if successor_index == primary_end_anchor_index {
                                    println!("Checking anchor end");
                                } else {
                                    println!(
                                        "Checking anchor P-{successor_index}: {}",
                                        self.anchors.primary(successor_index),
                                    );
                                }
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            Some(generate_primary_successors(
                                successor_index,
                                predecessor,
                                cost,
                                primary_end_anchor_index,
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::StartToPrimary {
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::StartToSecondary { ts_kind, offset } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_jump_12_from_start_in_cost_order(ts_kind, offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                println!(
                                    "Checking anchor S{}-{successor_index}: {}",
                                    ts_kind.digits(),
                                    self.anchors.secondary(successor_index, ts_kind),
                                );
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            Some(generate_secondary_successors(
                                successor_index,
                                ts_kind,
                                successor_index,
                                predecessor,
                                cost,
                                self.ts_limits.length_23.contains(&self.k),
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::StartToSecondary {
                            ts_kind,
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::PrimaryToPrimary { index, offset } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_primary_in_cost_order(index, offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                if successor_index == primary_end_anchor_index {
                                    println!("Checking anchor end");
                                } else {
                                    println!(
                                        "Checking anchor P-{successor_index}: {}",
                                        self.anchors.primary(successor_index),
                                    );
                                }
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            Some(generate_primary_successors(
                                successor_index,
                                predecessor,
                                cost,
                                primary_end_anchor_index,
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::PrimaryToPrimary {
                            index,
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::PrimaryToSecondary {
                index,
                ts_kind,
                offset,
            } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_jump_12_in_cost_order(index, ts_kind, offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                println!(
                                    "Checking anchor S{}-{successor_index}: {}",
                                    ts_kind.digits(),
                                    self.anchors.secondary(successor_index, ts_kind),
                                );
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            Some(generate_secondary_successors(
                                successor_index,
                                ts_kind,
                                successor_index,
                                predecessor,
                                cost,
                                self.ts_limits.length_23.contains(&self.k),
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::PrimaryToSecondary {
                            index,
                            ts_kind,
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::SecondaryToSecondary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_secondary_in_cost_order(index, ts_kind, offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                println!(
                                    "Checking anchor S{}-{successor_index}: {}",
                                    ts_kind.digits(),
                                    self.anchors.secondary(successor_index, ts_kind),
                                );
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            let first_anchor =
                                &self.anchors.secondary(first_secondary_index, ts_kind);
                            let ts_length = first_anchor.ts_length_until(
                                self.anchors.secondary(successor_index, ts_kind),
                                self.k,
                            );

                            Some(generate_secondary_successors(
                                successor_index,
                                ts_kind,
                                first_secondary_index,
                                predecessor,
                                cost,
                                self.ts_limits.length_23.contains(&ts_length),
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::SecondaryToSecondary {
                            index,
                            ts_kind,
                            first_secondary_index,
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::SecondaryToPrimary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => {
                let mut iter = MaxSuccessorsIterator::new(
                    self.chaining_cost_function
                        .iter_jump_34_in_cost_order(index, ts_kind, offset),
                    self.max_successors,
                );

                output.extend(
                    iter.by_ref()
                        .map_while(|(successor_index, chaining_cost)| {
                            if DEBUG_CHAINER {
                                if successor_index == primary_end_anchor_index {
                                    println!("Checking anchor end");
                                } else {
                                    println!(
                                        "Checking anchor P-{successor_index}: {}",
                                        self.anchors.primary(successor_index),
                                    );
                                }
                            }

                            let cost = offset_zero_cost.checked_add(&chaining_cost)?;
                            if DEBUG_CHAINER {
                                println!("Cost: {}+{}", offset_zero_cost, cost - offset_zero_cost);
                            }

                            Some(generate_primary_successors(
                                successor_index,
                                predecessor,
                                cost,
                                primary_end_anchor_index,
                            ))
                        })
                        .flatten(),
                );

                if let Some(next_cost) = iter.next_cost() {
                    output.extend(iter::once(Node {
                        identifier: Identifier::SecondaryToPrimary {
                            index,
                            ts_kind,
                            first_secondary_index,
                            offset: offset + AnchorIndex::from(iter.successor_count()),
                        },
                        // Skip nodes with non-zero offset in predecessor pointers.
                        predecessor: node.predecessor.map(Identifier::with_offset_zero),
                        cost: next_cost,
                        offset_zero_cost,
                    }));
                }
            }

            Identifier::End => { /* Has no successors */ }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier == Identifier::End
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }

    fn is_label_setting(&self) -> bool {
        false
    }
}

/// Generate successors for a primary node.
///
/// These are all successors with an offset of zero, i.e. they are actual successors rather than copies of the same node with different offset.
fn generate_primary_successors<Cost: Copy>(
    index: AnchorIndex,
    predecessor: Option<Identifier>,
    cost: Cost,
    primary_end_anchor_index: AnchorIndex,
) -> impl Iterator<Item = Node<Cost>> {
    (index == primary_end_anchor_index)
        .then_some(Identifier::End)
        .into_iter()
        .chain(
            (index != primary_end_anchor_index)
                .then(move || {
                    [Identifier::PrimaryToPrimary {
                        index,
                        offset: AnchorIndex::zero(),
                    }]
                    .into_iter()
                    .chain(
                        TsKind::iter().map(move |ts_kind| Identifier::PrimaryToSecondary {
                            index,
                            ts_kind,
                            offset: AnchorIndex::zero(),
                        }),
                    )
                })
                .into_iter()
                .flatten(),
        )
        .map(move |identifier| Node {
            identifier,
            predecessor,
            cost,
            offset_zero_cost: cost,
        })
}

/// Generate successors for a secondary node.
///
/// These are all successors with an offset of zero, i.e. they are actual successors rather than copies of the same node with different offset.
fn generate_secondary_successors<Cost: Copy>(
    index: AnchorIndex,
    ts_kind: TsKind,
    first_secondary_index: AnchorIndex,
    predecessor: Option<Identifier>,
    cost: Cost,
    can_34_jump: bool,
) -> impl Iterator<Item = Node<Cost>> {
    [
        can_34_jump.then(|| Identifier::SecondaryToPrimary {
            index,
            ts_kind,
            first_secondary_index,
            offset: AnchorIndex::zero(),
        }),
        Some(Identifier::SecondaryToSecondary {
            index,
            ts_kind,
            first_secondary_index,
            offset: AnchorIndex::zero(),
        }),
    ]
    .into_iter()
    .flatten()
    .map(move |identifier| Node {
        identifier,
        predecessor,
        cost,
        offset_zero_cost: cost,
    })
}

impl<Cost: AStarCost> Reset for Context<'_, '_, '_, Cost> {
    fn reset(&mut self) {
        // Nothing to do.
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = Identifier;

    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        Cost::zero()
    }

    fn secondary_maximisable_score(&self) -> usize {
        0
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        self.predecessor
    }
}

impl Identifier {
    pub fn with_offset_zero(mut self) -> Self {
        match &mut self {
            Identifier::Start | Identifier::End => { /* Nothing to do */ }
            Identifier::StartToPrimary { offset }
            | Identifier::StartToSecondary { offset, .. }
            | Identifier::PrimaryToPrimary { offset, .. }
            | Identifier::PrimaryToSecondary { offset, .. }
            | Identifier::SecondaryToSecondary { offset, .. }
            | Identifier::SecondaryToPrimary { offset, .. } => *offset = Zero::zero(),
        }
        self
    }

    pub fn offset(&self) -> AnchorIndex {
        match self {
            Identifier::Start | Identifier::End => Zero::zero(),
            Identifier::StartToPrimary { offset }
            | Identifier::StartToSecondary { offset, .. }
            | Identifier::PrimaryToPrimary { offset, .. }
            | Identifier::PrimaryToSecondary { offset, .. }
            | Identifier::SecondaryToSecondary { offset, .. }
            | Identifier::SecondaryToPrimary { offset, .. } => *offset,
        }
    }
}

impl AStarIdentifier for Identifier {}

impl<Cost: Ord> Ord for Node<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.cmp(&other.cost)
    }
}

impl<Cost: Ord> PartialOrd for Node<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Identifier::Start => write!(f, "start"),
            Identifier::StartToPrimary { offset } => write!(f, "start-to-primary-{offset}"),
            Identifier::StartToSecondary { ts_kind, offset } => {
                write!(f, "start-to-secondary{}-{offset}", ts_kind.digits())
            }
            Identifier::PrimaryToPrimary { index, offset } => {
                write!(f, "primary-to-primary-{index}-{offset}")
            }
            Identifier::PrimaryToSecondary {
                index,
                ts_kind,
                offset,
            } => write!(
                f,
                "primary-to-secondary{}-{index}-{offset}",
                ts_kind.digits()
            ),
            Identifier::SecondaryToSecondary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => write!(
                f,
                "secondary{}-to-secondary{}-{first_secondary_index}-{index}-{offset}",
                ts_kind.digits(),
                ts_kind.digits()
            ),
            Identifier::SecondaryToPrimary {
                index,
                ts_kind,
                first_secondary_index,
                offset,
            } => write!(
                f,
                "secondary{}-to-primary-{first_secondary_index}-{index}-{offset}",
                ts_kind.digits()
            ),
            Identifier::End => write!(f, "end"),
        }
    }
}
