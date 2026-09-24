use tagged_vec::TaggedVec;

use crate::chaining_lower_bounds::gap_affine::inexact_algo::{
    history_alignment_operations::AlignmentHistoryOperation, history_vec::AlignmentHistory,
};

pub struct HistoryGraph<AlignmentHistoryVec> {
    nodes: TaggedVec<HistoryNodeIndex, HistoryNode<AlignmentHistoryVec>>,
}

pub struct HistoryNode<AlignmentHistoryVec> {
    id: HistoryNodeIndex,
    history: AlignmentHistoryVec,

    /// The successor of this node if the next alignment operation is a match.
    match_successor: HistoryNodeIndex,

    /// The successor of this node if the next alignment operation is a substitution.
    substitution_successor: HistoryNodeIndex,

    /// The successor of this node if the next alignment operation is a gap_a.
    gap_a_successor: HistoryNodeIndex,

    /// The successor of this node if the next alignment operation is a gap_b.
    gap_b_successor: HistoryNodeIndex,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct HistoryNodeIndex(pub usize);

impl<AlignmentHistoryVec: AlignmentHistory> HistoryGraph<AlignmentHistoryVec> {
    pub fn new() -> Self {
        let mut nodes = TaggedVec::new();
        nodes.push(HistoryNode {
            id: HistoryNodeIndex(0),
            history: AlignmentHistoryVec::default(),
            match_successor: HistoryNodeIndex(usize::MAX),
            substitution_successor: HistoryNodeIndex(usize::MAX),
            gap_a_successor: HistoryNodeIndex(usize::MAX),
            gap_b_successor: HistoryNodeIndex(usize::MAX),
        });
        Self { nodes }
    }

    /// Returns true if the alignment history vector can be extended with another alignment history operation without becoming a valid anchor.
    pub fn can_extend(
        &mut self,
        node_id: HistoryNodeIndex,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> bool {
        let node = &mut self.nodes[node_id];
        let cache = match alignment {
            AlignmentHistoryOperation::Match => &mut node.match_successor,
            AlignmentHistoryOperation::Substitution => &mut node.substitution_successor,
            AlignmentHistoryOperation::GapA => &mut node.gap_a_successor,
            AlignmentHistoryOperation::GapB => &mut node.gap_b_successor,
        };

        if cache.is_dead_end() {
            false
        } else if cache.is_unknown() {
            let can_extend = node
                .history
                .can_extend(alignment, anchor_k, max_anchor_mutations);
            if !can_extend {
                *cache = HistoryNodeIndex::new_dead_end();
            }
            can_extend
        } else {
            true
        }
    }

    /// Return the extension of this the alignment history vector with another alignment history operation.
    ///
    /// This method panics in debug mode if the alignment history cannot be extended with the given alignment history operation.
    /// In release mode, this error is silently ignored.
    pub fn extend(
        &mut self,
        node_id: HistoryNodeIndex,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> HistoryNodeIndex {
        debug_assert!(self.can_extend(node_id, alignment, anchor_k, max_anchor_mutations));

        let node = &self.nodes[node_id];
        let cache = match alignment {
            AlignmentHistoryOperation::Match => node.match_successor,
            AlignmentHistoryOperation::Substitution => node.substitution_successor,
            AlignmentHistoryOperation::GapA => node.gap_a_successor,
            AlignmentHistoryOperation::GapB => node.gap_b_successor,
        };

        if cache.is_unknown() {
            let history = node
                .history
                .extend(alignment, anchor_k, max_anchor_mutations);
            self.nodes
                .push_in_place(|index| HistoryNode::new(index, history))
        } else {
            cache
        }
    }
}

impl<AlignmentHistoryVec> HistoryNode<AlignmentHistoryVec> {
    pub fn new(id: HistoryNodeIndex, history: AlignmentHistoryVec) -> Self {
        Self {
            id,
            history,
            match_successor: HistoryNodeIndex::new_unknown(),
            substitution_successor: HistoryNodeIndex::new_unknown(),
            gap_a_successor: HistoryNodeIndex::new_unknown(),
            gap_b_successor: HistoryNodeIndex::new_unknown(),
        }
    }
}

impl HistoryNodeIndex {
    /// Creates a history node index that implies that the extendability into this node is unknown.
    pub fn new_unknown() -> Self {
        Self(usize::MAX)
    }

    /// Creates a history node index that implies that extending into this node would result in a valid anchor.
    pub fn new_dead_end() -> Self {
        Self(usize::MAX - 1)
    }

    /// Returns true if the extendability into this node is unknown.
    pub fn is_unknown(&self) -> bool {
        self.0 == usize::MAX
    }

    /// Returns true if extending into this node would result in a valid anchor.
    pub fn is_dead_end(&self) -> bool {
        self.0 == usize::MAX - 1
    }
}

impl From<usize> for HistoryNodeIndex {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<HistoryNodeIndex> for usize {
    fn from(value: HistoryNodeIndex) -> Self {
        value.0
    }
}
