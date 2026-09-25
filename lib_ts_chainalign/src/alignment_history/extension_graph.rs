use std::fmt::Display;

use tagged_vec::TaggedVec;

use crate::alignment_history::{
    history_alignment_operations::AlignmentHistoryOperation, history_vec::AlignmentHistory,
};

pub struct HistoryExtensionGraph<AlignmentHistoryVec> {
    nodes: TaggedVec<HistoryExtensionGraphNodeIndex, HistoryNode<AlignmentHistoryVec>>,
}

pub struct HistoryNode<AlignmentHistoryVec> {
    history: AlignmentHistoryVec,

    /// The successor of this node if the next alignment operation is a match.
    match_successor: HistoryExtensionGraphNodeIndex,

    /// The successor of this node if the next alignment operation is a substitution.
    substitution_successor: HistoryExtensionGraphNodeIndex,

    /// The successor of this node if the next alignment operation is a gap in a, i.e. sequence a misses a character.
    gap_in_a_successor: HistoryExtensionGraphNodeIndex,

    /// The successor of this node if the next alignment operation is a gap in b, i.e. sequence b misses a character.
    gap_in_b_successor: HistoryExtensionGraphNodeIndex,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct HistoryExtensionGraphNodeIndex(pub usize);

impl<AlignmentHistoryVec: AlignmentHistory> HistoryExtensionGraph<AlignmentHistoryVec> {
    pub fn new() -> Self {
        let mut nodes = TaggedVec::new();
        nodes.push(HistoryNode {
            history: AlignmentHistoryVec::default(),
            match_successor: HistoryExtensionGraphNodeIndex(usize::MAX),
            substitution_successor: HistoryExtensionGraphNodeIndex(usize::MAX),
            gap_in_a_successor: HistoryExtensionGraphNodeIndex(usize::MAX),
            gap_in_b_successor: HistoryExtensionGraphNodeIndex(usize::MAX),
        });
        Self { nodes }
    }

    /// Returns the id of the node that represents the empty history.
    pub fn empty_node_id(&self) -> HistoryExtensionGraphNodeIndex {
        HistoryExtensionGraphNodeIndex(0)
    }

    /// Returns true if the alignment history vector can be extended with another alignment history operation without becoming a valid anchor.
    #[allow(dead_code)]
    pub fn can_extend(
        &mut self,
        node_id: HistoryExtensionGraphNodeIndex,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> bool {
        self.try_extend(node_id, alignment, anchor_k, max_anchor_mutations)
            .is_some()
    }

    /// Return the extension of this the alignment history vector with another alignment history operation.
    ///
    /// This method panics if the alignment history cannot be extended with the given alignment history operation.
    #[allow(dead_code)]
    pub fn extend(
        &mut self,
        node_id: HistoryExtensionGraphNodeIndex,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> HistoryExtensionGraphNodeIndex {
        if let Some(extension_index) =
            self.try_extend(node_id, alignment, anchor_k, max_anchor_mutations)
        {
            extension_index
        } else {
            panic!(
                "Cannot extend the alignment history vector with the given alignment history operation."
            )
        }
    }

    /// Try to extend the alignment history vector with another alignment history operation.
    ///
    /// Returns `Some` with the index of the new node if the extension is possible, or `None` if it would result in a valid anchor.
    pub fn try_extend(
        &mut self,
        node_id: HistoryExtensionGraphNodeIndex,
        alignment: AlignmentHistoryOperation,
        anchor_k: u8,
        max_anchor_mutations: u8,
    ) -> Option<HistoryExtensionGraphNodeIndex> {
        let nodes_len = self.nodes.len();
        let node = &mut self.nodes[node_id];
        let cache = match alignment {
            AlignmentHistoryOperation::Match => &mut node.match_successor,
            AlignmentHistoryOperation::Substitution => &mut node.substitution_successor,
            AlignmentHistoryOperation::GapInA => &mut node.gap_in_a_successor,
            AlignmentHistoryOperation::GapInB => &mut node.gap_in_b_successor,
        };

        if cache.is_dead_end() {
            None
        } else if cache.is_unknown() {
            if let Some(extension) =
                node.history
                    .try_extend(alignment, anchor_k, max_anchor_mutations)
            {
                *cache = nodes_len.into();
                let new_node_id = self.nodes.push(HistoryNode::new(extension));

                debug_assert_eq!(
                    {
                        let node = &self.nodes[node_id];
                        match alignment {
                            AlignmentHistoryOperation::Match => node.match_successor,
                            AlignmentHistoryOperation::Substitution => node.substitution_successor,
                            AlignmentHistoryOperation::GapInA => node.gap_in_a_successor,
                            AlignmentHistoryOperation::GapInB => node.gap_in_b_successor,
                        }
                    },
                    new_node_id
                );

                Some(new_node_id)
            } else {
                *cache = HistoryExtensionGraphNodeIndex::new_dead_end();
                None
            }
        } else {
            Some(*cache)
        }
    }
}

impl<AlignmentHistoryVec> HistoryNode<AlignmentHistoryVec> {
    pub fn new(history: AlignmentHistoryVec) -> Self {
        Self {
            history,
            match_successor: HistoryExtensionGraphNodeIndex::new_unknown(),
            substitution_successor: HistoryExtensionGraphNodeIndex::new_unknown(),
            gap_in_a_successor: HistoryExtensionGraphNodeIndex::new_unknown(),
            gap_in_b_successor: HistoryExtensionGraphNodeIndex::new_unknown(),
        }
    }
}

impl HistoryExtensionGraphNodeIndex {
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

impl From<usize> for HistoryExtensionGraphNodeIndex {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<HistoryExtensionGraphNodeIndex> for usize {
    fn from(value: HistoryExtensionGraphNodeIndex) -> Self {
        value.0
    }
}

impl Display for HistoryExtensionGraphNodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<AlignmentHistoryVec: AlignmentHistory> Default for HistoryExtensionGraph<AlignmentHistoryVec> {
    fn default() -> Self {
        Self::new()
    }
}
