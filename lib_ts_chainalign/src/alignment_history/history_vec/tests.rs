use crate::alignment_history::history_vec::{AlignmentHistory, UnsignedIntAlignmentHistoryVec};

#[test]
fn required_capacity() {
    fn required_capacity(anchor_k: u8, max_anchor_mutations: u8) -> usize {
        UnsignedIntAlignmentHistoryVec::<u8>::required_capacity(anchor_k, max_anchor_mutations)
    }

    assert_eq!(required_capacity(8, 0), 8);
    assert_eq!(required_capacity(8, 1), 9);
    assert_eq!(required_capacity(8, 2), 9);
    assert_eq!(required_capacity(8, 3), 10);
    assert_eq!(required_capacity(8, 4), 10);

    assert_eq!(required_capacity(6, 0), 6);
    assert_eq!(required_capacity(6, 1), 7);
    assert_eq!(required_capacity(6, 2), 7);
    assert_eq!(required_capacity(6, 3), 8);
    assert_eq!(required_capacity(6, 4), 8);
}
