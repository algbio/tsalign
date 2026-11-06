use generic_a_star::cost::{AStarCost, U32Cost};
use ndarray::Array2;

use crate::{
    alignment_history::{
        history_alignment_operations::AlignmentHistoryOperation,
        history_vec::UnsignedIntAlignmentHistoryVec,
        test_util::{OPERATIONS, alignment_cost, ends_with_anchor},
    },
    chaining_lower_bounds::gap_affine::GapAffineLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::GapAffineCosts,
};

#[test]
fn exact_max_match_run_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 2;
    let lower_bounds = GapAffineLowerBounds::new_exact_anchors(max_n, 0, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            u32::MAX, 3, 4,
                   3, 2, 5,
                   4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_1() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds = GapAffineLowerBounds::new_exact_anchors(max_n, 1, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            u32::MAX, 3, 4, 5, 6,
                   3, 2, 5, 6, 7,
                   4, 5, 4, 5, 6,
                   5, 6, 5, 4, 7,
                   6, 7, 6, 7, 6,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_2() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds = GapAffineLowerBounds::new_exact_anchors(max_n, 2, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            u32::MAX, 3, 4, 5, 6, 7, 8,
                   3, 2, 5, 6, 7, 8, 9,
                   4, 5, 4, 5, 6, 7, 8,
                   5, 6, 5, 4, 5, 6, 7,
                   6, 7, 6, 5, 4, 7, 8,
                   7, 8, 7, 6, 7, 6, 7,
                   8, 9, 8, 7, 8, 7, 6,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_0_allow_start_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 2;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_start_match(max_n, 0, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4,
               3, 2, 5,
               4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_1_allow_start_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_start_match(max_n, 1, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6,
               3, 2, 3, 4, 5,
               4, 3, 2, 5, 6,
               5, 4, 5, 4, 5,
               6, 5, 6, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_2_allow_start_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_start_match(max_n, 2, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6, 7, 8,
               3, 2, 3, 4, 5, 6, 7,
               4, 3, 2, 3, 4, 5, 6,
               5, 4, 3, 2, 5, 6, 7,
               6, 5, 4, 5, 4, 5, 6,
               7, 6, 5, 6, 5, 4, 5,
               8, 7, 6, 7, 6, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_0_allow_end_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 2;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_end_match(max_n, 0, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4,
               3, 2, 5,
               4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_1_allow_end_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_end_match(max_n, 1, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6,
               3, 2, 3, 4, 5,
               4, 3, 2, 5, 6,
               5, 4, 5, 4, 5,
               6, 5, 6, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn exact_max_match_run_2_allow_end_match() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds =
        GapAffineLowerBounds::new_exact_anchors_allow_end_match(max_n, 2, &cost_table);

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6, 7, 8,
               3, 2, 3, 4, 5, 6, 7,
               4, 3, 2, 3, 4, 5, 6,
               5, 4, 3, 2, 5, 6, 7,
               6, 5, 4, 5, 4, 5, 6,
               7, 6, 5, 6, 5, 4, 5,
               8, 7, 6, 7, 6, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_max_match_run_1_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 2;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 1,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4,
               3, 0, 3,
               4, 3, 2,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_max_match_run_2_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 2,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6,
               3, 0, 3, 4, 5,
               4, 3, 0, 3, 4,
               5, 4, 3, 2, 3,
               6, 5, 4, 3, 2,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_max_match_run_3_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 3,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6, 7, 8,
               3, 0, 3, 4, 5, 6, 7,
               4, 3, 0, 3, 4, 5, 6,
               5, 4, 3, 0, 3, 4, 5,
               6, 5, 4, 3, 2, 3, 4,
               7, 6, 5, 4, 3, 2, 3,
               8, 7, 6, 5, 4, 3, 2,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_max_match_run_2_1() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 2,
            max_anchor_mutations: 1,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6,
               3, 0, 3, 4, 5,
               4, 3, 0, 3, 4,
               5, 4, 3, 2, 5,
               6, 5, 4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_max_match_run_3_1() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 3,
            max_anchor_mutations: 1,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0u32, 3, 4, 5, 6, 7, 8,
               3, 0, 3, 4, 5, 6, 7,
               4, 3, 0, 3, 4, 5, 6,
               5, 4, 3, 0, 3, 4, 5,
               6, 5, 4, 3, 2, 3, 4,
               7, 6, 5, 4, 3, 2, 5,
               8, 7, 6, 5, 4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_1_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 2;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 1,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0, 3, 4,
            3, 2, 5,
            4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_2_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 2,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0, 3, 4, 5, 6,
            3, 0, 3, 4, 5,
            4, 3, 2, 3, 4,
            5, 4, 3, 2, 5,
            6, 5, 4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_3_0() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 3,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0, 3, 4, 5, 6, 7, 8,
            3, 0, 3, 4, 5, 6, 7,
            4, 3, 0, 3, 4, 5, 6,
            5, 4, 3, 2, 3, 4, 5,
            6, 5, 4, 3, 2, 3, 4,
            7, 6, 5, 4, 3, 2, 5,
            8, 7, 6, 5, 4, 5, 4,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_2_1() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 4;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 2,
            max_anchor_mutations: 1,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0, 3, 4, 5, 6,
            3, 0, 5, 6, 7,
            4, 5, 4, 7, 8,
            5, 6, 7, 6, 9,
            6, 7, 8, 9, 8,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_3_1() {
    let cost_table = GapAffineCosts {
        substitution: U32Cost::from(2u8),
        gap_open: U32Cost::from(3u8),
        gap_extend: U32Cost::from(1u8),
    };
    let max_n = 6;
    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
        max_n,
        &ChainingLowerBoundConfig {
            anchor_k: 3,
            max_anchor_mutations: 1,
            inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
        },
        &cost_table,
    );

    #[rustfmt::skip]
    let exepcted_lower_bounds = Array2::from_shape_vec(
        (max_n+1, max_n+1),
        vec![
            0, 3, 4, 5, 6, 7, 8,
            3, 0, 3, 4, 5, 6, 7,
            4, 3, 0, 5, 4, 5, 6,
            5, 4, 5, 4, 5, 6, 7,
            6, 5, 4, 5, 4, 7, 8,
            7, 6, 5, 6, 7, 6, 9,
            8, 7, 6, 7, 8, 9, 8,
        ],
    )
    .unwrap();

    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                lower_bounds.lower_bound(a, b).as_primitive(),
                exepcted_lower_bounds[[a, b]],
                "lower bound({}, {})",
                a,
                b
            );
        }
    }
}

const WINDOWED_MIN_MUTATIONS_TEST_COSTS: [(u8, u8, u8); 6] = [
    (2, 3, 1),
    (1, 1, 1),
    (3, 2, 1),
    (1, 4, 1),
    (2, 2, 2),
    (2, 5, 0),
];

fn windowed_min_mutations_config(
    anchor_k: u8,
    max_anchor_mutations: u8,
) -> ChainingLowerBoundConfig {
    ChainingLowerBoundConfig {
        anchor_k: anchor_k.into(),
        max_anchor_mutations,
        inexact_lower_bound: InexactLowerBoundKind::WindowedMinMutations,
    }
}

fn assert_lower_bounds_eq(
    actual: &GapAffineLowerBounds<U32Cost>,
    expected: impl Fn(usize, usize) -> u32,
    max_n: usize,
    description: &str,
) {
    for a in 0..=max_n {
        for b in 0..=max_n {
            assert_eq!(
                actual.lower_bound(a, b).as_primitive(),
                expected(a, b),
                "lower bound({a}, {b}) with {description}"
            );
        }
        assert_eq!(
            actual.variable_gap2_lower_bound(a).as_primitive(),
            (0..=max_n).map(|b| expected(a, b)).min().unwrap(),
            "variable gap2 lower bound({a}) with {description}"
        );
    }
}

/// Computes the windowed-min-mutation lower bounds by enumerating all alignments that contain no anchor.
fn brute_force_windowed_min_mutations(
    max_n: usize,
    anchor_k: usize,
    max_anchor_mutations: usize,
    costs: &GapAffineCosts<u32>,
) -> Array2<u32> {
    fn recurse(
        alignment: &mut Vec<AlignmentHistoryOperation>,
        (a, b): (usize, usize),
        max_n: usize,
        anchor_k: usize,
        max_anchor_mutations: usize,
        costs: &GapAffineCosts<u32>,
        result: &mut Array2<u32>,
    ) {
        let cost = alignment_cost(alignment, costs);
        result[[a, b]] = result[[a, b]].min(cost);

        for operation in OPERATIONS {
            let next = (a + operation.len_a(), b + operation.len_b());
            if next.0 > max_n || next.1 > max_n {
                continue;
            }

            alignment.push(operation);
            if !ends_with_anchor(alignment, anchor_k, max_anchor_mutations) {
                recurse(
                    alignment,
                    next,
                    max_n,
                    anchor_k,
                    max_anchor_mutations,
                    costs,
                    result,
                );
            }
            alignment.pop();
        }
    }

    let mut result = Array2::from_elem((max_n + 1, max_n + 1), u32::MAX);
    recurse(
        &mut Vec::new(),
        (0, 0),
        max_n,
        anchor_k,
        max_anchor_mutations,
        costs,
        &mut result,
    );
    result
}

#[test]
fn inexact_windowed_min_mutations_brute_force() {
    for (substitution, gap_open, gap_extend) in WINDOWED_MIN_MUTATIONS_TEST_COSTS {
        let cost_table = GapAffineCosts::new(
            U32Cost::from(substitution),
            U32Cost::from(gap_open),
            U32Cost::from(gap_extend),
        );
        let primitive_cost_table = GapAffineCosts::new(
            u32::from(substitution),
            u32::from(gap_open),
            u32::from(gap_extend),
        );

        for (max_n, anchor_k, max_anchor_mutations) in [
            (4, 1, 0),
            (4, 2, 0),
            (4, 2, 1),
            (4, 3, 0),
            (4, 3, 1),
            (4, 3, 2),
            (5, 3, 1),
            (5, 4, 1),
            (5, 4, 2),
            (5, 4, 3),
            (5, 5, 2),
        ] {
            let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                &windowed_min_mutations_config(anchor_k, max_anchor_mutations),
                &cost_table,
            );
            let expected = brute_force_windowed_min_mutations(
                max_n,
                anchor_k.into(),
                max_anchor_mutations.into(),
                &primitive_cost_table,
            );
            assert_lower_bounds_eq(
                &lower_bounds,
                |a, b| expected[[a, b]],
                max_n,
                &format!(
                    "max_n = {max_n}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}"
                ),
            );
        }
    }
}

#[test]
fn inexact_windowed_min_mutations_pareto_equals_unpruned() {
    for (substitution, gap_open, gap_extend) in WINDOWED_MIN_MUTATIONS_TEST_COSTS {
        let cost_table = GapAffineCosts::new(
            U32Cost::from(substitution),
            U32Cost::from(gap_open),
            U32Cost::from(gap_extend),
        );

        for (max_n, anchor_k, max_anchor_mutations) in [
            (10, 1, 0),
            (10, 2, 0),
            (10, 2, 1),
            (10, 3, 0),
            (10, 3, 1),
            (10, 3, 2),
            (12, 4, 0),
            (12, 4, 1),
            (12, 4, 2),
            (12, 4, 3),
            (14, 5, 1),
            (14, 5, 2),
            (14, 6, 2),
        ] {
            let description = format!(
                "max_n = {max_n}, anchor_k = {anchor_k}, max_anchor_mutations = {max_anchor_mutations}, costs = {cost_table:?}"
            );
            let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
                max_n,
                &windowed_min_mutations_config(anchor_k, max_anchor_mutations),
                &cost_table,
            );
            let unpruned =
                GapAffineLowerBounds::compute_windowed_min_mutation_without_pareto_pruning::<
                    UnsignedIntAlignmentHistoryVec<u64>,
                >(max_n, anchor_k, max_anchor_mutations, &cost_table);
            assert_lower_bounds_eq(
                &lower_bounds,
                |a, b| unpruned.lower_bound(a, b).as_primitive(),
                max_n,
                &description,
            );
        }
    }
}

/// Checks properties of larger lower bounds that are too expensive to compute without Pareto pruning.
#[test]
fn inexact_windowed_min_mutations_properties() {
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    for (max_n, anchor_k, max_anchor_mutations) in [(64, 12, 3), (48, 16, 2), (40, 8, 7)] {
        let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(
            max_n,
            &windowed_min_mutations_config(anchor_k, max_anchor_mutations),
            &cost_table,
        );
        let anchor_k = usize::from(anchor_k);

        for a in 0..=max_n {
            for b in 0..=max_n {
                let lower_bound = lower_bounds.lower_bound(a, b).as_primitive();
                // Swapping the sequences swaps gaps in a and gaps in b.
                assert_eq!(lower_bound, lower_bounds.lower_bound(b, a).as_primitive());
                // At least the gap length difference has to be filled with gaps.
                if a != b {
                    assert!(lower_bound >= 3 + a.abs_diff(b) as u32 - 1);
                }
            }

            // Shorter than an anchor, so matches are possible.
            if a < anchor_k {
                assert_eq!(lower_bounds.lower_bound(a, a).as_primitive(), 0);
            } else {
                // The alignment contains `a / anchor_k` disjoint subalignments that each cover `anchor_k` characters of a.
                // Each needs more than `max_anchor_mutations` mutations, and each mutation costs at least `gap_extend = 1`.
                let windows = a / anchor_k;
                assert!(
                    lower_bounds.lower_bound(a, a).as_primitive()
                        >= windows as u32 * (u32::from(max_anchor_mutations) + 1),
                    "lower bound({a}, {a})"
                );
            }
        }
    }
}
