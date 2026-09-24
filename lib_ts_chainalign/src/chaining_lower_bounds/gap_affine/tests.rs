use generic_a_star::cost::{AStarCost, U32Cost};
use ndarray::Array2;

use crate::{chaining_lower_bounds::gap_affine::GapAffineLowerBounds, costs::GapAffineCosts};

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
