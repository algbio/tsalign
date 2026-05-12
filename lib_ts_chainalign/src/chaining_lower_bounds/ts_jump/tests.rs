use generic_a_star::cost::{AStarCost, U32Cost};

use crate::{
    chaining_lower_bounds::ts_jump::TsJumpLowerBounds,
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
};

#[test]
fn test_max_match_run_0() {
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts {
            substitution: U32Cost::from(2u8),
            gap_open: U32Cost::from(3u8),
            gap_extend: U32Cost::from(1u8),
        },
        secondary_costs: GapAffineCosts {
            substitution: U32Cost::from(4u8),
            gap_open: U32Cost::from(6u8),
            gap_extend: U32Cost::from(2u8),
        },
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..200,
            ancestor_gap: -100..100,
        },
    };

    let max_n = 2;
    let lower_bounds = TsJumpLowerBounds::new(max_n, 0, &cost_table);

    let expected_lower_bounds_12 = [2, 4, 6];
    let expected_lower_bounds_34 =
        expected_lower_bounds_12.map(|cost| cost - cost_table.ts_base_cost.min().as_primitive());

    for descendant_gap in 0..=max_n {
        assert_eq!(
            lower_bounds.lower_bound_12(descendant_gap).as_primitive(),
            expected_lower_bounds_12[descendant_gap],
            "lower_bound_12({})",
            descendant_gap
        );
        assert_eq!(
            lower_bounds.lower_bound_34(descendant_gap).as_primitive(),
            expected_lower_bounds_34[descendant_gap],
            "lower_bound_34({})",
            descendant_gap
        );
    }
}

#[test]
fn test_max_match_run_1() {
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts {
            substitution: U32Cost::from(2u8),
            gap_open: U32Cost::from(3u8),
            gap_extend: U32Cost::from(1u8),
        },
        secondary_costs: GapAffineCosts {
            substitution: U32Cost::from(4u8),
            gap_open: U32Cost::from(6u8),
            gap_extend: U32Cost::from(2u8),
        },
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..200,
            ancestor_gap: -100..100,
        },
    };

    let max_n = 8;
    let lower_bounds = TsJumpLowerBounds::new(max_n, 1, &cost_table);

    let expected_lower_bounds_12 = [2, 2, 2, 4, 4, 6, 6, 8, 8];
    let expected_lower_bounds_34 =
        expected_lower_bounds_12.map(|cost| cost - cost_table.ts_base_cost.min().as_primitive());

    for descendant_gap in 0..=max_n {
        assert_eq!(
            lower_bounds.lower_bound_12(descendant_gap).as_primitive(),
            expected_lower_bounds_12[descendant_gap],
            "lower_bound_12({})",
            descendant_gap
        );
        assert_eq!(
            lower_bounds.lower_bound_34(descendant_gap).as_primitive(),
            expected_lower_bounds_34[descendant_gap],
            "lower_bound_34({})",
            descendant_gap
        );
    }
}

#[test]
fn test_max_match_run_2() {
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts {
            substitution: U32Cost::from(2u8),
            gap_open: U32Cost::from(3u8),
            gap_extend: U32Cost::from(1u8),
        },
        secondary_costs: GapAffineCosts {
            substitution: U32Cost::from(4u8),
            gap_open: U32Cost::from(6u8),
            gap_extend: U32Cost::from(2u8),
        },
        ts_base_cost: U32Cost::from(2u8).into(),
        ts_limits: TsLimits {
            inter_jump_12: -100..100,
            intra_jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..200,
            ancestor_gap: -100..100,
        },
    };

    let max_n = 9;
    let lower_bounds = TsJumpLowerBounds::new(max_n, 2, &cost_table);

    let expected_lower_bounds_12 = [2, 2, 2, 2, 2, 4, 4, 4, 6, 6];
    let expected_lower_bounds_34 =
        expected_lower_bounds_12.map(|cost| cost - cost_table.ts_base_cost.min().as_primitive());

    for descendant_gap in 0..=max_n {
        assert_eq!(
            lower_bounds.lower_bound_12(descendant_gap).as_primitive(),
            expected_lower_bounds_12[descendant_gap],
            "lower_bound_12({})",
            descendant_gap
        );
        assert_eq!(
            lower_bounds.lower_bound_34(descendant_gap).as_primitive(),
            expected_lower_bounds_34[descendant_gap],
            "lower_bound_34({})",
            descendant_gap
        );
    }
}
