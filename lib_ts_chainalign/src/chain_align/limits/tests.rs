use generic_a_star::{AStarResult, cost::U32Cost};
use lib_tsalign::{alignment_geometry::AlignmentRange, alignment_result::AlignmentResult};

use crate::{
    chain_align::{
        limits::AlignmentLimits, performance_parameters::AlignmentPerformanceParameters,
    },
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    dna_rc_fn,
};

/// Costs with template switches of cost 2 and unlimited geometry.
fn ts_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    )
}

/// The optimal alignment contains a template switch (cost 8 with exact anchors and k = 4).
const TS_REFERENCE: &[u8] = b"TCTAACTCCAGTCAACTTCCGAGAGCACCA";
const TS_QUERY: &[u8] = b"TCTAACTCCAGGAAGTAAACGAGAGATCCA";

/// With exact anchors, k = 4 and the chaining cost function initialised with lower bounds only,
/// finding the optimal chain takes four chaining rounds.
const ROUNDS_REFERENCE: &[u8] = b"GGACATTCATCAAGGCCCAGTACCTAATTC";
const ROUNDS_QUERY: &[u8] = b"GAACATTCATCTGGGCCTTGTGCCTAATTC";

fn align(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
    costs: &AlignmentCosts<U32Cost>,
    limits: &AlignmentLimits<U32Cost>,
) -> AlignmentResult<U32Cost> {
    let max_n = reference.len().max(query.len()).max(1).next_power_of_two();
    let lower_bounds = ChainingLowerBounds::new(
        max_n,
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        costs.clone(),
    );
    // Initialise the cost function with lower bounds only, such that the chaining takes several rounds.
    let performance_parameters = AlignmentPerformanceParameters {
        max_exact_cost_function_cost: U32Cost::from(0u32),
        ..Default::default()
    };
    crate::align_with_limits(
        reference.to_vec(),
        query.to_vec(),
        AlignmentRange::new_complete(reference.len(), query.len()),
        &performance_parameters,
        limits,
        &dna_rc_fn,
        "reference",
        "query",
        &lower_bounds,
    )
}

fn cost_limit(cost_limit: u32) -> AlignmentLimits<U32Cost> {
    AlignmentLimits {
        cost_limit: Some(cost_limit.into()),
        ..AlignmentLimits::none()
    }
}

fn iteration_limit(iteration_limit: u64) -> AlignmentLimits<U32Cost> {
    AlignmentLimits {
        iteration_limit: Some(iteration_limit),
        ..AlignmentLimits::none()
    }
}

fn chain_cost(result: &AlignmentResult<U32Cost>) -> u32 {
    result.statistics().chain_cost.unwrap().raw() as u32
}

/// The cost limit is inclusive: a limit equal to the optimal chain cost returns the same alignment as no limit,
/// and a limit one below reports that the limit was exceeded, without an alignment.
#[test]
fn cost_limit_is_inclusive() {
    let costs = ts_costs();
    for (reference, query, k, max_anchor_mutations) in [
        (TS_REFERENCE, TS_QUERY, 4, 0),
        (TS_REFERENCE, TS_QUERY, 6, 1),
        (ROUNDS_REFERENCE, ROUNDS_QUERY, 4, 0),
        (ROUNDS_REFERENCE, ROUNDS_QUERY, 6, 1),
    ] {
        let description = format!(
            "k = {k}, max_anchor_mutations = {max_anchor_mutations}, {} vs {}",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(query)
        );
        let unlimited = align(
            reference,
            query,
            k,
            max_anchor_mutations,
            &costs,
            &AlignmentLimits::none(),
        );
        assert!(
            matches!(unlimited, AlignmentResult::WithTarget { .. }),
            "{description}"
        );
        let optimum = chain_cost(&unlimited);
        assert!(optimum > 0, "{description}");

        for limit in [optimum, optimum + 1, u32::MAX] {
            let limited = align(
                reference,
                query,
                k,
                max_anchor_mutations,
                &costs,
                &cost_limit(limit),
            );
            assert!(
                matches!(limited, AlignmentResult::WithTarget { .. }),
                "{description}, limit {limit}: {limited}"
            );
            assert_eq!(chain_cost(&limited), optimum, "{description}");
            assert_eq!(
                limited.statistics().cost,
                unlimited.statistics().cost,
                "{description}"
            );
            assert_eq!(limited.cigar(), unlimited.cigar(), "{description}");
        }

        for limit in [0, optimum / 2, optimum - 1] {
            let limited = align(
                reference,
                query,
                k,
                max_anchor_mutations,
                &costs,
                &cost_limit(limit),
            );
            let AlignmentResult::WithoutTarget { statistics } = &limited else {
                panic!("{description}, limit {limit}: found an alignment: {limited}");
            };
            assert_eq!(
                statistics.result,
                AStarResult::ExceededCostLimit {
                    cost_limit: U32Cost::from(limit)
                },
                "{description}, limit {limit}"
            );
            assert_eq!(statistics.exceeded_iteration_limit, None);
            assert_eq!(statistics.chain_cost, None);
            assert!(
                limited.to_string().contains("Exceeded cost limit"),
                "{limited}"
            );
        }
    }
}

/// Identical sequences have cost zero, which a cost limit of zero allows.
#[test]
fn cost_limit_zero_for_identical_sequences() {
    let result = align(
        ROUNDS_REFERENCE,
        ROUNDS_REFERENCE,
        4,
        0,
        &ts_costs(),
        &cost_limit(0),
    );
    assert!(matches!(result, AlignmentResult::WithTarget { .. }));
    assert_eq!(result.statistics().cost.raw(), 0.0);
}

/// The cost limit applies to the chain cost, which tschainalign minimises.
/// For inexact anchors, the cost of the alignment may be lower than its chain cost (issue #221),
/// here 5 instead of 7 (see `double_gap_open_between_gap_and_inexact_anchor` in the crate tests).
#[test]
fn cost_limit_applies_to_the_chain_cost() {
    let costs = AlignmentCosts::new(
        GapAffineCosts::new(4u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([10_000u32; 4].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    );
    let reference = b"CATGCATTACGGGACT";
    let query = b"CATGCATTACACT";

    let result = align(reference, query, 4, 1, &costs, &cost_limit(7));
    assert!(matches!(result, AlignmentResult::WithTarget { .. }));
    assert_eq!(
        (result.statistics().cost.raw(), chain_cost(&result)),
        (5.0, 7)
    );

    for limit in [5, 6] {
        let result = align(reference, query, 4, 1, &costs, &cost_limit(limit));
        assert_eq!(
            result.statistics().result,
            AStarResult::ExceededCostLimit {
                cost_limit: U32Cost::from(limit)
            },
            "limit {limit}"
        );
    }
}

/// The iteration limit is the maximum number of chaining rounds.
/// If the optimal chain is confirmed within the limit, then the alignment is the same as without a limit.
/// Otherwise, no alignment is returned, and the result reports the iteration limit
/// and a lower bound of the chain cost, which does not decrease with the limit.
#[test]
fn iteration_limit_counts_chaining_rounds() {
    let costs = ts_costs();
    let unlimited = align(
        ROUNDS_REFERENCE,
        ROUNDS_QUERY,
        4,
        0,
        &costs,
        &AlignmentLimits::none(),
    );
    assert!(matches!(unlimited, AlignmentResult::WithTarget { .. }));
    let optimum = chain_cost(&unlimited);

    let mut previous_lower_bound = 0.0;
    for limit in 0..4 {
        let result = align(
            ROUNDS_REFERENCE,
            ROUNDS_QUERY,
            4,
            0,
            &costs,
            &iteration_limit(limit),
        );
        let AlignmentResult::WithoutTarget { statistics } = &result else {
            panic!("limit {limit}: found an alignment: {result}");
        };
        assert_eq!(statistics.exceeded_iteration_limit, Some(limit));
        assert_eq!(statistics.result, AStarResult::NoTarget);
        let lower_bound = statistics.cost.raw();
        assert!(
            previous_lower_bound <= lower_bound && lower_bound <= f64::from(optimum),
            "limit {limit}: lower bound {lower_bound}, previous {previous_lower_bound}, optimum {optimum}"
        );
        previous_lower_bound = lower_bound;
        assert!(
            result
                .to_string()
                .contains(&format!("Exceeded iteration limit of {limit}")),
            "{result}"
        );
    }
    // Zero rounds give no information about the cost.
    assert_eq!(
        align(
            ROUNDS_REFERENCE,
            ROUNDS_QUERY,
            4,
            0,
            &costs,
            &iteration_limit(0)
        )
        .statistics()
        .cost
        .raw(),
        0.0
    );

    for limit in [4, 5, u64::MAX] {
        let result = align(
            ROUNDS_REFERENCE,
            ROUNDS_QUERY,
            4,
            0,
            &costs,
            &iteration_limit(limit),
        );
        assert!(
            matches!(result, AlignmentResult::WithTarget { .. }),
            "limit {limit}: {result}"
        );
        assert_eq!(result.statistics().exceeded_iteration_limit, None);
        assert_eq!(chain_cost(&result), optimum);
        assert_eq!(result.cigar(), unlimited.cigar());
    }
}

/// Both limits together: the first limit that is reached is reported.
#[test]
fn cost_and_iteration_limit() {
    let costs = ts_costs();
    let optimum = chain_cost(&align(
        ROUNDS_REFERENCE,
        ROUNDS_QUERY,
        4,
        0,
        &costs,
        &AlignmentLimits::none(),
    ));

    let both = |cost_limit: u32, iteration_limit: u64| {
        align(
            ROUNDS_REFERENCE,
            ROUNDS_QUERY,
            4,
            0,
            &costs,
            &AlignmentLimits {
                cost_limit: Some(cost_limit.into()),
                iteration_limit: Some(iteration_limit),
            },
        )
    };

    assert!(matches!(
        both(optimum, 4),
        AlignmentResult::WithTarget { .. }
    ));
    assert_eq!(
        both(optimum, 3).statistics().exceeded_iteration_limit,
        Some(3)
    );
    // The searches of the first rounds find chains with a lower bound of at most the limit,
    // so the cost limit is reached in a later round.
    let result = both(optimum - 1, 100);
    assert_eq!(
        result.statistics().result,
        AStarResult::ExceededCostLimit {
            cost_limit: U32Cost::from(optimum - 1)
        }
    );
    assert_eq!(result.statistics().exceeded_iteration_limit, None);
}
