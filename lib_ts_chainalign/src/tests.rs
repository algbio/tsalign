use generic_a_star::cost::U32Cost;

use crate::{
    alignment::{
        AlignmentType,
        coordinates::{PrimaryAlignmentCoordinates, range::PrimaryAlignmentRange},
        sequences::AlignmentSequences,
    },
    anchors::Anchors,
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::{ChainingLowerBounds, gap_affine::GapAffineLowerBounds},
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, GapAffineCosts, TsBaseCost, TsLimits},
    max_match_run_chaining::gap_affine::GapAffineAligner,
    panic_on_extend::PanicOnExtend,
};

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => panic!("Unsupported character: {c}"),
    }
}

#[test]
fn lower_bound_violation() {
    let sequences = AlignmentSequences::new_complete(
        b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".into(),
        b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT".into(),
    );
    let config = ChainingLowerBoundConfig {
        anchor_k: 10,
        max_anchor_mutations: 1,
        inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
    };
    let costs = AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    );

    let lower_bounds = GapAffineLowerBounds::new_inexact_anchors(32, &config, &costs.primary_costs);
    assert_eq!(lower_bounds.lower_bound(10, 10), U32Cost::from(0u32));

    let lower_bounds = ChainingLowerBounds::new(32, config, costs.clone());
    assert_eq!(
        lower_bounds.primary_lower_bound(10, 10),
        U32Cost::from(0u32)
    );

    let mut gap_affine_aligner =
        GapAffineAligner::new(&sequences, &costs.primary_costs, &rc_fn, 9, 10);
    let mut additional_primary_targets_output = Vec::new();
    let (cost, alignment) = gap_affine_aligner.align(
        PrimaryAlignmentCoordinates::new(12, 12),
        PrimaryAlignmentCoordinates::new(22, 22),
        &mut additional_primary_targets_output,
        &mut PanicOnExtend,
    );
    assert_eq!(cost, U32Cost::from(0u32));
    assert_eq!(alignment.alignment, vec![(10, AlignmentType::Match)]);

    additional_primary_targets_output.clear();
    gap_affine_aligner.align_until_cost_limit(
        PrimaryAlignmentCoordinates::new(12, 12),
        U32Cost::from(0u32),
        &mut additional_primary_targets_output,
        &mut PanicOnExtend,
    );

    let anchors = Anchors::new_inexact(&sequences, 10, 1, &costs, &rc_fn);
    println!(
        "Anchors: {:?}",
        anchors
            .enumerate_primaries()
            .map(|(_, anchor)| anchor)
            .collect::<Vec<_>>()
    );
    let from_anchor = anchors
        .primary_index_from_range(PrimaryAlignmentRange::new_from_ranges(2..12, 2..12))
        .unwrap();
    let to_anchor = anchors
        .primary_index_from_range(PrimaryAlignmentRange::new_from_ranges(22..32, 22..32))
        .unwrap();

    let cost_function = ChainingCostFunction::new_from_lower_bounds(
        &lower_bounds,
        &anchors,
        &sequences,
        0u32.into(),
        &rc_fn,
    );
    assert_eq!(
        cost_function.primary(from_anchor, to_anchor),
        U32Cost::from(0u32)
    );
}
