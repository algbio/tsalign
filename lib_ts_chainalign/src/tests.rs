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

/// Align with tschainalign and return the cost and the number of template switches.
fn chain_align(
    reference: &[u8],
    query: &[u8],
    k: u32,
    max_anchor_mutations: u8,
    costs: &AlignmentCosts<U32Cost>,
) -> (u32, u32) {
    use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;

    use crate::chain_align::performance_parameters::AlignmentPerformanceParameters;

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
    let performance_parameters = AlignmentPerformanceParameters {
        max_exact_cost_function_cost: U32Cost::from(0u32),
        ..Default::default()
    };
    let result = crate::align::<DnaAlphabet>(
        reference.to_vec(),
        query.to_vec(),
        AlignmentRange::new_complete(reference.len(), query.len()),
        max_anchor_mutations,
        &performance_parameters,
        &rc_fn,
        "reference",
        "query",
        &lower_bounds,
    );
    let statistics = result.statistics();
    (
        statistics.cost.raw() as u32,
        statistics.template_switch_amount.raw() as u32,
    )
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().rev().copied().map(rc_fn).collect()
}

/// Costs with template switches of cost 2 and unlimited geometry.
fn ts_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([2u32, 2, 2, 2].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    )
}

/// Costs where template switches are never worth it.
fn no_ts_costs() -> AlignmentCosts<U32Cost> {
    AlignmentCosts::new(
        GapAffineCosts::new(2u32, 3, 1).into_costs::<U32Cost>(),
        GapAffineCosts::new(4u32, 3, 2).into_costs::<U32Cost>(),
        TsBaseCost::new([10_000u32; 4].map(U32Cost::from)),
        TsLimits::new_unlimited(),
    )
}

/// The optimal gap-affine alignment cost with substitution cost 2, gap open cost 3 and gap extend cost 1,
/// where the gap open cost is paid for the first character of a gap.
fn gap_affine_cost(a: &[u8], b: &[u8]) -> u32 {
    const INFINITY: u32 = u32::MAX / 4;
    let (n, m) = (a.len(), b.len());
    // Best cost ending with a match/substitution (or at the origin), a deletion, or an insertion.
    let mut diagonal = vec![vec![INFINITY; m + 1]; n + 1];
    let mut deletion = vec![vec![INFINITY; m + 1]; n + 1];
    let mut insertion = vec![vec![INFINITY; m + 1]; n + 1];
    diagonal[0][0] = 0;
    for i in 0..=n {
        for j in 0..=m {
            if i > 0 && j > 0 {
                let best = diagonal[i - 1][j - 1]
                    .min(deletion[i - 1][j - 1])
                    .min(insertion[i - 1][j - 1]);
                diagonal[i][j] = best + if a[i - 1] == b[j - 1] { 0 } else { 2 };
            }
            if i > 0 {
                deletion[i][j] =
                    (diagonal[i - 1][j].min(insertion[i - 1][j]) + 3).min(deletion[i - 1][j] + 1);
            }
            if j > 0 {
                insertion[i][j] =
                    (diagonal[i][j - 1].min(deletion[i][j - 1]) + 3).min(insertion[i][j - 1] + 1);
            }
        }
    }
    diagonal[n][m].min(deletion[n][m]).min(insertion[n][m])
}

/// Sequences that are shorter than k, where the optimal alignment starts or ends with matches.
/// The lower bounds for gaps from the start or to the end of the alignment used to assume
/// that the gap starts and ends with a non-match, which made them too large here.
#[test]
fn short_sequences_with_matches_at_start_and_end() {
    let costs = no_ts_costs();
    for (reference, query, expected) in [
        (&b"CCCTAGTT"[..], &b"CCCAAGTT"[..], 2),
        (b"A", b"TTTTTTTTTATTTTTTTTTT", 23),
        (b"ACGTACGTA", b"ACGTACGTA", 0),
        (b"ACGTACGTAC", b"ACGTACGTA", 3),
    ] {
        assert_eq!(gap_affine_cost(reference, query), expected);
        for k in [4, 6, 9, 12] {
            let (cost, _) = chain_align(reference, query, k, 0, &costs);
            assert_eq!(
                cost,
                expected,
                "k = {k}, reference = {}, query = {}",
                String::from_utf8_lossy(reference),
                String::from_utf8_lossy(query)
            );
        }
    }
}

/// A template switch right after matches at the start, or right before matches at the end of the alignment.
/// The optimal cost is 2 (one template switch), since the sequences differ and every operation costs at least 2.
/// The lower bounds for 12-jumps from the start and 34-jumps to the end used to assume
/// that the primary part of the jump starts (or ends) with a non-match.
///
/// If the primary part between the start and the template switch (or between the template switch and the end)
/// is shorter than k, tschainalign used to shorten the template switch
/// and align e.g. `3=1X[TS]` with cost 4 instead of `2=[TS]` with cost 2,
/// because the 12-jump (34-jump) aligner required a non-match in every jump with a non-empty descendant gap,
/// also when the jump starts at the start (ends at the end) of the alignment, which is no anchor
/// (see `ts_12_jump::tests::matches_from_alignment_start_then_jump` and
/// `ts_34_jump::tests::jump_then_matches_to_alignment_end`).
#[test]
fn template_switch_next_to_matches_at_start_and_end() {
    let costs = ts_costs();
    let inner = b"GATTACAGGCTTAC";
    let flank = b"CC";

    for (reference, query) in [
        // Matches at the start, template switch until the end.
        (
            [&flank[..], inner].concat(),
            [&flank[..], &reverse_complement(inner)].concat(),
        ),
        // Template switch from the start, matches at the end.
        (
            [&inner[..], flank].concat(),
            [&reverse_complement(inner)[..], flank].concat(),
        ),
        // Matches at both ends.
        (
            [&flank[..], inner, flank].concat(),
            [&flank[..], &reverse_complement(inner), flank].concat(),
        ),
    ] {
        for k in [4, 6] {
            let (cost, template_switches) = chain_align(&reference, &query, k, 0, &costs);
            assert_eq!(
                (cost, template_switches),
                (2, 1),
                "k = {k}, reference = {}, query = {}",
                String::from_utf8_lossy(&reference),
                String::from_utf8_lossy(&query)
            );
        }
    }
}

/// Random short sequences without template switches must be aligned optimally.
#[test]
fn random_short_sequences_without_template_switches() {
    let costs = no_ts_costs();
    let mut state = 0x5eed_u64;
    let mut next = |limit: usize| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) % limit as u64) as usize
    };

    let mut failures = Vec::new();
    for _ in 0..60 {
        let alphabet: &[u8] = if next(2) == 0 { b"AC" } else { b"ACGT" };
        let reference: Vec<u8> = (0..next(14) + 1)
            .map(|_| alphabet[next(alphabet.len())])
            .collect();
        let query: Vec<u8> = reference
            .iter()
            .flat_map(|&c| match next(8) {
                0 => vec![alphabet[next(alphabet.len())]],
                1 => vec![],
                2 => vec![c, alphabet[next(alphabet.len())]],
                _ => vec![c],
            })
            .collect();
        if query.is_empty() {
            continue;
        }
        let expected = gap_affine_cost(&reference, &query);

        for k in [3, 5, 8, 12] {
            let result = std::panic::catch_unwind(|| chain_align(&reference, &query, k, 0, &costs));
            match result {
                Ok((cost, _)) if cost == expected => {}
                Ok((cost, _)) => failures.push(format!(
                    "k = {k}, {} vs {}: cost {cost}, expected {expected}",
                    String::from_utf8_lossy(&reference),
                    String::from_utf8_lossy(&query)
                )),
                Err(_) => failures.push(format!(
                    "k = {k}, {} vs {}: panicked",
                    String::from_utf8_lossy(&reference),
                    String::from_utf8_lossy(&query)
                )),
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Align with tschainalign and return the reported alignment.
fn chain_align_result(
    reference: &[u8],
    query: &[u8],
    k: u32,
    costs: &AlignmentCosts<U32Cost>,
) -> lib_tsalign::a_star_aligner::alignment_result::alignment::Alignment<
    lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType,
> {
    use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
    use lib_tsalign::a_star_aligner::{
        alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    };

    use crate::chain_align::performance_parameters::AlignmentPerformanceParameters;

    let max_n = reference.len().max(query.len()).max(1).next_power_of_two();
    let lower_bounds = ChainingLowerBounds::new(
        max_n,
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations: 0,
            inexact_lower_bound: InexactLowerBoundKind::MaxMatchRun,
        },
        costs.clone(),
    );
    let performance_parameters = AlignmentPerformanceParameters {
        max_exact_cost_function_cost: U32Cost::from(0u32),
        ..Default::default()
    };
    match crate::align::<DnaAlphabet>(
        reference.to_vec(),
        query.to_vec(),
        AlignmentRange::new_complete(reference.len(), query.len()),
        0,
        &performance_parameters,
        &rc_fn,
        "reference",
        "query",
        &lower_bounds,
    ) {
        AlignmentResult::WithTarget { alignment, .. } => alignment,
        AlignmentResult::WithoutTarget { .. } => panic!("No alignment found"),
    }
}

/// The amount of reference and query characters consumed by an alignment,
/// following the semantics of lib_tsalign: after a template switch, the descendant continues after the characters
/// consumed by the template switch, and the anti-descendant continues at its entrance index plus the anti-descendant gap.
fn consumed_lengths(
    alignment: &lib_tsalign::a_star_aligner::alignment_result::alignment::Alignment<
        lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType,
    >,
) -> (usize, usize) {
    use lib_tsalign::a_star_aligner::template_switch_distance::{
        AlignmentType as T, TemplateSwitchDescendant,
    };

    let (mut reference, mut query) = (0isize, 0isize);
    let mut descendant_is_reference = false;
    let mut entrance_anti_descendant = 0isize;
    for alignment_type in alignment.iter_flat_cloned() {
        match alignment_type {
            T::PrimaryMatch | T::PrimarySubstitution => {
                reference += 1;
                query += 1;
            }
            T::PrimaryInsertion => query += 1,
            T::PrimaryDeletion => reference += 1,
            T::TemplateSwitchEntrance { descendant, .. } => {
                descendant_is_reference = descendant == TemplateSwitchDescendant::Reference;
                entrance_anti_descendant = if descendant_is_reference {
                    query
                } else {
                    reference
                };
            }
            T::SecondaryMatch | T::SecondarySubstitution | T::SecondaryInsertion => {
                if descendant_is_reference {
                    reference += 1;
                } else {
                    query += 1;
                }
            }
            T::SecondaryDeletion => {}
            T::TemplateSwitchExit {
                anti_descendant_gap,
            } => {
                if descendant_is_reference {
                    query = entrance_anti_descendant + anti_descendant_gap;
                } else {
                    reference = entrance_anti_descendant + anti_descendant_gap;
                }
            }
            other => panic!("Unexpected alignment type {other:?}"),
        }
    }
    (
        usize::try_from(reference).unwrap(),
        usize::try_from(query).unwrap(),
    )
}

/// Intra template switches (the ancestor is the descendant sequence) that start off the main diagonal.
/// The anti-descendant gap of such a template switch used to be computed relative to the descendant
/// instead of the anti-descendant, such that the reported alignment did not consume the sequences completely.
#[test]
fn intra_template_switches_off_the_main_diagonal() {
    let costs = ts_costs();
    let mut state = 0x1257_u64;
    let mut next = |limit: usize| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) % limit as u64) as usize
    };

    let mut failures = Vec::new();
    for _ in 0..30 {
        let reference: Vec<u8> = (0..30 + next(10)).map(|_| b"ACGT"[next(4)]).collect();
        // Insert characters before the template switch, such that it starts off the main diagonal.
        let insertion: Vec<u8> = (0..1 + next(4)).map(|_| b"ACGT"[next(4)]).collect();
        let split = 5 + next(5);
        let mut query = [&reference[..split], &insertion, &reference[split..]].concat();
        // Replace a segment of the query by the reverse complement of an earlier segment of the query.
        let length = 8 + next(5);
        let target = query.len() - length - next(4);
        let source = next(target - length);
        let inner = reverse_complement(&query[source..source + length]);
        query[target..target + length].copy_from_slice(&inner);

        for (reference, query) in [(&reference, &query), (&query, &reference)] {
            for k in [4, 5] {
                let result =
                    std::panic::catch_unwind(|| chain_align_result(reference, query, k, &costs));
                let context = format!(
                    "k = {k}, reference = {}, query = {}",
                    String::from_utf8_lossy(reference),
                    String::from_utf8_lossy(query)
                );
                match result {
                    Ok(alignment) => {
                        let lengths = consumed_lengths(&alignment);
                        if lengths != (reference.len(), query.len()) {
                            failures.push(format!(
                                "{context}: alignment {alignment} consumes {lengths:?}"
                            ));
                        }
                    }
                    Err(_) => failures.push(format!("{context}: panicked")),
                }
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// A template switch directly after a primary prefix and/or directly before a primary suffix,
/// for all prefix/suffix lengths from 1 to k + 1.
///
/// The prefix and suffix are equal in both sequences and the inner part of the query is the reverse complement
/// of the inner part of the reference, so a single template switch over the inner part aligns the sequences
/// with cost 2 (the template switch base cost).
/// This is optimal, because the sequences differ and every operation costs at least 2.
///
/// If the prefix (suffix) is at least k long, then there is a primary anchor directly before (after) the template switch
/// and the template switch is entered (exited) by a jump between anchors.
/// If it is shorter than k, there is no primary anchor, and the template switch is entered by a 12-jump
/// from the start of the alignment that consists of matches only (or exited by a 34-jump to the end of the alignment
/// that consists of matches only).
/// The 12-/34-jump aligners used to require a non-match in such jumps,
/// so the cases with a prefix or suffix of length 1 to k - 2 failed (cost 4 with one such flank, 6 with two),
/// both for exact anchors and for anchors with one mutation.
/// Length k - 1 worked by accident: the inner part and its reverse complement both start with `G` and end with `C`,
/// so the primary matches next to the template switch can be extended to length k, which gives a primary anchor
/// (with a template switch that is one shorter, but still costs 2).
#[test]
fn template_switch_next_to_flanks_shorter_than_k() {
    let costs = ts_costs();
    let inner = b"GATTACAGGCTTAC";
    let flank_source = b"ACCTGACT";

    let mut failures = Vec::new();
    for (k, max_anchor_mutations) in [(4u32, 0u8), (6, 0), (4, 1)] {
        for flank_length in 1..=k as usize + 1 {
            let flank = &flank_source[..flank_length];
            for (prefix, suffix) in [(flank, &b""[..]), (b"", flank), (flank, flank)] {
                let reference = [prefix, inner, suffix].concat();
                let query = [prefix, &reverse_complement(inner), suffix].concat();
                let (cost, template_switches) =
                    chain_align(&reference, &query, k, max_anchor_mutations, &costs);
                if (cost, template_switches) != (2, 1) {
                    failures.push(format!(
                        "k = {k}, m = {max_anchor_mutations}, prefix = {}, suffix = {}: cost {cost}, {template_switches} template switches",
                        String::from_utf8_lossy(prefix),
                        String::from_utf8_lossy(suffix),
                    ));
                }
            }
        }
    }

    assert!(failures.is_empty(), "Failures:\n{}", failures.join("\n"));
}
