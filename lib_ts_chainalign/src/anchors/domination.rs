//! PROTOTYPE: statistics and pruning of dominated inexact anchors.
//!
//! Controlled by environment variables:
//! * `TSALIGN_ANCHOR_PRUNE`: `none` (default, only statistics), `start`, `end` or `both`.
//! * `TSALIGN_ANCHOR_PRUNE_SECONDARY`: `skip` (default, never prune secondary anchors),
//!   `desc` (prune secondary anchors only if the descendant coordinate of the side that differs is equal,
//!   so the TS length check on anchor spans is unaffected), or `all`.
//! * `TSALIGN_ANCHOR_PRUNE_STRICT`: if set to `1`, use `<` instead of `<=` in the domination relation.
//!
//! Relations (anchors are given in the coordinates of the generator, where both sequences are forward):
//! * R1 (same start): A1, A2 share the start, end(A1) <= end(A2) componentwise, end(A1) != end(A2),
//!   and cost(A1) + D(end(A1), end(A2)) <= cost(A2), where D is the unrestricted gap-affine alignment cost.
//!   Then A2 is dominated.
//! * R2 (same end): mirror image with starts.

use std::collections::HashMap;

use generic_a_star::cost::AStarCost;
use inexact_alignment_anchors::anchor::Anchor;
use log::info;

use crate::costs::GapAffineCosts;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PruneMode {
    None,
    Start,
    End,
    Both,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SecondaryPolicy {
    Skip,
    DescendantFixed,
    All,
}

pub fn config_from_env() -> (PruneMode, SecondaryPolicy, bool) {
    let mode = match std::env::var("TSALIGN_ANCHOR_PRUNE").as_deref() {
        Ok("start") => PruneMode::Start,
        Ok("end") => PruneMode::End,
        Ok("both") => PruneMode::Both,
        _ => PruneMode::None,
    };
    let secondary = match std::env::var("TSALIGN_ANCHOR_PRUNE_SECONDARY").as_deref() {
        Ok("desc") => SecondaryPolicy::DescendantFixed,
        Ok("all") => SecondaryPolicy::All,
        _ => SecondaryPolicy::Skip,
    };
    let strict = std::env::var("TSALIGN_ANCHOR_PRUNE_STRICT").as_deref() == Ok("1");
    (mode, secondary, strict)
}

/// Unrestricted global gap-affine alignment cost of `x` vs `y`.
fn gap_affine_distance<Cost: AStarCost>(x: &[u8], y: &[u8], costs: &GapAffineCosts<Cost>) -> Cost {
    let n = x.len();
    let m = y.len();
    let inf = Cost::max_value();
    let add = |a: Cost, b: Cost| if a == inf { inf } else { a + b };
    // M: last op diagonal (or start), I: last op gap in x (consume y), D: last op gap in y (consume x).
    let mut mm = vec![vec![inf; m + 1]; n + 1];
    let mut ii = vec![vec![inf; m + 1]; n + 1];
    let mut dd = vec![vec![inf; m + 1]; n + 1];
    mm[0][0] = Cost::zero();
    for i in 0..=n {
        for j in 0..=m {
            if i > 0 && j > 0 {
                let c = if x[i - 1] == y[j - 1] {
                    Cost::zero()
                } else {
                    costs.substitution
                };
                let best = mm[i - 1][j - 1].min(ii[i - 1][j - 1]).min(dd[i - 1][j - 1]);
                mm[i][j] = add(best, c);
            }
            if j > 0 {
                let open = mm[i][j - 1].min(dd[i][j - 1]);
                ii[i][j] = add(open, costs.gap_open).min(add(ii[i][j - 1], costs.gap_extend));
            }
            if i > 0 {
                let open = mm[i - 1][j].min(ii[i - 1][j]);
                dd[i][j] = add(open, costs.gap_open).min(add(dd[i - 1][j], costs.gap_extend));
            }
        }
    }
    mm[n][m].min(ii[n][m]).min(dd[n][m])
}

#[derive(Default, Debug)]
struct Stats {
    total: usize,
    distinct_starts: usize,
    distinct_ends: usize,
    in_shared_start: usize,
    in_shared_end: usize,
    max_start_group: usize,
    max_end_group: usize,
    zero_cost: usize,
    r1_le: usize,
    r1_lt: usize,
    r2_le: usize,
    r2_lt: usize,
    either_le: usize,
    pruned: usize,
}

/// Compute statistics and possibly prune dominated anchors.
pub fn process<Cost: AStarCost>(
    label: &str,
    is_secondary: bool,
    anchors: Vec<Anchor<Cost>>,
    seq_a: &[u8],
    seq_b: &[u8],
    costs: &GapAffineCosts<Cost>,
) -> Vec<Anchor<Cost>> {
    let (mode, secondary_policy, strict) = config_from_env();
    let mut stats = Stats {
        total: anchors.len(),
        ..Default::default()
    };

    let mut by_start: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    let mut by_end: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, a) in anchors.iter().enumerate() {
        by_start
            .entry((a.offset_a, a.offset_b))
            .or_default()
            .push(i);
        by_end.entry((a.limit_a, a.limit_b)).or_default().push(i);
        if a.cost == Cost::zero() {
            stats.zero_cost += 1;
        }
    }
    stats.distinct_starts = by_start.len();
    stats.distinct_ends = by_end.len();
    for group in by_start.values() {
        stats.max_start_group = stats.max_start_group.max(group.len());
        if group.len() > 1 {
            stats.in_shared_start += group.len();
        }
    }
    for group in by_end.values() {
        stats.max_end_group = stats.max_end_group.max(group.len());
        if group.len() > 1 {
            stats.in_shared_end += group.len();
        }
    }

    let mut r1_le = vec![false; anchors.len()];
    let mut r1_lt = vec![false; anchors.len()];
    let mut r2_le = vec![false; anchors.len()];
    let mut r2_lt = vec![false; anchors.len()];
    // Whether the domination keeps the descendant coordinate of the differing side.
    let mut r1_le_desc = vec![false; anchors.len()];
    let mut r1_lt_desc = vec![false; anchors.len()];
    let mut r2_le_desc = vec![false; anchors.len()];
    let mut r2_lt_desc = vec![false; anchors.len()];

    for group in by_start.values() {
        for &i1 in group {
            for &i2 in group {
                let (a1, a2) = (&anchors[i1], &anchors[i2]);
                if i1 == i2
                    || a1.limit_a > a2.limit_a
                    || a1.limit_b > a2.limit_b
                    || (a1.limit_a == a2.limit_a && a1.limit_b == a2.limit_b)
                {
                    continue;
                }
                let d = gap_affine_distance(
                    &seq_a[a1.limit_a..a2.limit_a],
                    &seq_b[a1.limit_b..a2.limit_b],
                    costs,
                );
                let c = a1.cost + d;
                let desc = a1.limit_b == a2.limit_b;
                if c <= a2.cost {
                    r1_le[i2] = true;
                    r1_le_desc[i2] |= desc;
                }
                if c < a2.cost {
                    r1_lt[i2] = true;
                    r1_lt_desc[i2] |= desc;
                }
            }
        }
    }
    for group in by_end.values() {
        for &i1 in group {
            for &i2 in group {
                let (a1, a2) = (&anchors[i1], &anchors[i2]);
                if i1 == i2
                    || a1.offset_a < a2.offset_a
                    || a1.offset_b < a2.offset_b
                    || (a1.offset_a == a2.offset_a && a1.offset_b == a2.offset_b)
                {
                    continue;
                }
                let d = gap_affine_distance(
                    &seq_a[a2.offset_a..a1.offset_a],
                    &seq_b[a2.offset_b..a1.offset_b],
                    costs,
                );
                let c = a1.cost + d;
                let desc = a1.offset_b == a2.offset_b;
                if c <= a2.cost {
                    r2_le[i2] = true;
                    r2_le_desc[i2] |= desc;
                }
                if c < a2.cost {
                    r2_lt[i2] = true;
                    r2_lt_desc[i2] |= desc;
                }
            }
        }
    }
    stats.r1_le = r1_le.iter().filter(|b| **b).count();
    stats.r1_lt = r1_lt.iter().filter(|b| **b).count();
    stats.r2_le = r2_le.iter().filter(|b| **b).count();
    stats.r2_lt = r2_lt.iter().filter(|b| **b).count();
    stats.either_le = (0..anchors.len()).filter(|&i| r1_le[i] || r2_le[i]).count();

    let (r1, r1_desc, r2, r2_desc) = if strict {
        (&r1_lt, &r1_lt_desc, &r2_lt, &r2_lt_desc)
    } else {
        (&r1_le, &r1_le_desc, &r2_le, &r2_le_desc)
    };
    let prune: Vec<bool> = (0..anchors.len())
        .map(|i| {
            let use_r1 = matches!(mode, PruneMode::Start | PruneMode::Both);
            let use_r2 = matches!(mode, PruneMode::End | PruneMode::Both);
            let (d1, d2) = if !is_secondary {
                (r1[i], r2[i])
            } else {
                match secondary_policy {
                    SecondaryPolicy::Skip => (false, false),
                    SecondaryPolicy::DescendantFixed => (r1_desc[i], r2_desc[i]),
                    SecondaryPolicy::All => (r1[i], r2[i]),
                }
            };
            (use_r1 && d1) || (use_r2 && d2)
        })
        .collect();
    stats.pruned = prune.iter().filter(|b| **b).count();

    info!(
        "ANCHOR_STATS kind={label} total={} distinct_starts={} distinct_ends={} in_shared_start={} in_shared_end={} max_start_group={} max_end_group={} zero_cost={} r1_le={} r1_lt={} r2_le={} r2_lt={} either_le={} pruned={} mode={mode:?} secondary={secondary_policy:?} strict={strict}",
        stats.total,
        stats.distinct_starts,
        stats.distinct_ends,
        stats.in_shared_start,
        stats.in_shared_end,
        stats.max_start_group,
        stats.max_end_group,
        stats.zero_cost,
        stats.r1_le,
        stats.r1_lt,
        stats.r2_le,
        stats.r2_lt,
        stats.either_le,
        stats.pruned,
    );

    anchors
        .into_iter()
        .zip(prune)
        .filter_map(|(a, p)| (!p).then_some(a))
        .collect()
}
