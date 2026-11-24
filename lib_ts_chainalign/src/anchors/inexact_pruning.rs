//! Inexact (heuristic) pruning of inexact anchors.
//!
//! The exact pruning in [`super::domination`] removes only anchors whose removal provably does not change the cost of
//! the optimal chain, which restricts it to primary anchors with `max_anchor_mutations == 1`.
//! This module implements the inexact pruning mode ([`AnchorPruning::Inexact`](crate::chain_align::performance_parameters::AnchorPruning::Inexact)),
//! which removes more anchors, also secondary anchors and anchors with `max_anchor_mutations >= 2`,
//! at the risk of increasing the cost of the optimal chain and losing template switches (see issue #236).
//!
//! Removing anchors only removes chains, so the chain cost never decreases.
//!
//! # Anchor space
//!
//! All rules work in the space in which the anchors are generated: an anchor is an alignment between a range
//! `a = [a_start, a_end)` of a sequence `A` and a range `b = [b_start, b_end)` of a sequence `B`,
//! one of which has length `k`, with at most `max_anchor_mutations` substitutions and gap characters.
//! For primary anchors, `A` is the first and `B` the second sequence.
//! For secondary anchors, `A` is the reverse complement of the ancestor and `B` is the descendant.
//! The secondary anchors of each template switch kind (ancestor and descendant) are pruned separately,
//! so a secondary anchor is only compared to anchors of the same kind.
//! There are only reverse complement template switches, so all anchors of a kind have the same orientation.
//!
//! The *diagonal* of a position `(i, j)` is `i - j`.
//! An *exact anchor* is an anchor with zero cost and both sides of length `k`.
//!
//! # Rules
//!
//! Each rule exists for primary and for secondary anchors, and can be switched on and off individually
//! (see [`InexactPruningRule`]). The names are stable, since experiments refer to them.
//!
//! * `prefix-domination` (`primary-prefix-domination`, `secondary-prefix-domination`):
//!   remove an anchor `A2` if there is an anchor `A1` with the same start such that `A2` is `A1` plus one column
//!   `c` at its end (a gap character in either sequence, a match or a substitution), and
//!   `cost(A1) + cost(c) <= cost(A2)`, where a gap column costs `gap_open`.
//!   This is the rule of the exact pruning, without its restriction to anchors with a side of length `k + 1`,
//!   so dominators may be dominated themselves (cascade).
//! * `suffix-domination`: the mirrored rule, `A1` has the same end, and `c` is the first column of `A2`.
//! * `drop-guards`: a modifier of the two domination rules.
//!   Without it, the domination rules keep the cases that the proof of the exact pruning excludes:
//!   anchors with zero cost, all anchors if `gap_extend > gap_open`,
//!   and all anchors if the inexact lower bound is not max-match-run.
//!   With it, these cases are pruned as well. Alone (without a domination rule), it removes nothing.
//! * `edge-substitution`: remove an anchor whose sides both have length `k`, whose alignment without gaps
//!   has between one and `max_anchor_mutations` substitutions and costs exactly as much as the anchor,
//!   and has a substitution in its first or last column.
//!   (For `max_anchor_mutations == 1`, these are runs of `k - 1` matches plus a substitution at one end.)
//! * `exact-overlap`: remove an anchor that is not exact, if there is an exact anchor `E` such that
//!   the anchor and `E` overlap in both sequences (share at least one position of `A` and one position of `B`),
//!   and the diagonal of `E` has a distance of at most `exact_overlap_diagonal_distance` to the diagonal range of the
//!   anchor, which is the range between the diagonal of its start and the diagonal of its end.
//!   With the default distance of 1, these are anchors on the same or a neighbouring diagonal.
//!
//! There is no guard for the length of template switches: pruning secondary anchors may shorten or remove
//! the secondary anchors of short template switches, which then become infeasible.
//!
//! # Determinism and cascades
//!
//! All rules are evaluated against the unpruned set of anchors (of the same type and template switch kind),
//! and all anchors that match at least one enabled rule are removed at once.
//! The result does not depend on the order of the anchors or of the rules.
//! A dominator is strictly shorter than the anchor it dominates, so every anchor removed by a domination rule has an
//! anchor at the end of its chain of dominators that is not dominated, which is kept unless another rule removes it.
//!
//! # Statistics
//!
//! For each rule, two numbers are reported (see [`RuleCounts`]):
//! the number of anchors *matched* by the rule (independently of the other rules, so these overlap),
//! and the number of anchors *removed* by the rule, where an anchor that matches several rules is attributed to the
//! first of them in the order `prefix-domination`, `suffix-domination`, `drop-guards`, `edge-substitution`,
//! `exact-overlap`. An anchor that is removed by a domination rule only because `drop-guards` is enabled
//! is attributed to `drop-guards`. The removed numbers sum up to the number of removed anchors.

use std::{collections::BTreeMap, fmt::Display, str::FromStr};

use clap::ValueEnum;
use generic_a_star::cost::AStarCost;
use rustc_hash::FxHashMap;

use crate::costs::GapAffineCosts;

#[cfg(test)]
mod tests;

/// The type of an anchor: primary or secondary.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum AnchorType {
    Primary,
    Secondary,
}

/// A rule of the inexact pruning, independent of the anchor type.
///
/// The order of the variants is the order in which removed anchors are attributed to rules.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub enum RuleKind {
    PrefixDomination,
    SuffixDomination,
    DropGuards,
    EdgeSubstitution,
    ExactOverlap,
}

impl RuleKind {
    pub const ALL: [Self; 5] = [
        Self::PrefixDomination,
        Self::SuffixDomination,
        Self::DropGuards,
        Self::EdgeSubstitution,
        Self::ExactOverlap,
    ];

    fn index(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::PrefixDomination => "prefix-domination",
            Self::SuffixDomination => "suffix-domination",
            Self::DropGuards => "drop-guards",
            Self::EdgeSubstitution => "edge-substitution",
            Self::ExactOverlap => "exact-overlap",
        }
    }
}

/// A rule of the inexact pruning mode, for primary or secondary anchors.
///
/// See the [module documentation](self) for the precise semantics.
/// The names (in kebab case, e.g. `primary-prefix-domination`) are stable.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, PartialOrd, Ord, ValueEnum)]
pub enum InexactPruningRule {
    /// P1: also prune with the domination rules where the exact pruning does not (zero-cost anchors,
    /// gap_extend > gap_open, other lower bounds than max-match-run). Modifies the primary domination rules.
    PrimaryDropGuards,
    /// P2: prune primary anchors that are a shorter anchor with the same start plus one column.
    PrimaryPrefixDomination,
    /// P2 (mirrored): prune primary anchors that are a shorter anchor with the same end plus one column.
    PrimarySuffixDomination,
    /// P3: prune primary substitution anchors with a substitution in their first or last column.
    PrimaryEdgeSubstitution,
    /// P4: prune inexact primary anchors that overlap an exact primary anchor on the same or a neighbouring diagonal.
    PrimaryExactOverlap,
    /// S1: like P1, for secondary anchors.
    SecondaryDropGuards,
    /// S2: like P2, for secondary anchors.
    SecondaryPrefixDomination,
    /// S2 (mirrored): like P2 (mirrored), for secondary anchors.
    SecondarySuffixDomination,
    /// S3: like P3, for secondary anchors.
    SecondaryEdgeSubstitution,
    /// S4: like P4, for secondary anchors of the same template switch kind.
    SecondaryExactOverlap,
}

impl InexactPruningRule {
    pub const ALL: [Self; 10] = [
        Self::PrimaryDropGuards,
        Self::PrimaryPrefixDomination,
        Self::PrimarySuffixDomination,
        Self::PrimaryEdgeSubstitution,
        Self::PrimaryExactOverlap,
        Self::SecondaryDropGuards,
        Self::SecondaryPrefixDomination,
        Self::SecondarySuffixDomination,
        Self::SecondaryEdgeSubstitution,
        Self::SecondaryExactOverlap,
    ];

    pub fn new(anchor_type: AnchorType, kind: RuleKind) -> Self {
        *Self::ALL
            .iter()
            .find(|rule| rule.anchor_type() == anchor_type && rule.kind() == kind)
            .unwrap()
    }

    pub fn anchor_type(self) -> AnchorType {
        match self {
            Self::PrimaryDropGuards
            | Self::PrimaryPrefixDomination
            | Self::PrimarySuffixDomination
            | Self::PrimaryEdgeSubstitution
            | Self::PrimaryExactOverlap => AnchorType::Primary,
            _ => AnchorType::Secondary,
        }
    }

    pub fn kind(self) -> RuleKind {
        match self {
            Self::PrimaryDropGuards | Self::SecondaryDropGuards => RuleKind::DropGuards,
            Self::PrimaryPrefixDomination | Self::SecondaryPrefixDomination => {
                RuleKind::PrefixDomination
            }
            Self::PrimarySuffixDomination | Self::SecondarySuffixDomination => {
                RuleKind::SuffixDomination
            }
            Self::PrimaryEdgeSubstitution | Self::SecondaryEdgeSubstitution => {
                RuleKind::EdgeSubstitution
            }
            Self::PrimaryExactOverlap | Self::SecondaryExactOverlap => RuleKind::ExactOverlap,
        }
    }

    /// The stable name of the rule, e.g. `primary-prefix-domination`.
    pub fn name(self) -> String {
        format!(
            "{}-{}",
            match self.anchor_type() {
                AnchorType::Primary => "primary",
                AnchorType::Secondary => "secondary",
            },
            self.kind().name()
        )
    }

    fn bit(self) -> u16 {
        1 << Self::ALL.iter().position(|rule| *rule == self).unwrap()
    }
}

impl Display for InexactPruningRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// A set of [`InexactPruningRule`]s.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct InexactPruningRules {
    bits: u16,
}

impl InexactPruningRules {
    pub fn all() -> Self {
        Self::from_iter(InexactPruningRule::ALL)
    }

    pub fn none() -> Self {
        Self { bits: 0 }
    }

    pub fn contains(&self, rule: InexactPruningRule) -> bool {
        self.bits & rule.bit() != 0
    }

    pub fn insert(&mut self, rule: InexactPruningRule) {
        self.bits |= rule.bit();
    }

    pub fn remove(&mut self, rule: InexactPruningRule) {
        self.bits &= !rule.bit();
    }

    pub fn is_empty(&self) -> bool {
        self.bits == 0
    }

    pub fn iter(&self) -> impl Iterator<Item = InexactPruningRule> + '_ {
        InexactPruningRule::ALL
            .into_iter()
            .filter(|rule| self.contains(*rule))
    }

    /// Returns true if any rule for the given anchor type is enabled.
    pub fn any_for(&self, anchor_type: AnchorType) -> bool {
        self.iter().any(|rule| rule.anchor_type() == anchor_type)
    }

    pub(crate) fn rule_set(&self, anchor_type: AnchorType) -> RuleSet {
        let mut enabled = [false; 5];
        for kind in RuleKind::ALL {
            enabled[kind.index()] = self.contains(InexactPruningRule::new(anchor_type, kind));
        }
        RuleSet { enabled }
    }
}

impl Default for InexactPruningRules {
    fn default() -> Self {
        Self::all()
    }
}

impl FromIterator<InexactPruningRule> for InexactPruningRules {
    fn from_iter<T: IntoIterator<Item = InexactPruningRule>>(iter: T) -> Self {
        let mut rules = Self::none();
        for rule in iter {
            rules.insert(rule);
        }
        rules
    }
}

impl Display for InexactPruningRules {
    /// Writes the rules as a comma-separated list, or `none` if the set is empty.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() {
            return write!(f, "none");
        }
        let names: Vec<_> = self.iter().map(InexactPruningRule::name).collect();
        write!(f, "{}", names.join(","))
    }
}

impl FromStr for InexactPruningRules {
    type Err = String;

    /// Parses a comma-separated list of rule names.
    ///
    /// Besides the rule names, the following keywords are accepted:
    /// `all` (all rules), `none` (no rule), `primary` (all primary rules), `secondary` (all secondary rules).
    /// A name prefixed with `-` removes the rule(s) from the set built so far, e.g. `all,-primary-exact-overlap`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut rules = Self::none();
        for item in s.split(',').map(str::trim).filter(|item| !item.is_empty()) {
            let (negate, name) = match item.strip_prefix('-') {
                Some(name) => (true, name),
                None => (false, item),
            };
            let selected: Vec<InexactPruningRule> = match name {
                "all" => InexactPruningRule::ALL.to_vec(),
                "none" => Vec::new(),
                "primary" | "secondary" => {
                    let anchor_type = if name == "primary" {
                        AnchorType::Primary
                    } else {
                        AnchorType::Secondary
                    };
                    InexactPruningRule::ALL
                        .into_iter()
                        .filter(|rule| rule.anchor_type() == anchor_type)
                        .collect()
                }
                name => vec![
                    InexactPruningRule::ALL
                        .into_iter()
                        .find(|rule| rule.name() == name)
                        .ok_or_else(|| {
                            format!(
                                "unknown inexact pruning rule '{name}', expected one of all, none, primary, secondary, {}",
                                InexactPruningRule::ALL.map(InexactPruningRule::name).join(", ")
                            )
                        })?,
                ],
            };
            for rule in selected {
                if negate {
                    rules.remove(rule);
                } else {
                    rules.insert(rule);
                }
            }
        }
        Ok(rules)
    }
}

/// Parameters of the inexact pruning mode.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct InexactPruningParameters {
    /// The enabled rules. Default: all.
    pub rules: InexactPruningRules,
    /// The maximum distance between the diagonal of an exact anchor and the diagonal range of an inexact anchor
    /// for the `exact-overlap` rules. Default: 1 (same or neighbouring diagonal).
    pub exact_overlap_diagonal_distance: usize,
}

impl Default for InexactPruningParameters {
    fn default() -> Self {
        Self {
            rules: InexactPruningRules::all(),
            exact_overlap_diagonal_distance: 1,
        }
    }
}

/// The enabled rules for one anchor type, indexed by [`RuleKind`].
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(crate) struct RuleSet {
    enabled: [bool; 5],
}

impl RuleSet {
    fn contains(&self, kind: RuleKind) -> bool {
        self.enabled[kind.index()]
    }

    fn any(&self) -> bool {
        self.enabled.iter().any(|enabled| *enabled)
    }
}

/// The number of anchors matched and removed by each rule, for one anchor type.
///
/// See the [module documentation](self#statistics).
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub struct RuleCounts {
    /// Indexed by [`RuleKind`].
    pub matched: [usize; 5],
    /// Indexed by [`RuleKind`].
    pub removed: [usize; 5],
}

impl RuleCounts {
    pub fn matched(&self, kind: RuleKind) -> usize {
        self.matched[kind.index()]
    }

    pub fn removed(&self, kind: RuleKind) -> usize {
        self.removed[kind.index()]
    }

    pub fn total_removed(&self) -> usize {
        self.removed.iter().sum()
    }

    pub fn add(&mut self, other: &Self) {
        for index in 0..self.matched.len() {
            self.matched[index] += other.matched[index];
            self.removed[index] += other.removed[index];
        }
    }

    /// Adds the counts to the given maps, keyed by the full rule names (e.g. `primary-prefix-domination`).
    pub fn insert_into(
        &self,
        anchor_type: AnchorType,
        matched: &mut BTreeMap<String, u64>,
        removed: &mut BTreeMap<String, u64>,
    ) {
        for kind in RuleKind::ALL {
            let name = InexactPruningRule::new(anchor_type, kind).name();
            *matched.entry(name.clone()).or_default() += self.matched(kind) as u64;
            *removed.entry(name).or_default() += self.removed(kind) as u64;
        }
    }
}

/// An anchor in the space in which it was generated, see the [module documentation](self#anchor-space).
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub(crate) struct RawAnchor<Cost> {
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
    pub cost: Cost,
}

impl<Cost> RawAnchor<Cost> {
    fn key(&self) -> (usize, usize, usize, usize) {
        (self.a_start, self.a_end, self.b_start, self.b_end)
    }

    fn len_a(&self) -> usize {
        self.a_end - self.a_start
    }

    fn len_b(&self) -> usize {
        self.b_end - self.b_start
    }
}

/// The input of the inexact pruning for one set of anchors (primary, or secondary of one template switch kind).
pub(crate) struct PruningContext<'a, Cost> {
    /// The sequence `A` of the anchor space.
    pub seq_a: &'a [u8],
    /// The sequence `B` of the anchor space.
    pub seq_b: &'a [u8],
    pub k: usize,
    pub max_mutations: usize,
    /// The costs with which the anchors were generated.
    pub costs: &'a GapAffineCosts<Cost>,
    /// True if the inexact lower bound is max-match-run.
    pub max_match_run_lower_bound: bool,
    pub rules: RuleSet,
    pub exact_overlap_diagonal_distance: usize,
}

/// Removes the anchors that match an enabled rule and returns the statistics.
///
/// `to_raw` converts an anchor into the [anchor space](self#anchor-space).
/// The order of the remaining anchors is preserved.
pub(crate) fn prune<Anchor, Cost: AStarCost>(
    anchors: &mut Vec<Anchor>,
    to_raw: impl Fn(&Anchor) -> RawAnchor<Cost>,
    context: &PruningContext<'_, Cost>,
) -> RuleCounts {
    let mut counts = RuleCounts::default();
    if !context.rules.any() || anchors.is_empty() {
        return counts;
    }

    let raw: Vec<_> = anchors.iter().map(&to_raw).collect();
    let remove = matching_rules(&raw, context, &mut counts);

    let mut index = 0;
    anchors.retain(|_| {
        let keep = !remove[index];
        index += 1;
        keep
    });
    counts
}

/// Returns for each anchor whether it is removed, and fills `counts`.
fn matching_rules<Cost: AStarCost>(
    anchors: &[RawAnchor<Cost>],
    context: &PruningContext<'_, Cost>,
    counts: &mut RuleCounts,
) -> Vec<bool> {
    let rules = context.rules;
    let any_domination =
        rules.contains(RuleKind::PrefixDomination) || rules.contains(RuleKind::SuffixDomination);

    let mut cost_by_range: FxHashMap<(usize, usize, usize, usize), Cost> = FxHashMap::default();
    if any_domination {
        for anchor in anchors {
            cost_by_range
                .entry(anchor.key())
                .and_modify(|cost| *cost = (*cost).min(anchor.cost))
                .or_insert(anchor.cost);
        }
    }

    let exact_index = if rules.contains(RuleKind::ExactOverlap) {
        ExactAnchorIndex::new(anchors, context.k)
    } else {
        ExactAnchorIndex::default()
    };

    // Guards that do not depend on the anchor.
    let global_guards_hold =
        context.costs.gap_extend <= context.costs.gap_open && context.max_match_run_lower_bound;

    anchors
        .iter()
        .map(|anchor| {
            let mut matched = [false; 5];

            if any_domination {
                let guards_hold = global_guards_hold && !anchor.cost.is_zero();
                let prefix = rules.contains(RuleKind::PrefixDomination)
                    && is_dominated(anchor, true, &cost_by_range, context);
                let suffix = rules.contains(RuleKind::SuffixDomination)
                    && is_dominated(anchor, false, &cost_by_range, context);
                matched[RuleKind::PrefixDomination.index()] = prefix && guards_hold;
                matched[RuleKind::SuffixDomination.index()] = suffix && guards_hold;
                matched[RuleKind::DropGuards.index()] =
                    rules.contains(RuleKind::DropGuards) && (prefix || suffix) && !guards_hold;
            }

            matched[RuleKind::EdgeSubstitution.index()] =
                rules.contains(RuleKind::EdgeSubstitution) && is_edge_substitution(anchor, context);
            matched[RuleKind::ExactOverlap.index()] =
                rules.contains(RuleKind::ExactOverlap) && exact_index.overlaps(anchor, context);

            for (index, matched) in matched.iter().enumerate() {
                if *matched {
                    counts.matched[index] += 1;
                }
            }
            if let Some(first) = matched.iter().position(|matched| *matched) {
                counts.removed[first] += 1;
                true
            } else {
                false
            }
        })
        .collect()
}

fn substitution_cost<Cost: AStarCost>(
    context: &PruningContext<'_, Cost>,
    a: usize,
    b: usize,
) -> Cost {
    if context.seq_a[a] == context.seq_b[b] {
        Cost::zero()
    } else {
        context.costs.substitution
    }
}

/// Returns true if the anchor is dominated by an anchor without its last column (`at_end`) or its first column.
///
/// This ignores the guards.
fn is_dominated<Cost: AStarCost>(
    anchor: &RawAnchor<Cost>,
    at_end: bool,
    cost_by_range: &FxHashMap<(usize, usize, usize, usize), Cost>,
    context: &PruningContext<'_, Cost>,
) -> bool {
    let (len_a, len_b) = (anchor.len_a(), anchor.len_b());
    let both_cost = if len_a == 0 || len_b == 0 {
        None
    } else if at_end {
        Some(substitution_cost(
            context,
            anchor.a_end - 1,
            anchor.b_end - 1,
        ))
    } else {
        Some(substitution_cost(context, anchor.a_start, anchor.b_start))
    };

    [
        (1, 0, Some(context.costs.gap_open)),
        (0, 1, Some(context.costs.gap_open)),
        (1, 1, both_cost),
    ]
    .into_iter()
    .any(|(delta_a, delta_b, column_cost)| {
        let Some(column_cost) = column_cost else {
            return false;
        };
        if len_a < delta_a || len_b < delta_b {
            return false;
        }
        let dominator = if at_end {
            (
                anchor.a_start,
                anchor.a_end - delta_a,
                anchor.b_start,
                anchor.b_end - delta_b,
            )
        } else {
            (
                anchor.a_start + delta_a,
                anchor.a_end,
                anchor.b_start + delta_b,
                anchor.b_end,
            )
        };
        cost_by_range
            .get(&dominator)
            .and_then(|dominator_cost| dominator_cost.checked_add(&column_cost))
            .is_some_and(|dominated_cost| dominated_cost <= anchor.cost)
    })
}

/// Returns true if the anchor is a gap-free anchor with a substitution in its first or last column,
/// see the [module documentation](self#rules).
fn is_edge_substitution<Cost: AStarCost>(
    anchor: &RawAnchor<Cost>,
    context: &PruningContext<'_, Cost>,
) -> bool {
    let length = anchor.len_a();
    if length != context.k || anchor.len_b() != context.k || length == 0 {
        return false;
    }

    let is_substitution = |offset: usize| {
        context.seq_a[anchor.a_start + offset] != context.seq_b[anchor.b_start + offset]
    };
    if !is_substitution(0) && !is_substitution(length - 1) {
        return false;
    }

    let substitutions = (0..length)
        .filter(|offset| is_substitution(*offset))
        .count();
    if substitutions > context.max_mutations {
        return false;
    }
    let mut diagonal_cost = Cost::zero();
    for _ in 0..substitutions {
        diagonal_cost += context.costs.substitution;
    }
    diagonal_cost == anchor.cost
}

fn is_exact<Cost: AStarCost>(anchor: &RawAnchor<Cost>, k: usize) -> bool {
    anchor.cost.is_zero() && anchor.len_a() == k && anchor.len_b() == k
}

/// The exact anchors by diagonal, with their sorted start positions in `A`.
#[derive(Default)]
struct ExactAnchorIndex {
    by_diagonal: FxHashMap<isize, Vec<usize>>,
}

impl ExactAnchorIndex {
    fn new<Cost: AStarCost>(anchors: &[RawAnchor<Cost>], k: usize) -> Self {
        let mut by_diagonal: FxHashMap<isize, Vec<usize>> = FxHashMap::default();
        for anchor in anchors.iter().filter(|anchor| is_exact(anchor, k)) {
            by_diagonal
                .entry(anchor.a_start as isize - anchor.b_start as isize)
                .or_default()
                .push(anchor.a_start);
        }
        for starts in by_diagonal.values_mut() {
            starts.sort_unstable();
            starts.dedup();
        }
        Self { by_diagonal }
    }

    /// Returns true if the anchor is not exact and overlaps an exact anchor, see the [module documentation](self#rules).
    fn overlaps<Cost: AStarCost>(
        &self,
        anchor: &RawAnchor<Cost>,
        context: &PruningContext<'_, Cost>,
    ) -> bool {
        let k = context.k;
        if is_exact(anchor, k) || self.by_diagonal.is_empty() {
            return false;
        }

        let start_diagonal = anchor.a_start as isize - anchor.b_start as isize;
        let end_diagonal = anchor.a_end as isize - anchor.b_end as isize;
        let distance = context.exact_overlap_diagonal_distance as isize;
        let first_diagonal = start_diagonal.min(end_diagonal) - distance;
        let last_diagonal = start_diagonal.max(end_diagonal) + distance;

        (first_diagonal..=last_diagonal).any(|diagonal| {
            let Some(starts) = self.by_diagonal.get(&diagonal) else {
                return false;
            };
            // Exact anchors that overlap the anchor in `A` start in (a_start - k, a_end).
            let first = starts.partition_point(|&start| start + k <= anchor.a_start);
            starts[first..]
                .iter()
                .take_while(|&&start| start < anchor.a_end)
                .any(|&start| {
                    let b_start = start as isize - diagonal;
                    b_start < anchor.b_end as isize
                        && b_start + k as isize > anchor.b_start as isize
                })
        })
    }
}
