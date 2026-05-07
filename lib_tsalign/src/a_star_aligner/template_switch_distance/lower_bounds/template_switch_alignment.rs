use std::{
    collections::HashMap,
    fmt::{Display, Write},
    iter,
};

use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{AStar, AStarNode, AStarResult, cost::AStarCost};
use log::{debug, info, trace};
use ndarray::Array2;
use traitsequence::interface::Sequence;

use crate::{
    a_star_aligner::{
        alignment_geometry::AlignmentRange,
        alignment_result::IAlignmentType,
        template_switch_distance::{
            Context, Identifier,
            context::{DynamicStrategies, Memory},
            strategies::{
                AlignmentStrategySelection,
                chaining::NoChainingStrategy,
                descendant::AnyTemplateSwitchDescendantStrategy,
                node_ord::CostOnlyNodeOrdStrategy,
                primary_match::{
                    MaxConsecutivePrimaryMatchMemory, MaxConsecutivePrimaryMatchStrategy,
                },
                primary_range::NoPrunePrimaryRangeStrategy,
                secondary_deletion::AllowSecondaryDeletionStrategy,
                shortcut::TemplateSwitchLowerBoundShortcutStrategy,
                template_switch_count::NoTemplateSwitchCountStrategy,
                template_switch_min_length::NoTemplateSwitchMinLengthStrategy,
                template_switch_total_length::NoTemplateSwitchTotalLengthStrategy,
            },
        },
    },
    config::{BaseCost, TemplateSwitchConfig},
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
};

use super::template_switch::TemplateSwitchLowerBoundMatrix;

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchAlignmentLowerBoundMatrix<Cost> {
    matrix: Array2<Cost>,
}

type TSALBAlignmentStrategies<AlphabetType, Cost> = AlignmentStrategySelection<
    AlphabetType,
    Cost,
    CostOnlyNodeOrdStrategy,
    NoTemplateSwitchMinLengthStrategy<Cost>,
    NoChainingStrategy<Cost>,
    NoTemplateSwitchCountStrategy,
    AllowSecondaryDeletionStrategy,
    TemplateSwitchLowerBoundShortcutStrategy<Cost>,
    MaxConsecutivePrimaryMatchStrategy,
    NoPrunePrimaryRangeStrategy,
    NoTemplateSwitchTotalLengthStrategy,
    AnyTemplateSwitchDescendantStrategy,
>;

impl<Cost: AStarCost> TemplateSwitchAlignmentLowerBoundMatrix<Cost> {
    pub fn new<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
        tslb_matrix: &TemplateSwitchLowerBoundMatrix<Cost>,
        reference_length: usize,
        query_length: usize,
        max_consecutive_primary_matches: usize,
        max_consecutive_primary_matches_at_start_and_end: usize,
    ) -> Self {
        info!("Computing TS alignment lower bound matrix...");
        let lower_bound_config = generate_template_switch_alignment_lower_bound_config(config);

        let mut closed_lower_bounds = HashMap::new();
        let genome_length = reference_length.max(query_length);
        let target_min_available_primary_matches =
            max_consecutive_primary_matches - max_consecutive_primary_matches_at_start_and_end;

        debug!("Max consecutive primary matches: {max_consecutive_primary_matches}");
        debug!(
            "Max consecutive primary matches at start and end: {max_consecutive_primary_matches_at_start_and_end}"
        );
        debug!("Min available primary matches at target: {target_min_available_primary_matches}");

        debug!("Using genome length {genome_length}");
        assert!(genome_length < usize::try_from(isize::MAX).unwrap() / 2);
        let genome = VectorGenome::<AlphabetType>::from_iter(iter::repeat_n(
            AlphabetType::iter().next().unwrap(),
            genome_length,
        ));
        let mut a_star = AStar::<_>::new(
            Context::<_, TSALBAlignmentStrategies<AlphabetType, Cost>>::new(
                genome.as_genome_subsequence(),
                genome.as_genome_subsequence(),
                "",
                "",
                AlignmentRange::new_complete(genome.len(), genome.len()),
                lower_bound_config.clone(),
                Memory {
                    template_switch_min_length: (),
                    chaining: (),
                    template_switch_count: (),
                    shortcut: tslb_matrix.clone(),
                    primary_match: MaxConsecutivePrimaryMatchMemory {
                        max_consecutive_primary_matches,
                        root_available_primary_matches:
                            max_consecutive_primary_matches_at_start_and_end,
                        fake_substitution_cost: lower_bound_config
                            .primary_edit_costs
                            .min_substitution_cost(),
                    },
                },
                DynamicStrategies {},
                None,
                None,
                false,
            ),
        );
        a_star.initialise();

        for (target_reference_index, target_query_index) in
            (0..=reference_length).flat_map(|reference_index| {
                (0..=query_length).map(move |query_index| (reference_index, query_index))
            })
        {
            if closed_lower_bounds.contains_key(&(target_reference_index, target_query_index)) {
                continue;
            }

            trace!("Searching for target ({target_reference_index}, {target_query_index})");

            match a_star.search_until(|_, node| match *node.identifier() {
                Identifier::Primary {
                    reference_index,
                    query_index,
                    flank_index,
                    data,
                    ..
                }
                | Identifier::PrimaryReentry {
                    reference_index,
                    query_index,
                    flank_index,
                    data,
                    ..
                } => {
                    if flank_index == 0
                        && reference_index <= reference_length
                        && query_index <= query_length
                        && (data >= target_min_available_primary_matches || (reference_index == 0 && query_index == 0))
                    {
                        if let Some(previous) = closed_lower_bounds.get(&(reference_index, query_index)) {
                            debug_assert!(*previous <= node.cost(), "Search may find the same node thrice due to gap types, but never with a lower cost.");
                        } else {
                            closed_lower_bounds.insert((reference_index, query_index), node.cost());

                            // Found target for the first time.
                            trace!("Found target {} at cost {}", node.identifier(), node.cost());
                        }

                        reference_index == target_reference_index && query_index == target_query_index
                    } else {
                        false
                    }
                }
                Identifier::TemplateSwitchEntrance { .. }
                | Identifier::Secondary { .. }
                | Identifier::TemplateSwitchExit { .. } => {
                    unreachable!()
                }
            }) {
                AStarResult::FoundTarget { identifier, cost } => {
                    trace!("Search terminated with target {identifier} at cost {cost}");

                    if let Identifier::Primary { .. } | Identifier::PrimaryReentry { .. } =
                        identifier
                    {
                        debug_assert_eq!(
                            closed_lower_bounds.get(&(target_reference_index, target_query_index)),
                            Some(cost).as_ref()
                        );
                    } else {
                        trace!("{:?}", {
                            let mut alignment = Vec::new();

                            // Backtrack.
                            for alignment_type in
                                a_star.backtrack()
                            {
                                if let Some((count, previous_alignment_type)) = alignment.last_mut()
                                {
                                    if alignment_type.is_repeated(previous_alignment_type) {
                                        *count += 1;
                                    } else {
                                        alignment.push((1, alignment_type));
                                    }
                                } else {
                                    alignment.push((1, alignment_type));
                                }
                            }

                            alignment.reverse();
                            alignment
                        });
                    }
                }
                AStarResult::NoTarget => {
                    unreachable!("Search terminated without target for target reference index {target_reference_index} and target query index {target_query_index}");
                }
                AStarResult::ExceededCostLimit { .. }
                | AStarResult::ExceededMemoryLimit { .. } => unreachable!("No limits set"),
            }
        }

        let (max_x, max_y) = closed_lower_bounds
            .iter()
            .fold((usize::MIN, usize::MIN), |(max_x, max_y), (&(x, y), _)| {
                (max_x.max(x), max_y.max(y))
            });

        let mut matrix = Array2::from_elem((max_x + 1, max_y + 1), Cost::max_value());
        for ((x, y), cost) in matrix.indexed_iter_mut() {
            *cost = *closed_lower_bounds
                .get(&(x, y))
                .unwrap_or_else(|| unreachable!("Missing matrix entry for ({x}, {y})"));
            assert_ne!(*cost, Cost::max_value());
        }

        Self { matrix }
    }

    pub fn cost(&self, delta_reference: usize, delta_query: usize) -> Cost {
        self.matrix[(delta_reference, delta_query)]
    }
}

fn generate_template_switch_alignment_lower_bound_config<
    AlphabetType: Alphabet,
    Cost: AStarCost,
>(
    config: &TemplateSwitchConfig<AlphabetType, Cost>,
) -> TemplateSwitchConfig<AlphabetType, Cost> {
    // Template switches are handled by the TS lower bound matrix, hence we disable them in the search.
    TemplateSwitchConfig {
        left_flank_length: config.left_flank_length,
        right_flank_length: config.right_flank_length,
        template_switch_min_length: usize::MAX,

        base_cost: BaseCost::new_max(),

        primary_edit_costs: config.primary_edit_costs.clone(),
        secondary_forward_edit_costs: GapAffineAlignmentCostTable::new_max(),
        secondary_reverse_edit_costs: GapAffineAlignmentCostTable::new_max(),
        left_flank_edit_costs: config.left_flank_edit_costs.clone(),
        right_flank_edit_costs: config.right_flank_edit_costs.clone(),

        rq_qr_offset_costs: CostFunction::new_max(),
        rr_qq_offset_costs: CostFunction::new_max(),
        length_costs: CostFunction::new_max(),
        length_difference_costs: CostFunction::new_max(),
        forward_anti_primary_gap_costs: CostFunction::new_max(),
        reverse_anti_primary_gap_costs: CostFunction::new_max(),
    }
}

impl<Cost: AStarCost> Display for TemplateSwitchAlignmentLowerBoundMatrix<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let min_x = 0;
        let min_y = 0;
        let max_x = self.matrix.shape()[0] - 1;
        let max_y = self.matrix.shape()[1] - 1;

        writeln!(f, "TemplateSwitchAlignmentLowerBoundMatrix:")?;
        writeln!(f, "x range: [0, {max_x}]; y range: [0, {max_y}]",)?;

        let mut buffer = String::new();
        let mut column_widths: HashMap<_, _> = (min_x..=max_x).map(|x| (x, 1)).collect();
        for ((x, _), cost) in self.matrix.indexed_iter() {
            let current_column_width = if *cost == Cost::max_value() {
                1
            } else {
                buffer.clear();
                write!(buffer, "{cost}").unwrap();
                buffer.len()
            };

            let column_width = column_widths.get_mut(&x).unwrap();
            *column_width = current_column_width.max(*column_width);
        }

        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let column_width = *column_widths.get(&x).unwrap();
                buffer.clear();
                let cost = self.matrix[(x, y)];
                if cost == Cost::max_value() {
                    write!(buffer, "∞").unwrap();
                } else {
                    write!(buffer, "{cost}").unwrap();
                }
                for _ in 0..=column_width - buffer.chars().count() {
                    write!(f, " ")?;
                }
                write!(f, "{buffer}")?;
            }

            writeln!(f)?;
        }

        Ok(())
    }
}
