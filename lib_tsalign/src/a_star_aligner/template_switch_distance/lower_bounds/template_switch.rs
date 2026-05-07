use std::{
    collections::{HashMap, HashSet},
    fmt::{Display, Write},
    iter,
};

use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{AStar, AStarNode, AStarResult, cost::AStarCost};
use log::{debug, info, trace};
use traitsequence::interface::Sequence;

use crate::{
    a_star_aligner::{
        alignment_geometry::AlignmentRange,
        alignment_result::IAlignmentType,
        template_switch_distance::{
            Context, Identifier, Node,
            context::Memory,
            identifier::GapType,
            strategies::{
                AlignmentStrategySelection, chaining::NoChainingStrategy,
                node_ord::CostOnlyNodeOrdStrategy, primary_match::AllowPrimaryMatchStrategy,
                primary_range::NoPrunePrimaryRangeStrategy,
                secondary_deletion::ForbidSecondaryDeletionStrategy, shortcut::NoShortcutStrategy,
                template_switch_count::MaxTemplateSwitchCountStrategy,
                template_switch_min_length::NoTemplateSwitchMinLengthStrategy,
                template_switch_total_length::NoTemplateSwitchTotalLengthStrategy,
            },
        },
    },
    config::TemplateSwitchConfig,
    costs::gap_affine::GapAffineAlignmentCostTable,
};

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchLowerBoundMatrix<Cost> {
    entries: Vec<TSLBMatrixEntry<Cost>>,
    min_distance_between_two_template_switches: usize,
}

#[derive(Debug, Clone, Ord, PartialOrd, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TSLBMatrixEntry<Cost> {
    x: isize,
    y: isize,
    cost: Cost,
}

type TSLBAlignmentStrategies<AlphabetType, Cost> = AlignmentStrategySelection<
    AlphabetType,
    Cost,
    CostOnlyNodeOrdStrategy,
    NoTemplateSwitchMinLengthStrategy<Cost>,
    NoChainingStrategy<Cost>,
    MaxTemplateSwitchCountStrategy,
    ForbidSecondaryDeletionStrategy,
    NoShortcutStrategy<Cost>,
    AllowPrimaryMatchStrategy,
    NoPrunePrimaryRangeStrategy,
    NoTemplateSwitchTotalLengthStrategy,
>;

impl<Cost: AStarCost> TemplateSwitchLowerBoundMatrix<Cost> {
    pub fn new<AlphabetType: Alphabet>(config: &TemplateSwitchConfig<AlphabetType, Cost>) -> Self {
        info!("Computing TS lower bound matrix...");
        let lower_bound_config = generate_template_switch_lower_bound_config(config);
        assert!(
            lower_bound_config
                .secondary_reverse_edit_costs
                .min_gap_extend_cost()
                > Cost::zero(),
            "Secondary gap extend costs must be greater than zero for all alphabet characters."
        );

        let mut open_lower_bounds = HashSet::new();
        open_lower_bounds.insert((0isize, 0isize));
        let mut closed_lower_bounds = HashMap::new();
        let mut previous_closed_lower_bounds = HashMap::new();
        let mut genome_length = 1_000;

        'outer: loop {
            debug!("Using genome length {genome_length}");
            assert!(genome_length < usize::try_from(isize::MAX).unwrap() / 2);
            let genome = VectorGenome::<AlphabetType>::from_iter(iter::repeat_n(
                AlphabetType::iter().next().unwrap(),
                genome_length,
            ));
            let mut a_star = AStar::<_>::new(Context::<
                _,
                TSLBAlignmentStrategies<AlphabetType, Cost>,
            >::new(
                genome.as_genome_subsequence(),
                genome.as_genome_subsequence(),
                "",
                "",
                AlignmentRange::new_complete(genome.len(), genome.len()),
                lower_bound_config.clone(),
                Memory {
                    template_switch_min_length: (),
                    chaining: (),
                    template_switch_count: 1,
                    shortcut: (),
                    primary_match: (),
                },
                None,
                None,
                false,
            ));
            let root_xy = genome_length / 2;
            a_star.initialise_with(|context| Node::new_root_at(root_xy, root_xy, context).into());
            previous_closed_lower_bounds.extend(closed_lower_bounds.drain());

            let root_xy_isize = isize::try_from(root_xy).unwrap();
            let string_length_isize = isize::try_from(genome_length).unwrap();

            'inner: while let Some(coordinates) = open_lower_bounds.iter().next().copied() {
                let (x, y) = open_lower_bounds.take(&coordinates).unwrap();
                if closed_lower_bounds.contains_key(&(x, y)) {
                    continue 'inner;
                }

                trace!("Searching for target ({x}, {y})");

                let has_target = match a_star.search_until(|context, node| {
                    match *node.identifier() {
                        Identifier::Primary {
                            reference_index,
                            query_index,
                            ..
                        } => {
                            reference_index == 0
                                || query_index == 0
                                || reference_index == genome_length - 1
                                || query_index == genome_length - 1
                        }
                        Identifier::PrimaryReentry {
                            reference_index,
                            query_index,
                            gap_type,
                            flank_index,
                            ..
                        } => {
                            debug_assert_eq!(gap_type, GapType::None);
                            debug_assert_eq!(flank_index, 0);

                            let reentry_x =
                                isize::try_from(reference_index).unwrap() - root_xy_isize;
                            let reentry_y = isize::try_from(query_index).unwrap() - root_xy_isize;

                            trace!("Found target {} at cost {}", node.identifier(), node.cost());

                            let previous =
                                closed_lower_bounds.insert((reentry_x, reentry_y), node.cost());

                            // If previous exists, then it was inserted when searching with a shorter string.
                            debug_assert!(previous.is_none());
                            let previous_previous = previous_closed_lower_bounds
                                .get(&(reentry_x, reentry_y))
                                .copied();
                            debug_assert!(
                                previous_previous.is_none()
                                    || previous_previous == Some(node.cost())
                            );

                            enqueue_neighbours(
                                reentry_x,
                                reentry_y,
                                &mut closed_lower_bounds,
                                &mut open_lower_bounds,
                            );

                            /*trace!(
                                "Found primary reentry node {} at cost {}",
                                node.identifier(),
                                node.cost()
                            );*/
                            reentry_x == x && reentry_y == y
                        }
                        Identifier::TemplateSwitchEntrance {
                            entrance_reference_index,
                            entrance_query_index,
                            template_switch_first_offset,
                            ..
                        } => {
                            entrance_reference_index == 0
                                || entrance_query_index == 0
                                || entrance_reference_index == genome_length - 1
                                || entrance_query_index == genome_length - 1
                                || template_switch_first_offset + root_xy_isize == 0
                                || template_switch_first_offset + root_xy_isize
                                    == string_length_isize - 1
                        }
                        Identifier::Secondary {
                            primary_index,
                            secondary_index,
                            ..
                        } => {
                            primary_index == 0
                                || primary_index == genome_length - 1
                                || secondary_index == 0
                                || secondary_index == genome_length - 1
                        }
                        Identifier::TemplateSwitchExit {
                            template_switch_direction,
                            anti_primary_gap,
                            ..
                        } => {
                            let cost_increment = context
                                .config
                                .anti_primary_gap_costs(template_switch_direction)
                                .evaluate(&anti_primary_gap);
                            node.generate_primary_reentry_successor(context, cost_increment)
                                .is_none()
                        }
                    }
                }) {
                    AStarResult::FoundTarget { identifier, cost } => {
                        trace!("Search termianted with target {identifier} at cost {cost}");

                        if let Identifier::PrimaryReentry { .. } = identifier {
                            debug_assert_eq!(closed_lower_bounds.get(&(x, y)), Some(cost).as_ref());
                            true
                        } else {
                            trace!("{:?}", {
                                let mut alignment = Vec::new();

                                // Backtrack.
                                for alignment_type in a_star.backtrack() {
                                    if let Some((count, previous_alignment_type)) =
                                        alignment.last_mut()
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
                            open_lower_bounds.insert((x, y));
                            genome_length *= 2;
                            continue 'outer;
                        }
                    }
                    AStarResult::NoTarget => {
                        trace!("Search terminated without target");
                        let previous = closed_lower_bounds.insert((x, y), Cost::max_value());
                        debug_assert!(previous.is_none());
                        false
                    }
                    AStarResult::ExceededCostLimit { .. }
                    | AStarResult::ExceededMemoryLimit { .. } => unreachable!("No limits set"),
                };

                if has_target {
                    trace!("Inserting neighbours after search");
                    enqueue_neighbours(x, y, &mut closed_lower_bounds, &mut open_lower_bounds);
                }
            }

            break;
        }

        let entries = closed_lower_bounds
            .into_iter()
            .filter_map(|((x, y), cost)| {
                if cost != Cost::max_value() {
                    Some(TSLBMatrixEntry { x, y, cost })
                } else {
                    None
                }
            })
            .collect();
        let min_distance_between_two_template_switches =
            usize::try_from(config.left_flank_length + config.right_flank_length).unwrap();

        Self {
            entries,
            min_distance_between_two_template_switches,
        }
    }

    pub fn min_distance_between_two_template_switches(&self) -> usize {
        self.min_distance_between_two_template_switches
    }

    pub fn iter(&self) -> impl Iterator<Item = &TSLBMatrixEntry<Cost>> {
        self.entries.iter()
    }
}

impl<Cost> TSLBMatrixEntry<Cost> {
    pub fn x(&self) -> isize {
        self.x
    }

    pub fn y(&self) -> isize {
        self.y
    }
}

impl<Cost: Copy> TSLBMatrixEntry<Cost> {
    pub fn cost(&self) -> Cost {
        self.cost
    }
}

fn generate_template_switch_lower_bound_config<AlphabetType: Alphabet, Cost: AStarCost>(
    config: &TemplateSwitchConfig<AlphabetType, Cost>,
) -> TemplateSwitchConfig<AlphabetType, Cost> {
    TemplateSwitchConfig {
        left_flank_length: 0,
        right_flank_length: 0,
        template_switch_min_length: config.template_switch_min_length,

        base_cost: config.base_cost.clone(),

        primary_edit_costs: GapAffineAlignmentCostTable::new_max(),
        secondary_forward_edit_costs: config
            .secondary_forward_edit_costs
            .clone()
            .into_match_agnostic_lower_bound(),
        secondary_reverse_edit_costs: config
            .secondary_reverse_edit_costs
            .clone()
            .into_match_agnostic_lower_bound(),
        left_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),
        right_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),

        // The offset only affects which part of the secondary string is being compared against, but otherwise does not change anything.
        // Hence we can ignore it for the lower bound, and simply choose its minimum.
        rq_qr_offset_costs: vec![
            (isize::MIN, Cost::max_value()),
            (0, config.rq_qr_offset_costs.min(..).unwrap()),
            (1, Cost::max_value()),
        ]
        .try_into()
        .unwrap(),
        rr_qq_offset_costs: vec![
            (isize::MIN, Cost::max_value()),
            (0, config.rr_qq_offset_costs.min(..).unwrap()),
            (1, Cost::max_value()),
        ]
        .try_into()
        .unwrap(),
        length_costs: config.length_costs.clone(),
        length_difference_costs: config.length_difference_costs.clone(),
        forward_anti_primary_gap_costs: config.forward_anti_primary_gap_costs.clone(),
        reverse_anti_primary_gap_costs: config.reverse_anti_primary_gap_costs.clone(),
    }
}

fn enqueue_neighbours<Cost>(
    x: isize,
    y: isize,
    closed_lower_bounds: &mut HashMap<(isize, isize), Cost>,
    open_lower_bounds: &mut HashSet<(isize, isize)>,
) {
    for (x, y) in [
        (x + 1, y + 1),
        (x + 1, y),
        (x + 1, y - 1),
        (x, y + 1),
        (x, y - 1),
        (x - 1, y + 1),
        (x - 1, y),
        (x - 1, y - 1),
    ] {
        if !closed_lower_bounds.contains_key(&(x, y)) {
            open_lower_bounds.insert((x, y));
        }
    }
}

impl<Cost: AStarCost> Display for TemplateSwitchLowerBoundMatrix<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (min_x, max_x, min_y, max_y) = self.entries.iter().fold(
            (isize::MAX, isize::MIN, isize::MAX, isize::MIN),
            |(min_x, max_x, min_y, max_y), &TSLBMatrixEntry { x, y, cost }| {
                if cost != Cost::max_value() {
                    (min_x.min(x), max_x.max(x), min_y.min(y), max_y.max(y))
                } else {
                    (min_x, max_x, min_y, max_y)
                }
            },
        );

        writeln!(f, "TemplateSwitchLowerBoundMatrix:")?;
        writeln!(
            f,
            "min distance: {}",
            self.min_distance_between_two_template_switches
        )?;
        writeln!(
            f,
            "x range: [{min_x}, {max_x}]; y range: [{min_y}, {max_y}]",
        )?;

        let mut buffer = String::new();
        let mut matrix: HashMap<_, _> = (min_x..=max_x)
            .flat_map(|x| (min_y..=max_y).map(move |y| ((x, y), Cost::max_value())))
            .collect();
        let mut column_widths: HashMap<_, _> = (min_x..=max_x).map(|x| (x, 1)).collect();
        for TSLBMatrixEntry { x, y, cost } in &self.entries {
            debug_assert_ne!(*cost, Cost::max_value());
            buffer.clear();
            write!(buffer, "{cost}").unwrap();

            let column_width = column_widths.get_mut(x).unwrap();
            *column_width = buffer.len().max(*column_width);
            *matrix.get_mut(&(*x, *y)).unwrap() = *cost;
        }

        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let column_width = *column_widths.get(&x).unwrap();
                buffer.clear();
                let cost = *matrix.get(&(x, y)).unwrap();
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
