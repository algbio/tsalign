use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};
use generic_a_star::cost::AStarCost;
use log::error;

use crate::{
    a_star_aligner::{
        alignment_result::{IAlignmentType, alignment::stream::AlignmentStream},
        template_switch_distance::{
            AlignmentType, TemplateSwitchAncestor, TemplateSwitchDescendant,
            TemplateSwitchDirection,
        },
    },
    config::TemplateSwitchConfig,
};

use super::Alignment;

impl Alignment<AlignmentType> {
    /// Moves the start of a template switch one character pair to the left if possible, moving the character pair into the template switch.
    ///
    /// Only works if the alignment preceding the template switch is a match or substitution, and not a flank.
    ///
    /// In case of success, it returns true, and otherwise false.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    #[allow(clippy::too_many_lines)]
    pub fn move_template_switch_start_backwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: &mut usize,
    ) -> bool {
        let Some(&(
            _,
            AlignmentType::TemplateSwitchEntrance {
                first_offset,
                descendant,
                ancestor,
                direction,
                ..
            },
        )) = self.alignment.get(*compact_index)
        else {
            error!(
                "(Programming) Error: compact index {compact_index} does not point to a template switch entrance"
            );
            return false;
        };

        if *compact_index > 0
            && matches!(
                self.alignment.get(*compact_index - 1),
                Some((
                    _,
                    AlignmentType::PrimaryMatch | AlignmentType::PrimarySubstitution
                ))
            )
        {
            // Compute TS inner first indices.
            let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
            stream.push_all(self.iter_compact_cloned().take(*compact_index));
            let ts_inner_descendant_index = match descendant {
                TemplateSwitchDescendant::Reference => stream.head_coordinates().reference(),
                TemplateSwitchDescendant::Query => stream.head_coordinates().query(),
            };
            if ts_inner_descendant_index == 0 {
                // We cannot extend more backwards.
                return false;
            }
            let Some(ts_inner_ancestor_index) = isize::try_from(match ancestor {
                TemplateSwitchAncestor::Reference => stream.head_coordinates().reference(),
                TemplateSwitchAncestor::Query => stream.head_coordinates().query(),
            })
            .ok()
            .and_then(|i| usize::try_from(i.checked_add(first_offset)?).ok()) else {
                error!(
                    "(Programming) Error: finding inner ancestor index -- integer math does not check out"
                );
                return false;
            };

            // Check if indices can be moved while staying in bounds.
            match direction {
                TemplateSwitchDirection::Forward if ts_inner_ancestor_index == 0 => return false,
                TemplateSwitchDirection::Reverse
                    if ts_inner_ancestor_index >= ancestor.get(reference, query).len() =>
                {
                    return false;
                }
                _ => {}
            }

            // Remove one match or substitution from before the TS.
            let multiplicity = &mut self.alignment[*compact_index - 1].0;
            if *multiplicity == 0 {
                error!("Invalid input alignment! Cannot have multiplicity of 0.");
                return false;
            }
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                *compact_index -= 1;
                self.alignment.remove(*compact_index);
            }

            // Check if the new inner pair is a match or a substitution.
            // ts_inner_descendant_index > 0, otherwise we would've returned false earlier
            let descendant_char =
                descendant.get(reference, query)[ts_inner_descendant_index - 1].clone();
            let ancestor_char = match direction {
                TemplateSwitchDirection::Forward => {
                    // ts_inner_ancestor_index > 0, otherwise we would've returned false earlier
                    ancestor.get(reference, query)[ts_inner_ancestor_index - 1].clone()
                }
                TemplateSwitchDirection::Reverse => {
                    ancestor.get(reference, query)[ts_inner_ancestor_index].complement()
                }
            };
            let inner_alignment_type = if descendant_char == ancestor_char {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            };

            // Insert new inner alignment.
            if self
                .alignment
                .get(*compact_index + 1)
                .is_some_and(|(_, ty)| *ty == inner_alignment_type)
            {
                self.alignment[*compact_index + 1].0 += 1;
            } else {
                self.alignment
                    .insert(*compact_index + 1, (1, inner_alignment_type));
            }

            // If reverse TS, then fix first offset.
            // (If forward TS, the changes to points 1 and 2 cancel out.)
            if direction == TemplateSwitchDirection::Reverse {
                let AlignmentType::TemplateSwitchEntrance { first_offset, .. } =
                    &mut self.alignment[*compact_index].1
                else {
                    unreachable!("reborrow of already borrowed element, definitely exists");
                };
                *first_offset += 2;
            }

            // Fix anti-descendant gap.
            let Some((
                _,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap,
                },
            )) = self.alignment[*compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!("There should be a TS exit..");
            };
            *anti_descendant_gap += 1;

            true
        } else {
            false
        }
    }

    /// Moves the start of a template switch one character pair to the right if possible, moving the character pair out of the template switch.
    ///
    /// Only works if the first alignment inside the template switch is a match or substitution.
    ///
    /// In case of success, it returns true, and otherwise false.
    ///
    /// Assumes that the template switch is not preceded by a flank.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_start_forwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: &mut usize,
    ) -> bool {
        let AlignmentType::TemplateSwitchEntrance { direction, .. } =
            self.alignment[*compact_index].1
        else {
            error!(
                "Error: compact index {compact_index} does not point to a template switch entrance"
            );
            return false;
        };

        // Assert that no flanks are involved.
        if *compact_index != 0
            && !self
                .alignment
                .get(*compact_index - 1)
                .is_none_or(|(_, alignment_type)| {
                    !matches!(
                        alignment_type,
                        AlignmentType::PrimaryFlankDeletion
                            | AlignmentType::PrimaryFlankInsertion
                            | AlignmentType::PrimaryFlankSubstitution
                            | AlignmentType::PrimaryFlankMatch
                    )
                })
        {
            error!("We did not want flanks to get involved...");
            return false;
        }

        if let Some((_, AlignmentType::SecondaryMatch | AlignmentType::SecondarySubstitution)) =
            self.alignment.get(*compact_index + 1)
        {
            // Compute TS outer first indices.
            let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
            stream.push_all(self.iter_compact_cloned().take(*compact_index));
            let ts_outer_reference_index = stream.head_coordinates().reference();
            let ts_outer_query_index = stream.head_coordinates().query();

            // Check if indices can be moved while staying in bounds.
            if ts_outer_reference_index == reference.len() || ts_outer_query_index == query.len() {
                return false;
            }

            // Remove one match or substitution from inside the TS.
            // compact_index + 1 is valid since the index points to the TS entrance, and that is certainly not the last entry
            let multiplicity = &mut self.alignment[*compact_index + 1].0;
            if *multiplicity == 0 {
                error!("Invalid input alignment! Cannot have multiplicity of 0.");
                return false;
            }
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                self.alignment.remove(*compact_index + 1);
            }

            // Check if the new outer pair is a match or a substitution.
            let reference_char = reference[ts_outer_reference_index].clone();
            let query_char = query[ts_outer_query_index].clone();
            let outer_alignment_type = if reference_char == query_char {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert new outer alignment.
            if *compact_index != 0
                && self.alignment.get(*compact_index - 1).map(|t| t.1) == Some(outer_alignment_type)
            {
                self.alignment[*compact_index - 1].0 += 1;
            } else {
                self.alignment
                    .insert(*compact_index, (1, outer_alignment_type));
                *compact_index += 1;
            }

            // If reverse TS, then fix first offset.
            // (If forward TS, the changes to points 1 and 2 cancel out.)
            if direction == TemplateSwitchDirection::Reverse {
                let AlignmentType::TemplateSwitchEntrance { first_offset, .. } =
                    &mut self.alignment[*compact_index].1
                else {
                    unreachable!("Merely a reborrow");
                };
                *first_offset -= 2;
            }

            // Fix anti-descendant gap.
            let Some((
                _,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap,
                },
            )) = self.alignment[*compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!();
            };
            *anti_descendant_gap -= 1;

            true
        } else {
            false
        }
    }

    /// Moves the end of a template switch one character pair to the right if possible, moving the character pair into the template switch.
    ///
    /// Only works if the alignment preceding the template switch is a match or substitution, and not a flank.
    ///
    /// In case of success, it returns true, and otherwise false.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    #[allow(clippy::too_many_lines)]
    pub fn move_template_switch_end_forwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: usize,
    ) -> bool {
        let AlignmentType::TemplateSwitchEntrance {
            first_offset,
            descendant,
            ancestor,
            direction,
            ..
        } = self.alignment[compact_index].1
        else {
            error!(
                "Error: compact index {compact_index} does not point to a template switch entrance"
            );
            return false;
        };
        let Some((exit_index_offset, _)) =
            self.alignment.iter().skip(compact_index).enumerate().find(
                |(_, (_, alignment_type))| {
                    matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. })
                },
            )
        else {
            error!("There should be a TS exit after the TS entrance");
            return false;
        };
        let mut exit_index = compact_index + exit_index_offset;
        let inner_ancestor_length = self
            .alignment
            .iter()
            .take(exit_index)
            .skip(compact_index + 1)
            .fold(
                0,
                |ancestor_length, (multiplicity, alignment_type)| match alignment_type {
                    AlignmentType::SecondaryInsertion => ancestor_length,
                    AlignmentType::SecondaryDeletion
                    | AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch => ancestor_length + *multiplicity,
                    _ => unreachable!(),
                },
            );

        if let Some((_, AlignmentType::PrimaryMatch | AlignmentType::PrimarySubstitution)) =
            self.alignment.get(exit_index + 1)
        {
            // Compute TS inner last indices.
            let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
            stream.push_all(self.iter_compact_cloned().take(compact_index));
            stream.clear();
            stream.push_all(
                self.iter_compact_cloned()
                    .take(exit_index + 1)
                    .skip(compact_index),
            );
            let ts_inner_descendant_index = match descendant {
                TemplateSwitchDescendant::Reference => stream.head_coordinates().reference(),
                TemplateSwitchDescendant::Query => stream.head_coordinates().query(),
            };

            let Some(ts_inner_ancestor_index) = isize::try_from(match ancestor {
                TemplateSwitchAncestor::Reference => stream.tail_coordinates().reference(),
                TemplateSwitchAncestor::Query => stream.tail_coordinates().query(),
            })
            .ok()
            .and_then(|i| usize::try_from(i.checked_add(first_offset)?).ok()) else {
                error!("Error finding inner ancestor index -- integer math does not check out");
                return false;
            };
            let ts_inner_ancestor_index = match direction {
                TemplateSwitchDirection::Forward => {
                    let Some(ts_inner_ancestor_index) =
                        ts_inner_ancestor_index.checked_add(inner_ancestor_length)
                    else {
                        error!(
                            "ts_inner_ancestor_index {ts_inner_ancestor_index} should be at least {inner_ancestor_length} away from any boundary"
                        );
                        return false;
                    };
                    // Check if indices can be moved while staying in bounds.
                    if ts_inner_ancestor_index >= ancestor.get(reference, query).len() {
                        return false;
                    }
                    ts_inner_ancestor_index
                }
                TemplateSwitchDirection::Reverse => {
                    let Some(ts_inner_ancestor_index) =
                        ts_inner_ancestor_index.checked_sub(inner_ancestor_length)
                    else {
                        error!(
                            "ts_inner_ancestor_index {ts_inner_ancestor_index} should be at least {inner_ancestor_length} away from any boundary"
                        );
                        return false;
                    };
                    // Check if indices can be moved while staying in bounds.
                    if ts_inner_ancestor_index == 0 {
                        return false;
                    }
                    ts_inner_ancestor_index
                }
            };

            // Remove one match or substitution from after the TS.
            let multiplicity = &mut self.alignment[exit_index + 1].0;
            if *multiplicity == 0 {
                error!("Invalid input alignment! Cannot have multiplicity of 0.");
                return false;
            }
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                self.alignment.remove(exit_index + 1);
            }

            // Check if the new inner pair is a match or a substitution.
            let descendant_char =
                descendant.get(reference, query)[ts_inner_descendant_index].clone();
            let ancestor_char = match direction {
                TemplateSwitchDirection::Forward => {
                    ancestor.get(reference, query)[ts_inner_ancestor_index].clone()
                }
                TemplateSwitchDirection::Reverse => {
                    // ts_inner_ancestor_index > 0 since we checked that earlier and would have returned otherwise.
                    ancestor.get(reference, query)[ts_inner_ancestor_index - 1].complement()
                }
            };
            let inner_alignment_type = if descendant_char == ancestor_char {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            };

            // Insert new inner alignment.
            if self.alignment[exit_index - 1].1 == inner_alignment_type {
                self.alignment[exit_index - 1].0 += 1;
            } else {
                self.alignment.insert(exit_index, (1, inner_alignment_type));
                exit_index += 1;
            }

            // Fix anti-descendant gap.
            let AlignmentType::TemplateSwitchExit {
                anti_descendant_gap,
            } = &mut self.alignment[exit_index].1
            else {
                unreachable!("Merely a reborrow");
            };
            *anti_descendant_gap += 1;

            true
        } else {
            false
        }
    }

    /// Moves the end of a template switch one character pair to the left if possible, moving the character pair out of the template switch.
    ///
    /// Only works if the last alignment inside the template switch is a match or substitution.
    ///
    /// In case of success, it returns true, and otherwise false.
    ///
    /// Assumes that the template switch is not succeeded by a flank.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_end_backwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: usize,
    ) -> bool {
        if !matches!(
            self.alignment[compact_index].1,
            AlignmentType::TemplateSwitchEntrance { .. }
        ) {
            error!(
                "Error: compact index {compact_index} does not point to a template switch entrance"
            );
            return false;
        }
        let Some((exit_index_offset, _)) =
            self.alignment.iter().skip(compact_index).enumerate().find(
                |(_, (_, alignment_type))| {
                    matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. })
                },
            )
        else {
            error!("There must be a TS exit after the TS entrance!");
            return false;
        };
        let mut exit_index = compact_index + exit_index_offset;

        // Assert that no flanks are involved.
        if !self
            .alignment
            .get(exit_index + 1)
            .is_none_or(|(_, alignment_type)| {
                !matches!(
                    alignment_type,
                    AlignmentType::PrimaryFlankDeletion
                        | AlignmentType::PrimaryFlankInsertion
                        | AlignmentType::PrimaryFlankSubstitution
                        | AlignmentType::PrimaryFlankMatch
                )
            })
        {
            error!("We have some unexpected flanks!");
            return false;
        }

        if let Some((_, AlignmentType::SecondaryMatch | AlignmentType::SecondarySubstitution)) =
            self.alignment.get(exit_index - 1)
        {
            // Compute TS outer first indices.
            let mut stream = AlignmentStream::new_with_offset(reference_offset, query_offset);
            stream.push_all(self.iter_compact_cloned().take(exit_index + 1));
            let ts_outer_reference_index = stream.head_coordinates().reference();
            let ts_outer_query_index = stream.head_coordinates().query();

            // Check if indices can be moved while staying in bounds.
            if ts_outer_reference_index == 0 || ts_outer_query_index == 0 {
                return false;
            }

            // Remove one match or substitution from inside the TS.
            let multiplicity = &mut self.alignment[exit_index - 1].0;
            if *multiplicity == 0 {
                error!("Invalid input alignment! Cannot have multiplicity of 0.");
                return false;
            }
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                exit_index -= 1;
                self.alignment.remove(exit_index);
            }

            // Check if the new outer pair is a match or a substitution.
            let reference_char = reference[ts_outer_reference_index - 1].clone();
            let query_char = query[ts_outer_query_index - 1].clone();
            let outer_alignment_type = if reference_char == query_char {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert new outer alignment.
            if self.alignment.get(exit_index + 1).map(|t| t.1) == Some(outer_alignment_type) {
                self.alignment[exit_index + 1].0 += 1;
            } else {
                self.alignment
                    .insert(exit_index + 1, (1, outer_alignment_type));
            }

            // Fix anti-descendant gap.
            let Some((
                _,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap,
                },
            )) = self.alignment[compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!("Just a reborrow");
            };
            *anti_descendant_gap -= 1;

            true
        } else {
            false
        }
    }

    /// Compute the cost of an alignment.
    ///
    /// Flanks are not supported.
    pub fn compute_cost<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        Cost: AStarCost,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
    ) -> Cost {
        let mut cost = Cost::zero();

        let mut last_alignment_type = None;
        let mut reference_index = reference_offset;
        let mut query_index = query_offset;
        let mut descendant_index = 0;
        let mut ancestor_index = 0;
        let mut descendant = TemplateSwitchDescendant::Reference;
        let mut ancestor = TemplateSwitchAncestor::Reference;
        let mut direction = TemplateSwitchDirection::Forward;
        for alignment_type in self.iter_flat_cloned() {
            let cost_increment = match alignment_type {
                AlignmentType::PrimaryInsertion => {
                    let cost_increment = if Some(alignment_type) == last_alignment_type {
                        config
                            .primary_edit_costs
                            .gap_extend_cost(query[query_index].clone())
                    } else {
                        config
                            .primary_edit_costs
                            .gap_open_cost(query[query_index].clone())
                    };
                    query_index += 1;
                    cost_increment
                }
                AlignmentType::PrimaryDeletion => {
                    let cost_increment = if Some(alignment_type) == last_alignment_type {
                        config
                            .primary_edit_costs
                            .gap_extend_cost(reference[reference_index].clone())
                    } else {
                        config
                            .primary_edit_costs
                            .gap_open_cost(reference[reference_index].clone())
                    };
                    reference_index += 1;
                    cost_increment
                }
                AlignmentType::PrimarySubstitution | AlignmentType::PrimaryMatch => {
                    let cost_increment = config.primary_edit_costs.match_or_substitution_cost(
                        reference[reference_index].clone(),
                        query[query_index].clone(),
                    );
                    reference_index += 1;
                    query_index += 1;
                    cost_increment
                }
                AlignmentType::PrimaryFlankInsertion
                | AlignmentType::PrimaryFlankDeletion
                | AlignmentType::PrimaryFlankSubstitution
                | AlignmentType::PrimaryFlankMatch => {
                    todo!("Flanks are not yet supported")
                }
                AlignmentType::SecondaryInsertion => {
                    let descendant_character = match descendant {
                        TemplateSwitchDescendant::Reference => reference[descendant_index].clone(),
                        TemplateSwitchDescendant::Query => query[descendant_index].clone(),
                    };
                    let cost_increment = if Some(alignment_type) == last_alignment_type {
                        config
                            .secondary_edit_costs(direction)
                            .gap_extend_cost(descendant_character)
                    } else {
                        config
                            .secondary_edit_costs(direction)
                            .gap_open_cost(descendant_character)
                    };
                    descendant_index += 1;
                    cost_increment
                }
                AlignmentType::SecondaryDeletion => {
                    let ancestor_character = match (ancestor, direction) {
                        (TemplateSwitchAncestor::Reference, TemplateSwitchDirection::Forward) => {
                            reference[ancestor_index].clone()
                        }
                        (TemplateSwitchAncestor::Reference, TemplateSwitchDirection::Reverse) => {
                            reference[ancestor_index - 1].complement()
                        }
                        (TemplateSwitchAncestor::Query, TemplateSwitchDirection::Forward) => {
                            query[ancestor_index].clone()
                        }
                        (TemplateSwitchAncestor::Query, TemplateSwitchDirection::Reverse) => {
                            query[ancestor_index - 1].complement()
                        }
                    };
                    let cost_increment = if Some(alignment_type) == last_alignment_type {
                        config
                            .secondary_edit_costs(direction)
                            .gap_extend_cost(ancestor_character)
                    } else {
                        config
                            .secondary_edit_costs(direction)
                            .gap_open_cost(ancestor_character)
                    };
                    match direction {
                        TemplateSwitchDirection::Forward => ancestor_index += 1,
                        TemplateSwitchDirection::Reverse => ancestor_index -= 1,
                    }
                    cost_increment
                }
                AlignmentType::SecondarySubstitution | AlignmentType::SecondaryMatch => {
                    let descendant_character = match descendant {
                        TemplateSwitchDescendant::Reference => reference[descendant_index].clone(),
                        TemplateSwitchDescendant::Query => query[descendant_index].clone(),
                    };
                    let ancestor_character = match (ancestor, direction) {
                        (TemplateSwitchAncestor::Reference, TemplateSwitchDirection::Forward) => {
                            reference[ancestor_index].clone()
                        }
                        (TemplateSwitchAncestor::Reference, TemplateSwitchDirection::Reverse) => {
                            reference[ancestor_index - 1].complement()
                        }
                        (TemplateSwitchAncestor::Query, TemplateSwitchDirection::Forward) => {
                            query[ancestor_index].clone()
                        }
                        (TemplateSwitchAncestor::Query, TemplateSwitchDirection::Reverse) => {
                            query[ancestor_index - 1].complement()
                        }
                    };
                    let cost_increment = config
                        .secondary_edit_costs(direction)
                        .match_or_substitution_cost(descendant_character, ancestor_character);
                    descendant_index += 1;
                    match direction {
                        TemplateSwitchDirection::Forward => ancestor_index += 1,
                        TemplateSwitchDirection::Reverse => ancestor_index -= 1,
                    }
                    cost_increment
                }
                AlignmentType::TemplateSwitchEntrance {
                    first_offset,
                    descendant: ts_descendant,
                    ancestor: ts_ancestor,
                    direction: ts_direction,
                    ..
                } => {
                    assert!(!matches!(
                        last_alignment_type,
                        Some(AlignmentType::TemplateSwitchEntrance { .. })
                    ));
                    descendant = ts_descendant;
                    ancestor = ts_ancestor;
                    direction = ts_direction;
                    let base_cost = config.base_cost.get(descendant, ancestor, direction);
                    let offset_cost = config
                        .offset_costs(descendant, ancestor)
                        .evaluate(&first_offset);
                    let Some(cost_increment) = base_cost.checked_add(&offset_cost) else {
                        return Cost::max_value();
                    };
                    descendant_index = match descendant {
                        TemplateSwitchDescendant::Reference => reference_index,
                        TemplateSwitchDescendant::Query => query_index,
                    };
                    ancestor_index = usize::try_from(
                        isize::try_from(match ancestor {
                            TemplateSwitchAncestor::Reference => reference_index,
                            TemplateSwitchAncestor::Query => query_index,
                        })
                        .unwrap()
                        .checked_add(first_offset)
                        .unwrap(),
                    )
                    .unwrap();
                    cost_increment
                }
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap,
                } => {
                    assert!(!matches!(
                        last_alignment_type,
                        Some(AlignmentType::TemplateSwitchExit { .. })
                    ));
                    let length = match descendant {
                        TemplateSwitchDescendant::Reference => {
                            let length = descendant_index - reference_index;
                            reference_index = descendant_index;
                            query_index = usize::try_from(
                                isize::try_from(query_index)
                                    .unwrap()
                                    .checked_add(anti_descendant_gap)
                                    .unwrap(),
                            )
                            .unwrap();
                            length
                        }
                        TemplateSwitchDescendant::Query => {
                            let length = descendant_index - query_index;
                            query_index = descendant_index;
                            reference_index = usize::try_from(
                                isize::try_from(reference_index)
                                    .unwrap()
                                    .checked_add(anti_descendant_gap)
                                    .unwrap(),
                            )
                            .unwrap();
                            length
                        }
                    };
                    let length_difference = anti_descendant_gap - isize::try_from(length).unwrap();
                    let cost_increment = config
                        .anti_descendant_gap_costs(direction)
                        .evaluate(&anti_descendant_gap);
                    let Some(cost_increment) =
                        cost_increment.checked_add(&config.length_costs.evaluate(&length))
                    else {
                        return Cost::max_value();
                    };
                    let Some(cost_increment) = cost_increment
                        .checked_add(&config.length_difference_costs.evaluate(&length_difference))
                    else {
                        return Cost::max_value();
                    };
                    cost_increment
                }
                AlignmentType::AlternativeStart {
                    reference_index: reference_start,
                    query_index: query_start,
                } => {
                    assert_eq!(reference_index, reference_start);
                    assert_eq!(query_index, query_start);
                    Cost::zero()
                }
                AlignmentType::PrimaryReentry {
                    reference_index: reference_reentry_index,
                    query_index: query_reentry_index,
                } => {
                    assert_eq!(reference_index, reference_reentry_index);
                    assert_eq!(query_index, query_reentry_index);
                    Cost::zero()
                }
                AlignmentType::Root | AlignmentType::SecondaryRoot => {
                    // Do nothing
                    Cost::zero()
                }
                AlignmentType::PrimaryShortcut { .. } => panic!("Not supported"),
            };

            cost = if let Some(cost) = cost.checked_add(&cost_increment) {
                cost
            } else {
                return Cost::max_value();
            };
            last_alignment_type = Some(alignment_type);
        }

        cost
    }
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use compact_genome::{
        implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
        interface::{
            alphabet::Alphabet,
            sequence::{GenomeSequence, OwnedGenomeSequence},
        },
    };
    use generic_a_star::cost::U64Cost;

    use crate::{
        a_star_aligner::{
            alignment_result::alignment::Alignment,
            template_switch_distance::{
                AlignmentType, EqualCostRange, TemplateSwitchAncestor, TemplateSwitchDescendant,
                TemplateSwitchDirection,
            },
        },
        config::{BaseCost, TemplateSwitchConfig},
        costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    };

    static START_REFERENCE: &[u8] = b"AGAGAGCTCTAA";
    static START_QUERY: &[u8] = b"AGAGAGCTTTAA";
    static START_COSTS: LazyLock<Vec<U64Cost>> = LazyLock::new(|| {
        [
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&-6)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&2)
                + CONFIG.length_costs.evaluate(&2)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&-4)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&3)
                + CONFIG.length_costs.evaluate(&3)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&-2)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&4)
                + CONFIG.length_costs.evaluate(&4)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&0)
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'G').unwrap(),
                    DnaAlphabet::ascii_to_character(b'T').unwrap(),
                )
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&5)
                + CONFIG.length_costs.evaluate(&5)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&2)
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'G').unwrap(),
                    DnaAlphabet::ascii_to_character(b'T').unwrap(),
                )
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'A').unwrap(),
                    DnaAlphabet::ascii_to_character(b'C').unwrap(),
                )
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&6)
                + CONFIG.length_costs.evaluate(&6)
                + CONFIG.length_difference_costs.evaluate(&0),
        ]
        .to_vec()
    });
    static START_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
        &[
            (6, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -6,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 2,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (5, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -4,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (3, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 3,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (4, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -2,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 4,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (3, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 0,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (1, AlignmentType::SecondarySubstitution),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 5,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (2, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 2,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondarySubstitution),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 6,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
    ];

    static END_REFERENCE: &[u8] = b"AACTCTAGAGAG";
    static END_QUERY: &[u8] = b"AATTCTAGAGAG";
    static END_COSTS: LazyLock<Vec<U64Cost>> = LazyLock::new(|| {
        [
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&10)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&2)
                + CONFIG.length_costs.evaluate(&2)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&10)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&3)
                + CONFIG.length_costs.evaluate(&3)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&10)
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&4)
                + CONFIG.length_costs.evaluate(&4)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&10)
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'A').unwrap(),
                    DnaAlphabet::ascii_to_character(b'C').unwrap(),
                )
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&5)
                + CONFIG.length_costs.evaluate(&5)
                + CONFIG.length_difference_costs.evaluate(&0),
            CONFIG.base_cost.rqr
                + CONFIG.rq_qr_offset_costs.evaluate(&10)
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'A').unwrap(),
                    DnaAlphabet::ascii_to_character(b'C').unwrap(),
                )
                + CONFIG.secondary_reverse_edit_costs.substitution_cost(
                    DnaAlphabet::ascii_to_character(b'G').unwrap(),
                    DnaAlphabet::ascii_to_character(b'T').unwrap(),
                )
                + CONFIG.reverse_anti_descendant_gap_costs.evaluate(&6)
                + CONFIG.length_costs.evaluate(&6)
                + CONFIG.length_difference_costs.evaluate(&0),
        ]
        .to_vec()
    });
    static END_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 2,
                },
            ),
            (6, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (3, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 3,
                },
            ),
            (5, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 4,
                },
            ),
            (4, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 5,
                },
            ),
            (3, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    descendant: TemplateSwitchDescendant::Reference,
                    ancestor: TemplateSwitchAncestor::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (2, AlignmentType::SecondarySubstitution),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_descendant_gap: 6,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
    ];

    static CONFIG: LazyLock<TemplateSwitchConfig<DnaAlphabet, U64Cost>> =
        LazyLock::new(|| TemplateSwitchConfig {
            left_flank_length: 0,
            right_flank_length: 0,
            template_switch_min_length: 3,
            base_cost: BaseCost {
                rrf: 10u64.into(),
                rqf: 100u64.into(),
                qrf: 1000u64.into(),
                qqf: 10000u64.into(),
                rrr: 100000u64.into(),
                rqr: 1000000u64.into(),
                qrr: 10000000u64.into(),
                qqr: 100000000u64.into(),
            },
            primary_edit_costs: GapAffineAlignmentCostTable::new(
                "",
                [0u64, 2, 2, 2, 2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 0]
                    .map(Into::into)
                    .to_vec(),
                [3u64, 3, 3, 3].map(Into::into).to_vec(),
                [1u64, 1, 1, 1].map(Into::into).to_vec(),
            ),
            secondary_forward_edit_costs: GapAffineAlignmentCostTable::new(
                "",
                [0u64, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0]
                    .map(Into::into)
                    .to_vec(),
                [3u64, 6, 6, 6].map(Into::into).to_vec(),
                [1u64, 1, 1, 1].map(Into::into).to_vec(),
            ),
            secondary_reverse_edit_costs: GapAffineAlignmentCostTable::new(
                "",
                [0u64, 5, 5, 5, 5, 0, 5, 5, 5, 5, 0, 5, 5, 5, 5, 0]
                    .map(Into::into)
                    .to_vec(),
                [3u64, 7, 7, 7].map(Into::into).to_vec(),
                [1u64, 1, 1, 1].map(Into::into).to_vec(),
            ),
            left_flank_edit_costs: GapAffineAlignmentCostTable::new_zero(),
            right_flank_edit_costs: GapAffineAlignmentCostTable::new_zero(),
            rq_qr_offset_costs: CostFunction::try_from(
                (-20..=20)
                    .map(|i| (i, U64Cost::from(17 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            rr_qq_offset_costs: CostFunction::try_from(
                (-20..=20)
                    .map(|i| (i, U64Cost::from(17 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            length_costs: CostFunction::try_from(
                (0..=20)
                    .map(|i| (i, U64Cost::from(19 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            length_difference_costs: CostFunction::try_from(
                (-20..=20)
                    .map(|i| (i, U64Cost::from(23 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            forward_anti_descendant_gap_costs: CostFunction::try_from(
                (-20..=20)
                    .map(|i| (i, U64Cost::from(29 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            reverse_anti_descendant_gap_costs: CostFunction::try_from(
                (-20..=20)
                    .map(|i| (i, U64Cost::from(31 * u64::try_from(i + 21).unwrap())))
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
        });

    #[test]
    fn move_template_switch_start_backwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
        let mut alignment = Alignment::from(START_ALIGNMENTS[0].to_vec());
        assert_eq!(
            alignment.compute_cost(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &CONFIG
            ),
            START_COSTS[0]
        );

        for (expected_alignment, expected_cost) in
            START_ALIGNMENTS[1..].iter().zip(&START_COSTS[1..])
        {
            assert!(alignment.move_template_switch_start_backwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &mut 1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
            assert_eq!(
                alignment.compute_cost(
                    reference.as_genome_subsequence(),
                    query.as_genome_subsequence(),
                    2,
                    2,
                    &CONFIG
                ),
                *expected_cost
            );
        }
    }

    #[test]
    fn move_template_switch_start_forwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
        let mut alignment = Alignment::from(START_ALIGNMENTS.last().unwrap().to_vec());
        assert_eq!(
            alignment.compute_cost(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &CONFIG
            ),
            *START_COSTS.last().unwrap()
        );

        for (expected_alignment, expected_cost) in START_ALIGNMENTS
            .iter()
            .zip(START_COSTS.iter())
            .rev()
            .skip(1)
        {
            assert!(alignment.move_template_switch_start_forwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                &mut 1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
            assert_eq!(
                alignment.compute_cost(
                    reference.as_genome_subsequence(),
                    query.as_genome_subsequence(),
                    2,
                    2,
                    &CONFIG
                ),
                *expected_cost
            );
        }
    }

    #[test]
    fn move_template_switch_end_backwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
        let mut alignment = Alignment::from(END_ALIGNMENTS.last().unwrap().to_vec());
        assert_eq!(
            alignment.compute_cost(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                &CONFIG
            ),
            *END_COSTS.last().unwrap()
        );

        for (expected_alignment, expected_cost) in
            END_ALIGNMENTS.iter().zip(END_COSTS.iter()).rev().skip(1)
        {
            assert!(alignment.move_template_switch_end_backwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
            assert_eq!(
                alignment.compute_cost(
                    reference.as_genome_subsequence(),
                    query.as_genome_subsequence(),
                    1,
                    1,
                    &CONFIG
                ),
                *expected_cost
            );
        }
    }

    #[test]
    fn move_template_switch_end_forwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
        let mut alignment = Alignment::from(END_ALIGNMENTS[0].to_vec());
        assert_eq!(
            alignment.compute_cost(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                &CONFIG
            ),
            *END_COSTS.first().unwrap()
        );

        for (expected_alignment, expected_cost) in END_ALIGNMENTS[1..].iter().zip(&END_COSTS[1..]) {
            assert!(alignment.move_template_switch_end_forwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
            assert_eq!(
                alignment.compute_cost(
                    reference.as_genome_subsequence(),
                    query.as_genome_subsequence(),
                    1,
                    1,
                    &CONFIG
                ),
                *expected_cost
            );
        }
    }
}
