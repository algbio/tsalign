use std::{collections::VecDeque, iter};

use crate::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{AlignmentType, TemplateSwitchDescendant},
};

#[derive(Debug, Clone)]
pub struct AlignmentStream {
    stream: VecDeque<(usize, AlignmentType)>,
    length: usize,

    /// The alignment coordinates of the first alignment after the head of the queue.
    head_coordinates: AlignmentStreamCoordinates,
    /// The alignment coordinates of the tailmost element in the queue.
    tail_coordinates: AlignmentStreamCoordinates,
}

#[derive(Debug, Clone, Copy)]
pub struct AlignmentStreamCoordinates {
    reference: usize,
    query: usize,
    template_switch_descendant: Option<TemplateSwitchDescendant>,
}

impl AlignmentStream {
    #[deprecated = "our alignments typically have offsets nowadays"]
    #[expect(clippy::new_without_default)]
    pub fn new() -> Self {
        Self::new_with_offset(0, 0)
    }

    pub fn new_with_offset(reference_offset: usize, query_offset: usize) -> Self {
        Self {
            stream: Default::default(),
            length: 0,
            head_coordinates: AlignmentStreamCoordinates::new_with_offset(
                reference_offset,
                query_offset,
            ),
            tail_coordinates: AlignmentStreamCoordinates::new_with_offset(
                reference_offset,
                query_offset,
            ),
        }
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn stream_iter(&self) -> impl use<'_> + DoubleEndedIterator<Item = (usize, AlignmentType)> {
        self.stream.iter().copied()
    }

    pub fn stream_iter_flat(&self) -> impl use<'_> + Iterator<Item = AlignmentType> {
        self.stream_iter()
            .flat_map(|(multiplicity, alignment_type)| iter::repeat_n(alignment_type, multiplicity))
    }

    pub fn stream_iter_flat_rev(&self) -> impl use<'_> + Iterator<Item = AlignmentType> {
        self.stream_iter()
            .rev()
            .flat_map(|(multiplicity, alignment_type)| iter::repeat_n(alignment_type, multiplicity))
    }

    fn stream_vec(&self) -> Vec<(usize, AlignmentType)> {
        self.stream_iter().collect()
    }

    pub fn stream_alignment(&self) -> Alignment<AlignmentType> {
        self.stream_vec().into()
    }

    pub fn head_coordinates(&self) -> AlignmentStreamCoordinates {
        self.head_coordinates
    }

    pub fn tail_coordinates(&self) -> AlignmentStreamCoordinates {
        self.tail_coordinates
    }

    pub fn push_until_full(
        &mut self,
        multiplicity: &mut usize,
        alignment_type: AlignmentType,
        requested_length: usize,
    ) {
        let available_length = requested_length - self.length;
        let push_length = *multiplicity * Self::stream_length(alignment_type);

        if available_length >= push_length {
            self.push(*multiplicity, alignment_type);
            *multiplicity = 0;
        } else {
            let push_multiplicity = available_length.div_ceil(Self::stream_length(alignment_type));
            self.push(push_multiplicity, alignment_type);
            *multiplicity -= push_multiplicity;
        }
    }

    pub fn clear(&mut self) {
        self.pop(0);
    }

    pub fn is_full(&self, requested_length: usize) -> bool {
        self.length >= requested_length
    }

    pub fn is_empty(&self) -> bool {
        self.stream_iter().next().is_none()
    }

    pub fn push(&mut self, multiplicity: usize, alignment_type: AlignmentType) {
        self.stream.push_back((multiplicity, alignment_type));
        self.head_coordinates.advance(multiplicity, alignment_type);
        self.length += multiplicity * Self::stream_length(alignment_type);
    }

    pub fn push_all(&mut self, iterator: impl IntoIterator<Item = (usize, AlignmentType)>) {
        for (multiplicity, alignment_type) in iterator {
            self.push(multiplicity, alignment_type);
        }
    }

    /// Pops one unit of length from the tail of the stream.
    pub fn pop_one(&mut self) {
        self.pop(self.length.saturating_sub(1));
    }

    pub fn pop(&mut self, requested_length: usize) {
        while self.length > requested_length {
            let requested_pop_length = self.length - requested_length;
            let (multiplicity, alignment_type) = self.stream.front_mut().unwrap();
            let front_length = *multiplicity * Self::stream_length(*alignment_type);

            if front_length <= requested_pop_length {
                self.tail_coordinates
                    .advance(*multiplicity, *alignment_type);
                self.stream.pop_front();
                self.length -= front_length;
            } else {
                let pop_multiplicity = requested_pop_length / Self::stream_length(*alignment_type);
                self.tail_coordinates
                    .advance(pop_multiplicity, *alignment_type);
                *multiplicity -= pop_multiplicity;
                self.length -= pop_multiplicity * Self::stream_length(*alignment_type);
                break;
            }
        }

        while let Some((multiplicity, alignment_type)) = self.stream.front().copied() {
            if Self::stream_length(alignment_type) == 0 {
                self.tail_coordinates.advance(multiplicity, alignment_type);
                self.stream.pop_front();
            } else {
                break;
            }
        }
    }

    fn stream_length(alignment_type: AlignmentType) -> usize {
        match alignment_type {
            AlignmentType::PrimaryInsertion
            | AlignmentType::PrimaryDeletion
            | AlignmentType::PrimarySubstitution
            | AlignmentType::PrimaryMatch
            | AlignmentType::PrimaryFlankInsertion
            | AlignmentType::PrimaryFlankDeletion
            | AlignmentType::PrimaryFlankSubstitution
            | AlignmentType::PrimaryFlankMatch
            | AlignmentType::SecondaryInsertion
            | AlignmentType::SecondaryDeletion
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryMatch => 1,
            AlignmentType::TemplateSwitchEntrance { .. }
            | AlignmentType::TemplateSwitchExit { .. }
            | AlignmentType::Root
            | AlignmentType::SecondaryRoot
            | AlignmentType::PrimaryReentry => 0,
            AlignmentType::PrimaryShortcut { .. } => {
                unreachable!("Shortcut alignments are not supported for show")
            }
        }
    }
}

impl AlignmentStreamCoordinates {
    pub fn new_with_offset(reference_offset: usize, query_offset: usize) -> Self {
        Self {
            reference: reference_offset,
            query: query_offset,
            template_switch_descendant: None,
        }
    }

    pub fn reference(&self) -> usize {
        self.reference
    }

    pub fn query(&self) -> usize {
        self.query
    }

    fn advance(&mut self, multiplicity: usize, alignment_type: AlignmentType) {
        let (reference_length, query_length) = match alignment_type {
            AlignmentType::PrimaryInsertion | AlignmentType::PrimaryFlankInsertion => (0, 1),
            AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => (1, 0),
            AlignmentType::PrimarySubstitution
            | AlignmentType::PrimaryMatch
            | AlignmentType::PrimaryFlankSubstitution
            | AlignmentType::PrimaryFlankMatch => (1, 1),
            AlignmentType::TemplateSwitchEntrance { descendant, .. } => {
                assert!(
                    self.template_switch_descendant.is_none(),
                    "Encountered template switch entrance within template switch"
                );
                self.template_switch_descendant = Some(descendant);
                (0, 0)
            }
            AlignmentType::SecondaryInsertion
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryMatch => match self.template_switch_descendant.unwrap() {
                TemplateSwitchDescendant::Reference => (1, 0),
                TemplateSwitchDescendant::Query => (0, 1),
            },
            AlignmentType::TemplateSwitchExit {
                anti_descendant_gap,
            } => {
                let Some(template_switch_descendant) = self.template_switch_descendant.take()
                else {
                    panic!(
                        "Encountered template switch exit without first encountering a template switch entrance"
                    )
                };
                match template_switch_descendant {
                    TemplateSwitchDescendant::Reference => {
                        self.query =
                            usize::try_from(self.query as isize + anti_descendant_gap).unwrap()
                    }
                    TemplateSwitchDescendant::Query => {
                        self.reference =
                            usize::try_from(self.reference as isize + anti_descendant_gap).unwrap()
                    }
                }
                (0, 0)
            }
            AlignmentType::SecondaryDeletion
            | AlignmentType::Root
            | AlignmentType::SecondaryRoot
            | AlignmentType::PrimaryReentry => (0, 0),
            AlignmentType::PrimaryShortcut { .. } => {
                unreachable!("Shortcut alignments are not supported for show")
            }
        };

        self.reference += multiplicity * reference_length;
        self.query += multiplicity * query_length;
    }
}
