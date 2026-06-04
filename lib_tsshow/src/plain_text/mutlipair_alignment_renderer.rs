use std::{
    collections::BTreeMap,
    fmt::{Debug, Display},
    iter, mem,
    ops::{Index, IndexMut},
};

use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;
use log::{debug, trace};

#[cfg(test)]
mod tests;

#[derive(Debug)]
pub struct MultipairAlignmentRenderer<SequenceName, CharacterData = NoCharacterData> {
    sequences: BTreeMap<SequenceName, MultipairAlignmentSequence<CharacterData>>,
}

#[derive(Debug)]
pub struct MultipairAlignmentSequence<CharacterData> {
    sequence: Vec<Character<CharacterData>>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Character<Data> {
    kind: CharacterKind,
    data: Data,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum CharacterKind {
    Char(char),
    Gap,
    Blank,
}

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq)]
pub struct NoCharacterData;

pub enum AlignmentRenderResult {
    Ok {
        /// Column of the first inserted query character or gap.
        query_offset_column: usize,
        /// Column after the last inserted query character or gap.
        query_limit_column: usize,
    },

    /// There was nothing to align.
    Empty,
}

impl<SequenceName: Eq + Ord, CharacterData>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    pub fn new(
        root_sequence_name: SequenceName,
        root_sequence: impl IntoIterator<Item = Character<CharacterData>>,
    ) -> Self {
        debug!("Adding root sequence");

        Self {
            sequences: [(root_sequence_name, root_sequence.into_iter().collect())]
                .into_iter()
                .collect(),
        }
    }

    pub fn new_empty() -> Self {
        debug!("Creating an empty renderer without root sequence");

        Self {
            sequences: Default::default(),
        }
    }

    pub fn sequence(
        &self,
        sequence_name: &SequenceName,
    ) -> &MultipairAlignmentSequence<CharacterData> {
        self.sequences.get(sequence_name).unwrap()
    }

    pub fn column_width(&self) -> usize {
        debug_assert!({
            let mut column_width = None;
            for sequence in self.sequences.values() {
                if let Some(column_width) = &column_width {
                    debug_assert_eq!(*column_width, sequence.len());
                } else {
                    column_width = Some(sequence.len());
                }
            }
            true
        });

        self.sequences.values().next().unwrap().len()
    }

    /// Append a sequence to the end of an existing rendered sequence.
    ///
    /// Existing gaps and blanks at the end of the existing rendered sequence will not be removed.
    /// The extension will be added after those.
    pub fn extend_sequence(
        &mut self,
        sequence_name: &SequenceName,
        extension: impl IntoIterator<Item = Character<CharacterData>>,
        mut blank_data_generator: impl FnMut() -> CharacterData,
    ) {
        debug!("Extending sequence");

        let sequence = self.sequences.get_mut(sequence_name).unwrap();
        trace!("Length before: {}", sequence.len());

        // Extend
        sequence.extend_with(extension);
        let new_length = sequence.len();
        trace!("Length after: {new_length}");

        // Add blanks to other sequences to make all sequences of the same length again
        for (other_sequence_name, other_sequence) in &mut self.sequences {
            if other_sequence_name != sequence_name {
                other_sequence.extend_with_blanks(&mut blank_data_generator, new_length);
            }
        }
    }

    /// Append a sequence to an existing sequence while aligning it to another sequence.
    #[expect(clippy::too_many_arguments)]
    pub fn extend_sequence_with_alignment(
        &mut self,
        reference_sequence_name: &SequenceName,
        query_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        extension: impl IntoIterator<Item = Character<CharacterData>>,
        blank_data_generator: impl FnMut() -> CharacterData,
        gap_data_generator: impl FnMut() -> CharacterData,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        debug!("Extending sequence with alignment");
        debug!("reference_sequence_offset: {reference_sequence_offset}");

        let rendered_reference_sequence = self.sequence(reference_sequence_name);
        let rendered_query_sequence = self.sequence(query_sequence_name);

        let mut rendered_sequence_offset = rendered_reference_sequence
            .translate_extension_offset(reference_sequence_offset)
            .unwrap_or_else(|| {
                panic!(
                    "sequence_offset {reference_sequence_offset} is out of bounds (character count: {})",
                    rendered_reference_sequence.character_count()
                )
            });

        while rendered_sequence_offset > 0 {
            if matches!(
                rendered_reference_sequence
                    .get(rendered_sequence_offset - 1)
                    .map(|c| c.kind()),
                Some(CharacterKind::Char(_))
            ) || !matches!(
                rendered_query_sequence
                    .get(rendered_sequence_offset - 1)
                    .map(|c| c.kind()),
                Some(CharacterKind::Blank)
            ) {
                break;
            }

            rendered_sequence_offset -= 1;
        }
        trace!("rendered_sequence_offset: {rendered_sequence_offset}");
        if let Some(reference_offset_character) =
            rendered_reference_sequence.get(rendered_sequence_offset)
        {
            trace!(
                "reference_offset_character: {:?}",
                reference_offset_character.as_char()
            );
        } else {
            trace!("reference_offset_character: offset is past the end of the reference");
        }

        self.sequences
            .get_mut(query_sequence_name)
            .unwrap()
            .prune_blanks(rendered_sequence_offset);
        self.extend_sequence_with_alignment_internal(
            reference_sequence_name,
            query_sequence_name,
            rendered_sequence_offset,
            extension,
            blank_data_generator,
            gap_data_generator,
            alignment,
            do_lowercasing,
            invert_alignment,
        );
    }

    #[expect(clippy::too_many_arguments)]
    fn extend_sequence_with_alignment_internal(
        &mut self,
        reference_sequence_name: &SequenceName,
        query_sequence_name: &SequenceName,
        rendered_sequence_offset: usize,
        extension: impl IntoIterator<Item = Character<CharacterData>>,
        mut blank_data_generator: impl FnMut() -> CharacterData,
        mut gap_data_generator: impl FnMut() -> CharacterData,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) -> AlignmentRenderResult {
        let alignment = alignment.into_iter();
        let mut reference_gaps = Vec::new();

        let (reference_sequence_name, mut translated_reference_sequence) = self
            .sequences
            .remove_entry(reference_sequence_name)
            .unwrap();
        trace!(
            "translated_reference_sequence.len(): {}",
            translated_reference_sequence.len()
        );
        trace!("translated_reference_sequence_offset: {rendered_sequence_offset}");
        trace!(
            "alignment_length: ({}, {})",
            alignment.size_hint().0,
            alignment
                .size_hint()
                .1
                .map_or(String::from("unknown"), |len| len.to_string())
        );

        let (query_sequence_name, mut translated_query_sequence) =
            self.sequences.remove_entry(query_sequence_name).unwrap();
        assert_eq!(translated_query_sequence.len(), rendered_sequence_offset);

        let mut index = rendered_sequence_offset;
        let mut extension = extension.into_iter();
        let mut query_offset_column = None;
        let mut query_limit_column = None;

        for alignment_type in alignment.map(|alignment_type| {
            if invert_alignment {
                alignment_type.inverted()
            } else {
                alignment_type
            }
        }) {
            trace!("alignment_type: {alignment_type}");

            while matches!(
                translated_reference_sequence
                    .get(index)
                    .map(Character::kind),
                Some(CharacterKind::Blank)
            ) {
                trace!("Skipping blank");
                translated_query_sequence.push(Character::new_blank(blank_data_generator()));
                index += 1;
            }

            match alignment_type {
                AlignmentType::PrimaryInsertion
                | AlignmentType::PrimaryFlankInsertion
                | AlignmentType::SecondaryInsertion => {
                    query_offset_column = query_offset_column.or(Some(index));

                    if matches!(
                        translated_reference_sequence
                            .get(index)
                            .map(Character::kind),
                        Some(CharacterKind::Gap)
                    ) {
                        index += 1;
                    } else {
                        reference_gaps.push(index);
                    }

                    translated_query_sequence.push(extension.next().unwrap());
                    query_limit_column = Some(index);
                }
                AlignmentType::PrimaryDeletion
                | AlignmentType::PrimaryFlankDeletion
                | AlignmentType::SecondaryDeletion => {
                    while matches!(
                        translated_reference_sequence
                            .get(index)
                            .map(Character::kind),
                        Some(CharacterKind::Gap | CharacterKind::Blank)
                    ) {
                        translated_query_sequence
                            .push(Character::new_blank(blank_data_generator()));
                        index += 1;
                    }

                    query_offset_column = query_offset_column.or(Some(index));
                    translated_query_sequence.push(Character::new_gap(gap_data_generator()));
                    index += 1;
                    query_limit_column = Some(index);
                }
                AlignmentType::PrimarySubstitution
                | AlignmentType::PrimaryFlankSubstitution
                | AlignmentType::SecondarySubstitution => {
                    while matches!(
                        translated_reference_sequence
                            .get(index)
                            .map(Character::kind),
                        Some(CharacterKind::Gap | CharacterKind::Blank)
                    ) {
                        translated_query_sequence
                            .push(Character::new_blank(blank_data_generator()));
                        index += 1;
                    }
                    let mut extension_character = extension.next().unwrap();
                    assert!(extension_character.is_char());

                    if do_lowercasing {
                        extension_character.make_ascii_lowercase();
                        translated_reference_sequence[index].make_ascii_lowercase();
                    }

                    query_offset_column = query_offset_column.or(Some(index));
                    translated_query_sequence.push(extension_character);
                    index += 1;
                    query_limit_column = Some(index);
                }
                AlignmentType::PrimaryMatch
                | AlignmentType::PrimaryFlankMatch
                | AlignmentType::SecondaryMatch => {
                    while matches!(
                        translated_reference_sequence
                            .get(index)
                            .map(Character::kind),
                        Some(CharacterKind::Gap | CharacterKind::Blank)
                    ) {
                        trace!("Skipping gap or blank before match");
                        translated_query_sequence
                            .push(Character::new_blank(blank_data_generator()));
                        index += 1;
                    }

                    query_offset_column = query_offset_column.or(Some(index));
                    translated_query_sequence.push(extension.next().unwrap());
                    index += 1;
                    query_limit_column = Some(index);
                }
                AlignmentType::Root
                | AlignmentType::PrimaryReentry { .. }
                | AlignmentType::TemplateSwitchEntrance { .. }
                | AlignmentType::TemplateSwitchExit { .. }
                | AlignmentType::AlternativeStart { .. }
                | AlignmentType::SecondaryRoot
                | AlignmentType::PrimaryShortcut { .. } => {
                    panic!("Not allowed in rendered alignment: {alignment_type:?}")
                }
            }

            assert!(
                index <= translated_reference_sequence.len(),
                "index: {index}, translated_reference_sequence.len(): {}",
                translated_reference_sequence.len(),
            );
        }

        assert!(extension.next().is_none());
        assert_eq!(query_offset_column.is_some(), query_limit_column.is_some());

        translated_query_sequence.extend_with_blanks(
            &mut blank_data_generator,
            translated_reference_sequence.len(),
        );
        translated_reference_sequence
            .insert_gaps(gap_data_generator, reference_gaps.iter().copied());

        for sequence in self.sequences.values_mut() {
            sequence.insert_blanks(&mut blank_data_generator, reference_gaps.iter().copied());
        }

        self.sequences
            .insert(reference_sequence_name, translated_reference_sequence);
        self.sequences
            .insert(query_sequence_name, translated_query_sequence);

        if let (Some(query_offset_column), Some(query_limit_column)) =
            (query_offset_column, query_limit_column)
        {
            AlignmentRenderResult::Ok {
                query_offset_column,
                query_limit_column,
            }
        } else {
            AlignmentRenderResult::Empty
        }
    }
}

impl<SequenceName: Eq + Ord + Clone, CharacterData>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    #[expect(clippy::too_many_arguments)]
    pub fn add_aligned_sequence(
        &mut self,
        reference_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        query_sequence_name: SequenceName,
        query_sequence: impl IntoIterator<Item = Character<CharacterData>>,
        mut blank_data_generator: impl FnMut() -> CharacterData,
        gap_data_generator: impl FnMut() -> CharacterData,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) -> AlignmentRenderResult {
        debug!("Adding aligned sequence");
        debug!("reference_offset: {reference_sequence_offset}");
        debug!("invert_alignment: {invert_alignment}");

        assert!(!self.sequences.contains_key(&query_sequence_name));

        let reference_sequence = self
            .sequences
            .get_mut(reference_sequence_name)
            .unwrap_or_else(|| panic!("Reference sequence does not exist"));
        let index = reference_sequence
            .translate_alignment_offset(reference_sequence_offset)
            .unwrap_or_else(|| {
                panic!("reference_sequence_offset {reference_sequence_offset} is out of bounds")
            });
        self.sequences.insert(
            query_sequence_name.clone(),
            iter::repeat_with(|| Character::new_blank(blank_data_generator()))
                .take(index)
                .collect(),
        );

        self.extend_sequence_with_alignment_internal(
            reference_sequence_name,
            &query_sequence_name,
            index,
            query_sequence,
            blank_data_generator,
            gap_data_generator,
            alignment,
            do_lowercasing,
            invert_alignment,
        )
    }

    pub fn add_independent_sequence(
        &mut self,
        sequence_name: SequenceName,
        sequence: impl IntoIterator<Item = Character<CharacterData>>,
    ) {
        assert!(!self.sequences.contains_key(&sequence_name));
        self.sequences
            .insert(sequence_name, sequence.into_iter().collect());
    }

    pub fn add_empty_independent_sequence(&mut self, sequence_name: SequenceName) {
        self.add_independent_sequence(sequence_name, None);
    }
}

impl<SequenceName: Eq + Ord, CharacterData: Default>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    pub fn extend_sequence_with_default_data(
        &mut self,
        sequence_name: &SequenceName,
        extension: impl IntoIterator<Item = char>,
    ) {
        self.extend_sequence(
            sequence_name,
            extension.into_iter().map(Character::new_char_with_default),
            Default::default,
        );
    }

    #[expect(clippy::too_many_arguments)]
    pub fn extend_sequence_with_alignment_and_default_data(
        &mut self,
        reference_sequence_name: &SequenceName,
        query_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        extension: impl IntoIterator<Item = char>,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        self.extend_sequence_with_alignment(
            reference_sequence_name,
            query_sequence_name,
            reference_sequence_offset,
            extension.into_iter().map(Character::new_char_with_default),
            Default::default,
            Default::default,
            alignment,
            do_lowercasing,
            invert_alignment,
        );
    }
}

impl<SequenceName: Eq + Ord + Clone, CharacterData: Default>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    #[expect(clippy::too_many_arguments)]
    pub fn add_aligned_sequence_with_default_data(
        &mut self,
        reference_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        query_sequence_name: SequenceName,
        query_sequence: impl IntoIterator<Item = char>,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) -> AlignmentRenderResult {
        self.add_aligned_sequence(
            reference_sequence_name,
            reference_sequence_offset,
            query_sequence_name,
            query_sequence
                .into_iter()
                .map(Character::new_char_with_default),
            Default::default,
            Default::default,
            alignment,
            do_lowercasing,
            invert_alignment,
        )
    }
}

impl<SequenceName: Eq + Ord + Clone> MultipairAlignmentRenderer<SequenceName> {
    pub fn new_without_data(
        root_sequence_name: SequenceName,
        root_sequence: impl IntoIterator<Item = char>,
    ) -> Self {
        Self::new(
            root_sequence_name,
            root_sequence
                .into_iter()
                .map(Character::new_char_with_default),
        )
    }

    #[expect(clippy::too_many_arguments)]
    pub fn add_aligned_sequence_without_data(
        &mut self,
        reference_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        query_sequence_name: SequenceName,
        query_sequence: impl IntoIterator<Item = char>,
        alignment: impl IntoIterator<Item = AlignmentType>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        self.add_aligned_sequence(
            reference_sequence_name,
            reference_sequence_offset,
            query_sequence_name,
            query_sequence
                .into_iter()
                .map(Character::new_char_with_default),
            || NoCharacterData,
            || NoCharacterData,
            alignment,
            do_lowercasing,
            invert_alignment,
        );
    }
}

impl<SequenceName: Eq + Ord + Display, CharacterData>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    pub fn render<'name>(
        &self,
        mut output: impl std::io::Write,
        names: impl IntoIterator<Item = &'name SequenceName>,
    ) -> Result<(), std::io::Error>
    where
        SequenceName: 'name,
    {
        let names: Vec<_> = names.into_iter().collect();
        let max_name_len = names
            .iter()
            .map(ToString::to_string)
            .map(|name| name.chars().count())
            .max()
            .unwrap();

        for name in names {
            let sequence = self
                .sequences
                .get(name)
                .unwrap_or_else(|| panic!("Render sequence name not found: {name}"));

            let name = name.to_string();
            write!(output, "{name}: ")?;
            for _ in name.len()..max_name_len {
                write!(output, " ")?;
            }

            writeln!(output, "{sequence}")?;
        }

        Ok(())
    }
}

impl<SequenceName: Eq + Ord, CharacterData>
    MultipairAlignmentRenderer<SequenceName, CharacterData>
{
    #[allow(unused)]
    pub fn render_without_names<'name>(
        &self,
        mut output: impl std::io::Write,
        names: impl IntoIterator<Item = &'name SequenceName>,
    ) -> Result<(), std::io::Error>
    where
        SequenceName: 'name,
    {
        let names: Vec<_> = names.into_iter().collect();

        for name in names {
            let sequence = self.sequences.get(name).unwrap();
            writeln!(output, "{sequence}")?;
        }

        Ok(())
    }
}

impl<CharacterData> MultipairAlignmentSequence<CharacterData> {
    /// Returns the smallest index that skips the first `offset` characters.
    ///
    /// Returns **`None`** if there are less than `offset` characters.
    ///
    /// # Example
    ///
    /// ```rust
    /// use lib_tsshow::plain_text::mutlipair_alignment_renderer::{Character, CharacterKind, MultipairAlignmentSequence, NoCharacterData};
    ///
    /// let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([CharacterKind::Blank, CharacterKind::Char('A'), CharacterKind::Gap, CharacterKind::Char('C'), CharacterKind::Gap]);
    /// assert_eq!(sequence.translate_alignment_offset(0), Some(0));
    /// assert_eq!(sequence.translate_alignment_offset(1), Some(2));
    /// assert_eq!(sequence.translate_alignment_offset(2), Some(4));
    /// assert_eq!(sequence.translate_alignment_offset(3), None);
    /// ```
    pub fn translate_alignment_offset(&self, offset: usize) -> Option<usize> {
        if offset == 0 {
            Some(0)
        } else {
            self.sequence
                .iter()
                .enumerate()
                .filter(|(_, character)| matches!(character.kind(), CharacterKind::Char { .. }))
                .nth(offset - 1)
                .map(|(index, _)| index + 1)
        }
    }

    /// Returns the largest index that skips only the first `offset` characters (but not the first `offset + 1` characters).
    ///
    /// Returns **`None`** if there are less than `offset` characters.
    ///
    /// # Example
    ///
    /// ```rust
    /// use lib_tsshow::plain_text::mutlipair_alignment_renderer::{Character, CharacterKind, MultipairAlignmentSequence, NoCharacterData};
    ///
    /// let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([CharacterKind::Blank, CharacterKind::Char('A'), CharacterKind::Gap, CharacterKind::Char('C'), CharacterKind::Gap]);
    /// assert_eq!(sequence.translate_extension_offset(0), Some(1));
    /// assert_eq!(sequence.translate_extension_offset(1), Some(3));
    /// assert_eq!(sequence.translate_extension_offset(2), Some(5));
    /// assert_eq!(sequence.translate_extension_offset(3), None);
    /// ```
    pub fn translate_extension_offset(&self, offset: usize) -> Option<usize> {
        self.translate_alignment_offset(offset).map(|offset| {
            self.sequence
                .iter()
                .enumerate()
                .skip(offset)
                .take_while(|(_, character)| !matches!(character.kind(), CharacterKind::Char(_)))
                .last()
                .map(|(index, _)| index + 1)
                .unwrap_or(offset)
        })
    }

    /// Returns the index of the first gap or character after skipping `offset` gaps and characters.
    ///
    /// Returns **`None`** if there are less than or equal to `offset` gaps and characters.
    ///
    /// # Example
    ///
    /// ```rust
    /// use lib_tsshow::plain_text::mutlipair_alignment_renderer::{Character, CharacterKind, MultipairAlignmentSequence, NoCharacterData};
    ///
    /// let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([CharacterKind::Blank, CharacterKind::Char('A'), CharacterKind::Gap, CharacterKind::Char('C'), CharacterKind::Blank]);
    /// assert_eq!(sequence.translate_offset_without_blanks(0), Some(1));
    /// assert_eq!(sequence.translate_offset_without_blanks(1), Some(2));
    /// assert_eq!(sequence.translate_offset_without_blanks(2), Some(3));
    /// assert_eq!(sequence.translate_offset_without_blanks(3), None);
    /// ```
    pub fn translate_offset_without_blanks(&self, offset: usize) -> Option<usize> {
        self.iter()
            .enumerate()
            .filter_map(|(index, character)| {
                if matches!(character.kind, CharacterKind::Blank) {
                    None
                } else {
                    Some(index)
                }
            })
            .nth(offset)
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn len_without_blanks(&self) -> usize {
        self.sequence
            .iter()
            .filter(|character| !matches!(character.kind, CharacterKind::Blank))
            .count()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn character_count(&self) -> usize {
        self.iter().filter(|c| c.is_char()).count()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Character<CharacterData>> {
        self.sequence.iter()
    }

    pub fn iter_characters(&self) -> impl use<'_, CharacterData> + Iterator<Item = char> {
        self.sequence.iter().map(Character::as_char)
    }

    /// Removes blanks from the back of the sequence until the desired length is reached.
    ///
    /// If the desired length is greater than or equal to the current length, then nothing happens.
    /// Panics if any removed character is not a blank.
    pub fn prune_blanks(&mut self, desired_length: usize) {
        while self.len() > desired_length {
            assert_eq!(
                self.sequence.pop().as_ref().map(Character::kind),
                Some(CharacterKind::Blank)
            );
        }
    }

    /// Adds blanks to the back of the sequence until the desired length is reached.
    ///
    /// If the desired length is lower than or equal to the current length, then nothing happens.
    pub fn extend_with_blanks(
        &mut self,
        mut blank_data_generator: impl FnMut() -> CharacterData,
        desired_length: usize,
    ) {
        while self.len() < desired_length {
            self.sequence
                .push(Character::new_blank(blank_data_generator()));
        }
    }

    /// Adds the given characters to the back of the sequence.
    pub fn extend_with(&mut self, extension: impl IntoIterator<Item = Character<CharacterData>>) {
        for character in extension {
            self.sequence.push(character);
        }
    }

    /// Adds the given character to the back of the sequence.
    pub fn push(&mut self, character: Character<CharacterData>) {
        self.sequence.push(character);
    }

    pub fn get(&self, index: usize) -> Option<&Character<CharacterData>> {
        self.sequence.get(index)
    }

    pub fn last(&self) -> Option<&Character<CharacterData>> {
        self.sequence.last()
    }

    pub fn insert_gaps(
        &mut self,
        mut gap_data_generator: impl FnMut() -> CharacterData,
        gaps: impl IntoIterator<Item = usize>,
    ) {
        self.multi_insert(
            iter::repeat_with(|| Character::new_gap(gap_data_generator())),
            gaps,
        );
    }

    pub fn insert_blanks(
        &mut self,
        mut blank_data_generator: impl FnMut() -> CharacterData,
        blanks: impl IntoIterator<Item = usize>,
    ) {
        self.multi_insert(
            iter::repeat_with(|| Character::new_blank(blank_data_generator())),
            blanks,
        );
    }

    pub fn multi_insert(
        &mut self,
        characters: impl IntoIterator<Item = Character<CharacterData>>,
        positions: impl IntoIterator<Item = usize>,
    ) {
        let original_sequence = mem::take(&mut self.sequence);
        let original_sequence_len = original_sequence.len();

        let mut characters = characters.into_iter();
        let mut positions = positions.into_iter().peekable();
        let original_characters = original_sequence.into_iter().enumerate();

        for (index, original_character) in original_characters {
            while let Some(position) = positions.peek().copied() {
                if position <= index {
                    self.sequence.push(characters.next().unwrap());
                    positions.next().unwrap();
                } else {
                    break;
                }
            }

            self.sequence.push(original_character);
        }

        for position in positions {
            debug_assert_eq!(position, original_sequence_len);
            self.sequence.push(characters.next().unwrap());
        }
    }
}

impl<Data> Character<Data> {
    fn new(kind: CharacterKind, data: Data) -> Self {
        Self { kind, data }
    }

    fn new_gap(data: Data) -> Self {
        Self::new(CharacterKind::Gap, data)
    }

    fn new_blank(data: Data) -> Self {
        Self::new(CharacterKind::Blank, data)
    }

    pub fn new_char(character: char, data: Data) -> Self {
        Self::new(CharacterKind::Char(character), data)
    }

    fn kind(&self) -> CharacterKind {
        self.kind
    }

    pub fn data(&self) -> &Data {
        &self.data
    }

    fn is_char(&self) -> bool {
        self.kind.is_char()
    }

    pub fn as_char(&self) -> char {
        self.kind.as_char()
    }

    fn make_ascii_lowercase(&mut self) {
        self.kind.make_ascii_lowercase()
    }
}

impl<Data: Default> Character<Data> {
    pub fn new_char_with_default(character: char) -> Self {
        Self::new_char(character, Default::default())
    }
}

impl CharacterKind {
    fn is_char(&self) -> bool {
        matches!(self, Self::Char(_))
    }

    fn as_char(&self) -> char {
        match self {
            Self::Char(character) => *character,
            Self::Gap => '-',
            Self::Blank => ' ',
        }
    }

    fn make_ascii_lowercase(&mut self) {
        if let Self::Char(character) = self {
            *character = character.to_ascii_lowercase()
        }
    }
}

impl<CharacterData> FromIterator<Character<CharacterData>>
    for MultipairAlignmentSequence<CharacterData>
{
    fn from_iter<T: IntoIterator<Item = Character<CharacterData>>>(iter: T) -> Self {
        Self {
            sequence: iter.into_iter().collect(),
        }
    }
}

impl<CharacterData: Default> FromIterator<CharacterKind>
    for MultipairAlignmentSequence<CharacterData>
{
    fn from_iter<T: IntoIterator<Item = CharacterKind>>(iter: T) -> Self {
        Self {
            sequence: iter
                .into_iter()
                .map(|kind| Character::new(kind, Default::default()))
                .collect(),
        }
    }
}

impl<CharacterData> Default for MultipairAlignmentSequence<CharacterData> {
    fn default() -> Self {
        Self {
            sequence: Default::default(),
        }
    }
}

impl<CharacterData> Display for MultipairAlignmentSequence<CharacterData> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for character in &self.sequence {
            write!(f, "{}", character.as_char())?;
        }

        Ok(())
    }
}

impl<CharacterData> From<Vec<Character<CharacterData>>>
    for MultipairAlignmentSequence<CharacterData>
{
    fn from(value: Vec<Character<CharacterData>>) -> Self {
        Self { sequence: value }
    }
}

impl<CharacterData> Index<usize> for MultipairAlignmentSequence<CharacterData> {
    type Output = <Vec<Character<CharacterData>> as Index<usize>>::Output;

    fn index(&self, index: usize) -> &Self::Output {
        self.sequence.index(index)
    }
}

impl<CharacterData> IndexMut<usize> for MultipairAlignmentSequence<CharacterData> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.sequence.index_mut(index)
    }
}

impl<Data: Display> Display for Character<Data> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.kind, self.data)
    }
}

impl Display for CharacterKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl Display for NoCharacterData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}
