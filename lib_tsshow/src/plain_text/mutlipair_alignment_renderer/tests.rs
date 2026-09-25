use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;

use crate::plain_text::mutlipair_alignment_renderer::{
    Character, CharacterKind, MultipairAlignmentSequence, NoCharacterData,
};

use super::MultipairAlignmentRenderer;

#[test]
fn parallel_gaps() {
    let mut renderer = MultipairAlignmentRenderer::new(
        "B".to_string(),
        "GGG"
            .chars()
            .map(|c| Character::new_char(c, NoCharacterData)),
    );
    renderer.add_aligned_sequence(
        &"B".to_string(),
        0,
        "A".to_string(),
        "GGGG"
            .chars()
            .map(|c| Character::new_char(c, NoCharacterData)),
        || NoCharacterData,
        || NoCharacterData,
        [
            AlignmentType::PrimaryMatch,
            AlignmentType::PrimaryMatch,
            AlignmentType::PrimaryInsertion,
            AlignmentType::PrimaryMatch,
        ],
        true,
        false,
    );
    renderer.add_aligned_sequence(
        &"B".to_string(),
        0,
        "C".to_string(),
        "GGGG"
            .chars()
            .map(|c| Character::new_char(c, NoCharacterData)),
        || NoCharacterData,
        || NoCharacterData,
        [
            AlignmentType::PrimaryMatch,
            AlignmentType::PrimaryMatch,
            AlignmentType::PrimaryInsertion,
            AlignmentType::PrimaryMatch,
        ],
        true,
        false,
    );

    let mut output = Vec::new();
    renderer
        .render(
            &mut output,
            [&"A".to_string(), &"B".to_string(), &"C".to_string()],
        )
        .unwrap();

    let output = String::from_utf8(output).unwrap();
    assert_eq!(output.trim(), "A: GGGG\nB: GG-G\nC: GGGG");
}

#[test]
fn translate_alignment_offset() {
    let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([
        CharacterKind::Blank,
        CharacterKind::Char('A'),
        CharacterKind::Gap,
        CharacterKind::Char('C'),
        CharacterKind::Gap,
    ]);
    assert_eq!(sequence.translate_alignment_offset(0), Some(0));
    assert_eq!(sequence.translate_alignment_offset(1), Some(2));
    assert_eq!(sequence.translate_alignment_offset(2), Some(4));
    assert_eq!(sequence.translate_alignment_offset(3), None);
}

#[test]
fn translate_extension_offset() {
    let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([
        CharacterKind::Blank,
        CharacterKind::Char('A'),
        CharacterKind::Gap,
        CharacterKind::Char('C'),
        CharacterKind::Gap,
    ]);
    assert_eq!(sequence.translate_extension_offset(0), Some(1));
    assert_eq!(sequence.translate_extension_offset(1), Some(3));
    assert_eq!(sequence.translate_extension_offset(2), Some(5));
    assert_eq!(sequence.translate_alignment_offset(3), None);

    let sequence = MultipairAlignmentSequence::<NoCharacterData>::from_iter([
        CharacterKind::Blank,
        CharacterKind::Char('A'),
        CharacterKind::Gap,
        CharacterKind::Char('C'),
    ]);
    assert_eq!(sequence.translate_extension_offset(0), Some(1));
    assert_eq!(sequence.translate_extension_offset(1), Some(3));
    assert_eq!(sequence.translate_extension_offset(2), Some(4));
    assert_eq!(sequence.translate_alignment_offset(3), None);
}
