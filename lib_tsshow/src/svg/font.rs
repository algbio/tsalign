use std::{
    borrow::Borrow,
    collections::{BTreeSet, HashSet},
    fmt::{Debug, Write},
};

use base64::{Engine, engine::general_purpose::STANDARD as BASE64};
use log::debug;
use svg::node::element::{Group, Style, Text};

use crate::{
    error::{Error, Result},
    plain_text::mutlipair_alignment_renderer::{Character, NoCharacterData},
};

use super::SvgLocation;

pub mod sans_serif;
pub mod sans_serif_mono;
pub mod typewriter;

/// A font that is embedded into the SVG, such that the SVG renders identically everywhere.
///
/// The font file is embedded as an `@font-face` rule (see [`embed_fonts`]),
/// and the characters are rendered as SVG text elements referring to that font family.
#[derive(Debug)]
pub struct Font {
    /// The complete font file, of which a subset is embedded into the SVG.
    file: &'static [u8],
    /// The font family, as declared by the font file.
    family: String,
    /// The copyright and licence notices, as declared by the font file.
    notice: String,
    /// The `font-size` used to render characters of this font.
    font_size: f32,
    /// The horizontal distance between the origins of two adjacent characters.
    pub character_width: f32,
    /// The vertical distance between the baselines of two adjacent lines.
    pub character_height: f32,
    /// The `letter-spacing` that makes the characters of a string exactly
    /// [`character_width`](Self::character_width) apart.
    ///
    /// This is `None` for fonts that are not monospaced.
    letter_spacing: Option<f32>,
}

/// The colour used to render characters without an explicit colour.
pub const DEFAULT_COLOR: Color = Color::opaque("black");

/// A fill colour with an opacity.
///
/// The opacity is kept separate from the colour, because Inkscape does not support
/// `#RRGGBBAA` colours.
#[derive(Debug, Clone, Copy)]
pub struct Color {
    pub fill: &'static str,
    pub fill_opacity: f32,
}

impl Color {
    pub const fn opaque(fill: &'static str) -> Self {
        Self {
            fill,
            fill_opacity: 1.0,
        }
    }

    pub const fn new(fill: &'static str, fill_opacity: f32) -> Self {
        Self { fill, fill_opacity }
    }
}

#[derive(Debug, Clone)]
pub struct CharacterData {
    pub color: Color,
}

impl Font {
    /// The complete font file.
    pub fn file(&self) -> &'static [u8] {
        self.file
    }

    /// The font family, as declared by the font file.
    pub fn family(&self) -> &str {
        &self.family
    }

    /// Creates a monospaced font from the given font file.
    ///
    /// The characters are rendered at a `font-size` of `character_height * scale`,
    /// and are placed `character_width * scale` apart.
    fn new_monospace(
        file: &'static [u8],
        character_width: f32,
        character_height: f32,
        scale: f32,
    ) -> Self {
        let face = parse_font_file(file);
        assert!(
            face.is_monospaced(),
            "Font {:?} is not monospaced",
            font_family(&face)
        );

        let font_size = character_height * scale;
        let character_width = character_width * scale;

        // Compensate for the difference between the advance width of the font
        // and the character width used to lay out the SVG.
        let advance_width = advance_width(&face, font_size);

        Self {
            file,
            family: font_family(&face),
            notice: font_notice(&face),
            font_size,
            character_width,
            character_height: font_size,
            letter_spacing: Some(character_width - advance_width),
        }
    }

    /// Creates a font that is not monospaced from the given font file.
    ///
    /// The width of a string rendered with this font depends on the string itself,
    /// hence [`character_width`](Self::character_width) is merely an estimate.
    fn new_proportional(file: &'static [u8], font_size: f32) -> Self {
        let face = parse_font_file(file);

        Self {
            file,
            family: font_family(&face),
            notice: font_notice(&face),
            font_size,
            character_width: advance_width(&face, font_size),
            character_height: font_size,
            letter_spacing: None,
        }
    }
}

/// Parses a font file that is embedded into this binary, and hence is known to be valid.
fn parse_font_file(file: &'static [u8]) -> ttf_parser::Face<'static> {
    ttf_parser::Face::parse(file, 0).expect("embedded font file is a valid font")
}

/// Returns the font family declared by the font file.
fn font_family(face: &ttf_parser::Face) -> String {
    font_name(face, ttf_parser::name_id::FAMILY).expect("font file declares a font family")
}

/// Returns the copyright and licence notices declared by the font file.
fn font_notice(face: &ttf_parser::Face) -> String {
    [
        ttf_parser::name_id::COPYRIGHT_NOTICE,
        ttf_parser::name_id::LICENSE,
        ttf_parser::name_id::LICENSE_URL,
    ]
    .into_iter()
    .filter_map(|name_id| font_name(face, name_id))
    .collect::<Vec<_>>()
    .join("\n\n")
}

/// Returns the unicode name with the given id declared by the font file.
fn font_name(face: &ttf_parser::Face, name_id: u16) -> Option<String> {
    face.names()
        .into_iter()
        .find(|name| name.name_id == name_id && name.is_unicode())
        .and_then(|name| name.to_string())
}

/// Returns the advance width of a character of the given font at the given font size.
fn advance_width(face: &ttf_parser::Face, font_size: f32) -> f32 {
    let glyph = face
        .glyph_index('0')
        .expect("font file contains the character '0'");
    let advance_width = face
        .glyph_hor_advance(glyph)
        .expect("font file declares advance widths");
    f32::from(advance_width) * font_size / f32::from(face.units_per_em())
}

/// Returns all fonts that may be used in an SVG.
pub fn fonts() -> impl Iterator<Item = &'static Font> {
    [
        &*typewriter::FONT,
        &*sans_serif_mono::FONT,
        &*sans_serif::FONT,
    ]
    .into_iter()
}

/// Returns the printable ASCII characters, which are always embedded into the fonts of an SVG.
pub fn printable_characters() -> impl Iterator<Item = char> {
    ' '..='~'
}

/// Returns a `<style>` element embedding the given fonts, reduced to the given characters.
///
/// Fonts occurring multiple times are embedded only once.
/// Characters not supported by a font are skipped.
pub fn embed_fonts<'font>(
    fonts: impl IntoIterator<Item = &'font Font>,
    characters: &BTreeSet<char>,
) -> Result<Style> {
    let characters: HashSet<char> = characters.iter().copied().collect();
    let mut embedded_families = Vec::new();
    let mut css = String::new();

    for font in fonts {
        if embedded_families.contains(&font.family) {
            continue;
        }
        embedded_families.push(font.family.clone());

        let subset = fontcull::subset_font_data(font.file, &characters, &[])
            .map_err(|error| Error::FontSubsetting(font.family.clone(), error.to_string()))?;
        debug!(
            "Embedding font {} reduced from {} to {} bytes",
            font.family,
            font.file.len(),
            subset.len()
        );

        // CSS comments cannot be nested, hence the notice must not terminate the comment.
        let notice = font.notice.replace("*/", "* /");
        writeln!(css, "/*\n * {}\n *", font.family).unwrap();
        for line in notice.lines() {
            if line.is_empty() {
                writeln!(css, " *").unwrap();
            } else {
                writeln!(css, " * {line}").unwrap();
            }
        }
        writeln!(css, " */").unwrap();
        writeln!(
            css,
            "@font-face {{ font-family: \"{}\"; font-style: normal; font-weight: normal; src: url(\"data:font/ttf;base64,{}\") format(\"truetype\"); }}",
            font.family,
            BASE64.encode(&subset),
        )
        .unwrap();
    }

    Ok(Style::new(css).set("type", "text/css"))
}

/// Renders a sequence of characters as one SVG text element per character.
///
/// This places the characters in a strict grid, and allows to colour each character individually.
pub fn svg_string<
    'character,
    Data: 'character + Debug,
    Item: 'character + Borrow<Character<Data>>,
>(
    string: impl IntoIterator<Item = Item>,
    location: &SvgLocation,
    font: &Font,
) -> Group
where
    CharacterData: for<'a> From<&'a Data>,
{
    let mut group = Group::new()
        .set("transform", location.as_transform())
        .set("font-family", font.family.clone())
        .set("font-size", font.font_size);

    for (index, character) in string.into_iter().enumerate() {
        let character = character.borrow();
        let c = character.as_char();

        // Whitespace renders as nothing, so it can be skipped altogether.
        if c.is_whitespace() {
            continue;
        }

        let data: CharacterData = character.data().into();
        let color = data.color;
        let mut text = Text::new(c.to_string())
            .set("x", index as f32 * font.character_width)
            .set("fill", color.fill);
        if color.fill_opacity < 1.0 {
            text = text.set("fill-opacity", color.fill_opacity);
        }
        group = group.add(text);
    }

    group
}

/// Renders a string as a single SVG text element.
///
/// If the font is monospaced, then the characters are placed
/// [`character_width`](Font::character_width) apart.
pub fn svg_phrase(
    string: impl AsRef<str>,
    color: Color,
    location: &SvgLocation,
    font: &Font,
) -> Text {
    let mut text = Text::new(string.as_ref())
        .set("x", location.x)
        .set("y", location.y)
        .set("font-family", font.family.clone())
        .set("font-size", font.font_size)
        .set("fill", color.fill);

    if color.fill_opacity < 1.0 {
        text = text.set("fill-opacity", color.fill_opacity);
    }

    if let Some(letter_spacing) = font.letter_spacing {
        text = text.set("letter-spacing", letter_spacing);
    }

    text
}

/// Renders text with a font that is not monospaced.
///
/// The width of the text depends on the text itself, and hence is not known exactly.
pub fn svg_text(text: impl AsRef<str>, location: &SvgLocation) -> Text {
    svg_phrase(text, DEFAULT_COLOR, location, &sans_serif::FONT)
}

impl CharacterData {
    pub fn new_colored(color: Color) -> Self {
        #[allow(clippy::needless_update)]
        Self { color }
    }
}

impl Default for CharacterData {
    fn default() -> Self {
        Self {
            color: DEFAULT_COLOR,
        }
    }
}

impl<'a> From<&'a NoCharacterData> for CharacterData {
    fn from(_: &'a NoCharacterData) -> Self {
        Self::default()
    }
}

impl<'a> From<&'a CharacterData> for CharacterData {
    fn from(value: &'a CharacterData) -> Self {
        value.clone()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::{embed_fonts, fonts, printable_characters};

    /// The characters embedded into a font must be rendered exactly like by the original font.
    #[test]
    fn test_embedded_fonts_contain_the_printable_characters() {
        let characters: BTreeSet<char> = printable_characters().collect();

        for font in fonts() {
            let subset =
                fontcull::subset_font_data(font.file, &characters.iter().copied().collect(), &[])
                    .unwrap();
            let original = ttf_parser::Face::parse(font.file, 0).unwrap();
            let embedded = ttf_parser::Face::parse(&subset, 0).unwrap();
            assert_eq!(original.units_per_em(), embedded.units_per_em());

            for character in &characters {
                let original_glyph = original.glyph_index(*character).unwrap_or_else(|| {
                    panic!("Font {} does not support {character:?}", font.family)
                });
                let embedded_glyph = embedded.glyph_index(*character).unwrap_or_else(|| {
                    panic!(
                        "Font {} was not embedded with support for {character:?}",
                        font.family
                    )
                });
                assert_eq!(
                    original.glyph_hor_advance(original_glyph),
                    embedded.glyph_hor_advance(embedded_glyph),
                    "Font {} renders {character:?} with a different advance width when embedded",
                    font.family,
                );
            }
        }
    }

    /// Each font must be embedded exactly once, under the family used to render it.
    #[test]
    fn test_embed_fonts_declares_each_font_once() {
        let characters: BTreeSet<char> = printable_characters().collect();
        let style = embed_fonts(fonts().chain(fonts()), &characters)
            .unwrap()
            .to_string();

        assert_eq!(style.matches("@font-face").count(), fonts().count());
        for font in fonts() {
            assert_eq!(
                style
                    .matches(&format!("font-family: \"{}\";", font.family))
                    .count(),
                1,
            );
        }

        // The copyright notices must not terminate the comment they are contained in.
        assert_eq!(style.matches("/*").count(), style.matches("*/").count());
    }
}
