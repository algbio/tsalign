use std::fmt::Display;

use svg::node::element::{Group, Path};

use crate::{svg::font::sans_serif_mono, ts_arrangement::index_types::ArrangementColumn};

use super::{
    SvgLocation,
    font::{DEFAULT_COLOR, svg_phrase, typewriter},
};

const BORDER_STROKE_WIDTH: f32 = 0.2;
const BORDER_PADDING: f32 = 0.4;
const BORDER_RADIUS: f32 = 1.2;

pub struct Number<Row> {
    number: String,
    column: ArrangementColumn,
    row: Row,
    align: NumberAlignment,
    scale: f32,
}

pub enum NumberAlignment {
    Left,
    Right,
}

impl<Row> Number<Row> {
    pub fn new(
        number: impl ToString,
        column: ArrangementColumn,
        row: Row,
        align: NumberAlignment,
        scale: f32,
    ) -> Self {
        Self {
            number: number.to_string(),
            column,
            row,
            align,
            scale,
        }
    }

    pub fn render(&self, mut rows: impl FnMut(&Row) -> f32) -> Group {
        let x = self.column.primitive() as f32 * typewriter::FONT.character_width
            + (BORDER_PADDING + 0.5 * BORDER_STROKE_WIDTH) * self.scale
            - match self.align {
                NumberAlignment::Left => 0.0,
                NumberAlignment::Right => self.width(),
            };
        let y = rows(&self.row);
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        let label_height = sans_serif_mono::FONT.character_height;

        let border_x = BORDER_RADIUS - BORDER_PADDING;
        let border_y = -BORDER_PADDING - label_height;
        let border_h = label_width + 2.0 * (BORDER_PADDING - BORDER_RADIUS);
        let border_v = label_height + 2.0 * (BORDER_PADDING - BORDER_RADIUS);

        let border_d = format!(
            "M {border_x},{border_y} a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 -{BORDER_RADIUS},{BORDER_RADIUS} v {border_v}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 {BORDER_RADIUS},{BORDER_RADIUS} h {border_h}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 {BORDER_RADIUS},-{BORDER_RADIUS} v -{border_v}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 -{BORDER_RADIUS},-{BORDER_RADIUS} h -{border_h}",
        );

        Group::new()
            .set(
                "transform",
                format!("translate({x},{y}) scale({})", self.scale),
            )
            .add(
                Path::new()
                    .set("fill", "none")
                    .set("stroke", "black")
                    .set("stroke-width", BORDER_STROKE_WIDTH)
                    .set("d", border_d),
            )
            .add(Group::new().add(svg_phrase(
                &self.number,
                DEFAULT_COLOR,
                &SvgLocation {
                    x: 0.0,
                    y: sans_serif_mono::FONT.character_height * -0.12,
                },
                &sans_serif_mono::FONT,
            )))
    }

    pub fn width(&self) -> f32 {
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        (label_width + BORDER_STROKE_WIDTH + 2.0 * BORDER_PADDING) * self.scale
    }
}

impl<Row: Display> Display for Number<Row> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} at {}/{}/{}",
            self.number, self.column, self.row, self.align
        )
    }
}

impl Display for NumberAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                NumberAlignment::Left => "L",
                NumberAlignment::Right => "R",
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{Number, NumberAlignment};
    use crate::{
        svg::font::{sans_serif_mono, typewriter},
        ts_arrangement::index_types::ArrangementColumn,
    };

    #[test]
    fn width_scales_with_the_label() {
        let one = Number::new(
            "A",
            ArrangementColumn::new(0usize),
            (),
            NumberAlignment::Left,
            1.0,
        );
        let two = Number::new(
            "A1",
            ArrangementColumn::new(0usize),
            (),
            NumberAlignment::Left,
            1.0,
        );
        let half = Number::new(
            "A1",
            ArrangementColumn::new(0usize),
            (),
            NumberAlignment::Left,
            0.5,
        );
        assert!((two.width() - one.width() - sans_serif_mono::FONT.character_width).abs() < 1e-4);
        assert!((half.width() - 0.5 * two.width()).abs() < 1e-4);
    }

    /// Left aligned numbers start at the column, right aligned numbers end at it.
    #[test]
    fn alignment() {
        let column = ArrangementColumn::new(5usize);
        let left = Number::new("A1", column, (), NumberAlignment::Left, 0.5);
        let right = Number::new("A1", column, (), NumberAlignment::Right, 0.5);
        let x = |number: &Number<()>| -> f32 {
            let rendered = number.render(|_| 7.0).to_string();
            let start = rendered.find("translate(").unwrap() + "translate(".len();
            let end = start + rendered[start..].find(',').unwrap();
            rendered[start..end].parse().unwrap()
        };

        let column_x = 5.0 * typewriter::FONT.character_width;
        assert!(x(&left) >= column_x);
        assert!((x(&left) - x(&right) - right.width()).abs() < 1e-4);
        assert!(left.render(|_| 7.0).to_string().contains(",7) scale(0.5)"));
        assert!(left.render(|_| 7.0).to_string().contains("A1"));
    }

    #[test]
    fn display() {
        let number = Number::new(
            "B4",
            ArrangementColumn::new(3usize),
            "R",
            NumberAlignment::Right,
            1.0,
        );
        assert_eq!(number.to_string(), "B4 at ArrangementColumn(3)/R/R");
    }
}
