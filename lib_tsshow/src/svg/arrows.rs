use std::fmt::Display;

use log::debug;
use svg::{
    Document,
    node::element::{Definitions, Group, Marker, Path},
};

use crate::{svg::font::typewriter, ts_arrangement::index_types::ArrangementColumn};

#[derive(Debug, PartialEq)]
pub struct Arrow<Row> {
    pub from: ArrowEndpoint<Row>,
    pub to: ArrowEndpoint<Row>,
    pub style: ArrowStyle,
}

#[derive(Debug, PartialEq)]
pub struct ArrowEndpoint<Row> {
    /// The column the endpoint is at, which will be translated to an `x`-coordinate.
    ///
    /// The column is interpreted as the left side of the bounding box of a character.
    pub column: ArrangementColumn,

    /// The row the endpoint is at, which will be translated to a `y`-coordinate.
    ///
    /// The row is interpreted as the horizontal center of the bounding box of a character.
    pub row: Row,

    /// The direction from which the arrow travels to this endpoint.
    ///
    /// **Note:** this has nothing to do with where the head and tail of the arrow are,
    /// but it decides where the curve the arrow is made of goes.
    pub direction: ArrowEndpointDirection,

    /// Extra offset affecting the `y`-coordinate of this endpoint.
    ///
    /// The offset is applied according to [`direction`](Self::direction).
    pub column_offset: f32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ArrowEndpointDirection {
    Forward,
    Backward,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ArrowStyle {
    Direct,
    Curved,
}

impl<Row> Arrow<Row> {
    #[expect(dead_code)]
    pub fn new_skip(from_column: ArrangementColumn, to_column: ArrangementColumn, row: Row) -> Self
    where
        Row: Clone,
    {
        Self {
            from: ArrowEndpoint {
                column: from_column,
                row: row.clone(),
                direction: if from_column <= to_column {
                    ArrowEndpointDirection::Forward
                } else {
                    ArrowEndpointDirection::Backward
                },
                column_offset: 0.0,
            },
            to: ArrowEndpoint {
                column: to_column,
                row,
                direction: if from_column <= to_column {
                    ArrowEndpointDirection::Backward
                } else {
                    ArrowEndpointDirection::Forward
                },
                column_offset: 0.0,
            },
            style: ArrowStyle::Direct,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_curved(
        from_column: ArrangementColumn,
        from_column_offset: f32,
        from_row: Row,
        from_direction: ArrowEndpointDirection,
        to_column: ArrangementColumn,
        to_column_offset: f32,
        to_row: Row,
        to_direction: ArrowEndpointDirection,
    ) -> Self {
        Self {
            from: ArrowEndpoint::new(from_column, from_row, from_direction, from_column_offset),
            to: ArrowEndpoint::new(to_column, to_row, to_direction, to_column_offset),
            style: ArrowStyle::Curved,
        }
    }

    pub fn render(&self, mut rows: impl FnMut(&Row) -> f32) -> Group
    where
        Row: Display,
    {
        match self.style {
            ArrowStyle::Direct => {
                debug!("Drawing direct arrow from {} to {}", self.from, self.to);

                let stroke_width = 0.15;
                let from_x = self.from.column.primitive() as f32 * typewriter::FONT.character_width
                    + self.from.column_offset * self.from.direction.sign();
                let from_y = rows(&self.from.row) - typewriter::FONT.character_height * 0.3;
                let to_x = self.to.column.primitive() as f32 * typewriter::FONT.character_width
                    + self.to.column_offset * self.to.direction.sign();
                let to_y = rows(&self.to.row) - typewriter::FONT.character_height * 0.3;

                Group::new().add(
                    Path::new()
                        .set("marker-end", "url(#arrow_head)")
                        .set("stroke-width", stroke_width)
                        .set("stroke", "black")
                        .set("fill", "none")
                        .set("stroke-dasharray", stroke_width)
                        .set("d", format!("M {from_x},{from_y} L {to_x},{to_y}")),
                )
            }

            ArrowStyle::Curved => {
                debug!("Drawing curved arrow from {} to {}", self.from, self.to);

                let stroke_width = 0.15;
                let from_x = self.from.column.primitive() as f32 * typewriter::FONT.character_width
                    + self.from.column_offset * self.from.direction.sign();
                let from_y = rows(&self.from.row) - typewriter::FONT.character_height * 0.3;
                let to_x = self.to.column.primitive() as f32 * typewriter::FONT.character_width
                    + self.to.column_offset * self.to.direction.sign();
                let to_y = rows(&self.to.row) - typewriter::FONT.character_height * 0.3;

                let (from_x_control, to_x_control) =
                    match (&self.from.direction, &self.to.direction) {
                        (ArrowEndpointDirection::Forward, ArrowEndpointDirection::Forward) => {
                            let x_control =
                                from_x.max(to_x) + 2.0 * typewriter::FONT.character_width;
                            (x_control, x_control)
                        }
                        (ArrowEndpointDirection::Backward, ArrowEndpointDirection::Backward) => {
                            let x_control =
                                from_x.min(to_x) - 2.0 * typewriter::FONT.character_width;
                            (x_control, x_control)
                        }
                        (ArrowEndpointDirection::Forward, ArrowEndpointDirection::Backward)
                        | (ArrowEndpointDirection::Backward, ArrowEndpointDirection::Forward) => {
                            let x_control_delta = ((from_x - to_x).abs() * 0.1)
                                .max(2.0 * typewriter::FONT.character_width);
                            (
                                from_x + self.from.direction.sign() * x_control_delta,
                                to_x + self.to.direction.sign() * x_control_delta,
                            )
                        }
                    };

                Group::new().add(
                    Path::new()
                        .set("marker-end", "url(#arrow_head_red)")
                        .set("stroke-width", stroke_width)
                        .set("stroke", "#CE2029")
                        .set("fill", "none")
                        .set("stroke-dasharray", stroke_width)
                        .set("d", format!("M {from_x},{from_y} C {from_x_control},{from_y} {to_x_control},{to_y} {to_x},{to_y}")),
                )
            }
        }
    }
}

impl<Row> ArrowEndpoint<Row> {
    pub fn new(
        column: ArrangementColumn,
        row: Row,
        direction: ArrowEndpointDirection,
        column_offset: f32,
    ) -> Self {
        Self {
            column,
            row,
            direction,
            column_offset,
        }
    }
}

impl ArrowEndpointDirection {
    /// Returns `1.0` if the direction is forward, and `-1.0` if the direction is backward.
    pub fn sign(&self) -> f32 {
        match self {
            ArrowEndpointDirection::Forward => 1.0,
            ArrowEndpointDirection::Backward => -1.0,
        }
    }
}

pub fn add_arrow_defs(svg: Document) -> Document {
    svg.add(
        Definitions::new()
            .add(
                Marker::new()
                    .set("id", "arrow_head")
                    .set("viewBox", "0 0 10 10")
                    .set("orient", "auto-start-reverse")
                    .set("markerWidth", 10)
                    .set("markerHeight", 10)
                    .set("refX", 10)
                    .set("refY", 5)
                    .add(
                        Path::new()
                            .set("d", "M 1 1 L 10 5 L 1 9")
                            .set("fill", "none")
                            .set("stroke", "black"),
                    ),
            )
            .add(
                Marker::new()
                    .set("id", "arrow_head_red")
                    .set("viewBox", "0 0 10 10")
                    .set("orient", "auto-start-reverse")
                    .set("markerWidth", 10)
                    .set("markerHeight", 10)
                    .set("refX", 10)
                    .set("refY", 5)
                    .add(
                        Path::new()
                            .set("d", "M 1 1 L 10 5 L 1 9")
                            .set("fill", "none")
                            .set("stroke", "#CE2029"),
                    ),
            ),
    )
}

impl<Row: Display> Display for Arrow<Row> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}({} -> {})",
            match self.style {
                ArrowStyle::Direct => "D",
                ArrowStyle::Curved => "C",
            },
            self.from,
            self.to
        )
    }
}

impl<Row: Display> Display for ArrowEndpoint<Row> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}/{}", self.column, self.row, self.direction)
    }
}

impl Display for ArrowEndpointDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ArrowEndpointDirection::Forward => "+",
                ArrowEndpointDirection::Backward => "-",
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use svg::Document;

    use super::{Arrow, ArrowEndpointDirection, add_arrow_defs};
    use crate::{svg::font::typewriter, ts_arrangement::index_types::ArrangementColumn};

    fn path(arrow: &Arrow<&str>) -> Vec<f32> {
        let rendered = arrow
            .render(|row| if *row == "top" { 10.0 } else { 20.0 })
            .to_string();
        let start = rendered.find(" d=\"").unwrap() + 4;
        let end = start + rendered[start..].find('"').unwrap();
        rendered[start..end]
            .split([' ', ','])
            .filter_map(|n| n.parse().ok())
            .collect()
    }

    #[test]
    fn sign() {
        assert_eq!(ArrowEndpointDirection::Forward.sign(), 1.0);
        assert_eq!(ArrowEndpointDirection::Backward.sign(), -1.0);
    }

    /// A curved arrow starts and ends at its columns, moved by the column offsets in the endpoint directions.
    #[test]
    fn curved_arrow_endpoints() {
        let arrow = Arrow::new_curved(
            ArrangementColumn::new(2usize),
            1.0,
            "top",
            ArrowEndpointDirection::Forward,
            ArrangementColumn::new(6usize),
            0.5,
            "bottom",
            ArrowEndpointDirection::Backward,
        );
        let width = typewriter::FONT.character_width;
        let height = typewriter::FONT.character_height;
        let numbers = path(&arrow);
        // M from C control1 control2 to
        assert_eq!(numbers.len(), 8);
        assert!((numbers[0] - (2.0 * width + 1.0)).abs() < 1e-4);
        assert!((numbers[1] - (10.0 - 0.3 * height)).abs() < 1e-4);
        assert!((numbers[6] - (6.0 * width - 0.5)).abs() < 1e-4);
        assert!((numbers[7] - (20.0 - 0.3 * height)).abs() < 1e-4);
        // Forward from, backward to: the curve leaves to the right and arrives from the left.
        assert!(numbers[2] > numbers[0]);
        assert!(numbers[4] < numbers[6]);
        assert_eq!(
            arrow.to_string(),
            "C(ArrangementColumn(2)/top/+ -> ArrangementColumn(6)/bottom/-)"
        );
    }

    /// Arrows with both endpoints in the same direction bend around the outermost endpoint.
    #[test]
    fn curved_arrow_same_directions() {
        for (direction, outer) in [
            (
                ArrowEndpointDirection::Forward,
                f32::max as fn(f32, f32) -> f32,
            ),
            (ArrowEndpointDirection::Backward, f32::min),
        ] {
            let arrow = Arrow::new_curved(
                ArrangementColumn::new(2usize),
                0.0,
                "top",
                direction,
                ArrangementColumn::new(6usize),
                0.0,
                "bottom",
                direction,
            );
            let numbers = path(&arrow);
            let control = outer(numbers[0], numbers[6])
                + direction.sign() * 2.0 * typewriter::FONT.character_width;
            assert!((numbers[2] - control).abs() < 1e-4);
            assert!((numbers[4] - control).abs() < 1e-4);
        }
    }

    #[test]
    fn arrow_heads_are_defined() {
        let document = add_arrow_defs(Document::new()).to_string();
        assert!(document.contains("id=\"arrow_head\""));
        assert!(document.contains("id=\"arrow_head_red\""));

        let arrow = Arrow::new_curved(
            ArrangementColumn::new(0usize),
            0.0,
            "top",
            ArrowEndpointDirection::Forward,
            ArrangementColumn::new(1usize),
            0.0,
            "top",
            ArrowEndpointDirection::Backward,
        );
        assert!(
            arrow
                .render(|_| 0.0)
                .to_string()
                .contains("url(#arrow_head_red)")
        );
    }
}
