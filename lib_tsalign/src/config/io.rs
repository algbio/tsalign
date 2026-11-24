use std::str::FromStr;

use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;
use log::{trace, warn};
use nom::{
    IResult, Parser,
    bytes::complete::{tag, take_while1},
    character::complete::{digit1, line_ending},
    combinator::opt,
};
use num_traits::{Bounded, PrimInt};

use crate::{
    config::BaseCost,
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    io::{
        parse_any_whitespace, parse_title, parse_whitespace, skip_any_whitespace,
        translate_nom_error,
    },
};

use super::TemplateSwitchConfig;

impl<AlphabetType: Alphabet, Cost: AStarCost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn read_plain(mut reader: impl std::io::Read) -> crate::error::Result<Self> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;
        input.parse()
    }

    /// Parses a config and checks its legacy flank sections.
    fn parse_and_verify(input: &str) -> crate::error::Result<Self> {
        let (result, legacy_flanks) = Self::parse_plain(input)
            .map(|(_, result)| result)
            .map_err(translate_nom_error)?;
        legacy_flanks.check()?;
        result.verify()?;
        Ok(result)
    }

    fn parse_plain(input: &str) -> IResult<&str, (Self, LegacyFlanks)> {
        trace!("Parsing legacy limits");
        let (input, limits) = opt(parse_specific_name("Limits")).parse(input)?;
        let (input, flank_lengths) = if limits.is_some() {
            let (input, left_flank_length) =
                parse_specific_equals_value("left_flank_length")(input)?;
            let (input, right_flank_length) =
                parse_specific_equals_value("right_flank_length")(input)?;
            (input, Some((left_flank_length, right_flank_length)))
        } else {
            (input, None)
        };

        trace!("Parsing base costs");
        let (input, ()) = parse_specific_name("Base Cost")(input)?;
        let (input, rrf_cost) = parse_specific_equals_value("rrf_cost")(input)?;
        let (input, rqf_cost) = parse_specific_equals_value("rqf_cost")(input)?;
        let (input, qrf_cost) = parse_specific_equals_value("qrf_cost")(input)?;
        let (input, qqf_cost) = parse_specific_equals_value("qqf_cost")(input)?;
        let (input, rrr_cost) = parse_specific_equals_value("rrr_cost")(input)?;
        let (input, rqr_cost) = parse_specific_equals_value("rqr_cost")(input)?;
        let (input, qrr_cost) = parse_specific_equals_value("qrr_cost")(input)?;
        let (input, qqr_cost) = parse_specific_equals_value("qqr_cost")(input)?;

        trace!("Parsing jump costs");
        let (input, ()) = parse_specific_name("Jump Costs")(input)?;
        let (input, rq_qr_offset_costs) = parse_named_cost_function("RQQROffset")(input)?;
        let (input, rr_qq_offset_costs) = parse_named_cost_function("RRQQOffset")(input)?;
        let (input, length_costs) = parse_named_cost_function("Length")(input)?;
        let (input, length_difference_costs) =
            parse_named_cost_function("LengthDifference")(input)?;
        let (input, forward_anti_descendant_gap_costs) =
            parse_named_cost_function("ForwardAntiDescendantGap")(input)?;
        let (input, reverse_anti_descendant_gap_costs) =
            parse_named_cost_function("ReverseAntiDescendantGap")(input)?;

        trace!("Parsing primary edit costs");
        let (input, primary_edit_costs) = parse_named_cost_table("Primary Edit Costs")(input)?;
        trace!("Parsing secondary forward edit costs");
        let (input, secondary_forward_edit_costs) =
            parse_named_cost_table("Secondary Forward Edit Costs")(input)?;
        trace!("Parsing secondary reverse edit costs");
        let (input, secondary_reverse_edit_costs) =
            parse_named_cost_table("Secondary Reverse Edit Costs")(input)?;

        trace!("Parsing legacy flank edit costs");
        let mut input = input;
        let mut flank_edit_cost_tables = false;
        for name in ["Left Flank Edit Costs", "Right Flank Edit Costs"] {
            if parse_title(input).is_ok_and(|(_, title)| title == name) {
                (input, _) = parse_named_cost_table::<AlphabetType, Cost>(name)(input)?;
                flank_edit_cost_tables = true;
            }
        }

        Ok((
            input,
            (
                Self {
                    template_switch_min_length: length_costs
                        .minimum_finite_input()
                        .unwrap_or(usize::MAX),

                    base_cost: BaseCost {
                        rrf: rrf_cost,
                        rqf: rqf_cost,
                        qrf: qrf_cost,
                        qqf: qqf_cost,
                        rrr: rrr_cost,
                        rqr: rqr_cost,
                        qrr: qrr_cost,
                        qqr: qqr_cost,
                    },

                    primary_edit_costs,
                    secondary_forward_edit_costs,
                    secondary_reverse_edit_costs,

                    rq_qr_offset_costs,
                    rr_qq_offset_costs,
                    length_costs,
                    length_difference_costs,
                    forward_anti_descendant_gap_costs,
                    reverse_anti_descendant_gap_costs,
                },
                LegacyFlanks {
                    flank_lengths,
                    flank_edit_cost_tables,
                },
            ),
        ))
    }
}

/// The flank sections of a config file from before flanks were removed (#275).
///
/// Flanks were never used: all known config files set both flank lengths to zero,
/// which makes the flank edit costs irrelevant.
/// Such files are still accepted, but the flank sections are ignored with a warning.
struct LegacyFlanks {
    /// The values of `left_flank_length` and `right_flank_length` in the `# Limits` section, if present.
    flank_lengths: Option<(isize, isize)>,
    /// True if the file contains a `Left Flank Edit Costs` or `Right Flank Edit Costs` table.
    flank_edit_cost_tables: bool,
}

impl LegacyFlanks {
    /// Returns an error if the flank lengths are non-zero, since then the flanks would change the alignment.
    /// Warns if the file contains zero-length flank sections, which are ignored.
    fn check(&self) -> crate::error::Result<()> {
        match self.flank_lengths {
            Some((left, right)) if left != 0 || right != 0 => {
                return Err(crate::error::Error::NonZeroFlankLength { left, right });
            }
            _ => {}
        }

        if self.flank_lengths.is_some() || self.flank_edit_cost_tables {
            warn!(
                "Ignoring the flank sections ('# Limits', flank edit costs) of the template switch cost file: \
                 flanks are not supported anymore. Remove these sections to silence this warning."
            );
        }

        Ok(())
    }
}

fn parse_specific_name(name: &str) -> impl '_ + FnMut(&str) -> IResult<&str, ()> {
    move |input| {
        (
            parse_any_whitespace,
            tag("#"),
            parse_whitespace,
            tag(name),
            parse_whitespace,
            line_ending,
            parse_any_whitespace,
        )
            .parse(input)
            .map(|(input, _)| (input, ()))
    }
}

fn parse_specific_equals_value<Value: FromStr + Bounded>(
    identifier: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, Value> {
    move |input| {
        let (input, (actual_identifier, value)) = parse_equals_value(input)?;
        if actual_identifier == identifier {
            Ok((input, value))
        } else {
            Err(nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            }))
        }
    }
}

fn parse_equals_value<Value: FromStr + Bounded>(input: &str) -> IResult<&str, (&str, Value)> {
    let input = skip_any_whitespace(input)?;
    let (input, identifier) = take_while1(|c: char| c.is_alphanumeric() || c == '_')(input)?;
    let (input, _) = (parse_whitespace, tag("="), parse_whitespace).parse(input)?;
    let (input, value) = parse_inf_value(input)?;

    Ok((input, (identifier, value)))
}

fn parse_named_cost_function<SourceType: FromStr + PrimInt, Cost: FromStr + Bounded>(
    name: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, CostFunction<SourceType, Cost>> {
    move |input| {
        let input = skip_any_whitespace(input)?;
        let (input, _) = tag(name)(input)?;
        CostFunction::parse_plain(input)
    }
}

fn parse_named_cost_table<AlphabetType: Alphabet, Cost: AStarCost>(
    name: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, GapAffineAlignmentCostTable<AlphabetType, Cost>> {
    move |input| {
        let (input, result) = GapAffineAlignmentCostTable::parse_plain(input)?;
        if result.name() == name {
            Ok((input, result))
        } else {
            Err(nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            }))
        }
    }
}

pub fn parse_inf_value<Output: FromStr + Bounded>(input: &str) -> IResult<&str, Output> {
    let mut length = 0;

    let negative = match input
        .chars()
        .next()
        .ok_or(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))? {
        '-' => {
            length += 1;
            true
        }
        '+' => {
            length += 1;
            false
        }
        _ => false,
    };

    if input[length..].starts_with("inf") {
        length += 3;

        if negative {
            Ok((&input[length..], Output::min_value()))
        } else {
            Ok((&input[length..], Output::max_value()))
        }
    } else {
        length += digit1(&input[length..])?.1.len();
        let result = Output::from_str(&input[..length]).map_err(|_| {
            nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            })
        })?;

        Ok((&input[length..], result))
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost> std::fmt::Display
    for TemplateSwitchConfig<AlphabetType, Cost>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "# Base Cost")?;
        writeln!(f, "rrf_cost = {}", self.base_cost.rrf)?;
        writeln!(f, "rqf_cost = {}", self.base_cost.rqf)?;
        writeln!(f, "qrf_cost = {}", self.base_cost.qrf)?;
        writeln!(f, "qqf_cost = {}", self.base_cost.qqf)?;
        writeln!(f, "rrr_cost = {}", self.base_cost.rrr)?;
        writeln!(f, "rqr_cost = {}", self.base_cost.rqr)?;
        writeln!(f, "qrr_cost = {}", self.base_cost.qrr)?;
        writeln!(f, "qqr_cost = {}", self.base_cost.qqr)?;

        writeln!(f, "# Jump Costs")?;
        writeln!(f, "RQQROffset")?;
        writeln!(f, "{}", self.rq_qr_offset_costs)?;
        writeln!(f, "RRQQOffset")?;
        writeln!(f, "{}", self.rr_qq_offset_costs)?;
        writeln!(f, "Length")?;
        writeln!(f, "{}", self.length_costs)?;
        writeln!(f, "LengthDifference")?;
        writeln!(f, "{}", self.length_difference_costs)?;
        writeln!(f, "ForwardAntiDescendantGap")?;
        writeln!(f, "{}", self.forward_anti_descendant_gap_costs)?;
        writeln!(f, "ReverseAntiDescendantGap")?;
        writeln!(f, "{}", self.reverse_anti_descendant_gap_costs)?;

        writeln!(f, "{}", self.primary_edit_costs)?;
        writeln!(f, "{}", self.secondary_forward_edit_costs)?;
        writeln!(f, "{}", self.secondary_reverse_edit_costs)?;

        Ok(())
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost> FromStr for TemplateSwitchConfig<AlphabetType, Cost> {
    type Err = crate::error::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse_and_verify(s)
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::alphabets::{dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN},
        interface::alphabet::Alphabet,
    };
    use generic_a_star::cost::U64Cost;

    use crate::{config::TemplateSwitchConfig, error::Error};

    #[test]
    fn display() {
        let config = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default();
        let string = config.to_string();
        assert_eq!(
            config,
            TemplateSwitchConfig::read_plain(string.as_bytes()).unwrap()
        );
    }

    const LEGACY_FLANK_EDIT_COSTS: &str = "
# Left Flank Edit Costs

SubstitutionCostTable
  |  A  C  G  T  N
--+---------------
A |  0  3  3  3  0
C |  3  0  3  3  0
G |  3  3  0  3  0
T |  3  3  3  0  0
N |  0  0  0  0  0

GapOpenCostVector
 A C G T N
 4 4 4 4 4

GapExtendCostVector
 A C G T N
 1 1 1 1 1

# Right Flank Edit Costs

SubstitutionCostTable
  |  A  C  G  T  N
--+---------------
A |  0  3  3  3  0
C |  3  0  3  3  0
G |  3  3  0  3  0
T |  3  3  3  0  0
N |  0  0  0  0  0

GapOpenCostVector
 A C G T N
 4 4 4 4 4

GapExtendCostVector
 A C G T N
 1 1 1 1 1
";

    /// A config in the format from before flanks were removed (#275).
    fn legacy_config(left_flank_length: isize, right_flank_length: isize) -> String {
        let config = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default();
        format!(
            "# Limits\n\nleft_flank_length = {left_flank_length}\nright_flank_length = {right_flank_length}\n\n{config}{LEGACY_FLANK_EDIT_COSTS}"
        )
    }

    #[test]
    fn legacy_zero_flanks_are_ignored() {
        let expected = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default();
        assert_eq!(
            expected,
            TemplateSwitchConfig::read_plain(legacy_config(0, 0).as_bytes()).unwrap()
        );

        // Only one of the two legacy sections.
        let limits_only = legacy_config(0, 0).replace(LEGACY_FLANK_EDIT_COSTS, "");
        assert_eq!(expected, limits_only.parse().unwrap());
        let tables_only = format!("{expected}{LEGACY_FLANK_EDIT_COSTS}");
        assert_eq!(expected, tables_only.parse().unwrap());
    }

    #[test]
    fn legacy_non_zero_flanks_are_rejected() {
        for (left, right) in [(5, 5), (1, 0), (0, 1)] {
            let result = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::read_plain(
                legacy_config(left, right).as_bytes(),
            );
            assert!(
                matches!(
                    result,
                    Err(Error::NonZeroFlankLength { left: l, right: r }) if l == left && r == right
                ),
                "{result:?}"
            );
        }
    }

    #[test]
    fn repository_configs_have_no_flanks() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let mut paths = vec![root.join("sample_tsa_config/config.tsa")];
        for entry in std::fs::read_dir(root.join("test_files/config")).unwrap() {
            paths.push(entry.unwrap().path().join("config.tsa"));
        }

        for path in paths {
            let input = std::fs::read_to_string(&path).unwrap();
            assert!(
                !input.to_lowercase().contains("flank"),
                "{} still contains flanks",
                path.display()
            );
            // `indel` uses an outdated format of the base costs (`rr_cost`, ...), which could not be parsed before #275 either.
            if path.ends_with("indel/config.tsa") {
                continue;
            }

            // Configs either use the DNA alphabet with or without N.
            if input.contains(" N\n") {
                assert_parses_and_round_trips::<DnaAlphabetOrN>(&input, &path);
            } else {
                assert_parses_and_round_trips::<DnaAlphabet>(&input, &path);
            }
        }
    }

    fn assert_parses_and_round_trips<AlphabetType: Alphabet + Eq + std::fmt::Debug>(
        input: &str,
        path: &std::path::Path,
    ) {
        let config: TemplateSwitchConfig<AlphabetType, U64Cost> = input
            .parse()
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        assert_eq!(config, config.to_string().parse().unwrap());
    }
}
