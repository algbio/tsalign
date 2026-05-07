use std::str::FromStr;

use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;
use log::trace;
use nom::{
    IResult, Parser,
    bytes::complete::{tag, take_while1},
    character::complete::{digit1, line_ending},
};
use num_traits::{Bounded, PrimInt};

use crate::{
    config::BaseCost,
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    io::{parse_any_whitespace, parse_whitespace, skip_any_whitespace, translate_nom_error},
};

use super::TemplateSwitchConfig;

impl<AlphabetType: Alphabet, Cost: AStarCost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn read_plain(mut reader: impl std::io::Read) -> crate::error::Result<Self> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;
        let input = input.as_str();
        let result = Self::parse_plain(input)
            .map(|(_, result)| result)
            .map_err(translate_nom_error)?;
        result.verify()?;
        Ok(result)
    }

    fn parse_plain(input: &str) -> IResult<&str, Self> {
        trace!("Parsing limits");
        let (input, ()) = parse_specific_name("Limits")(input)?;
        let (input, left_flank_length) = parse_specific_equals_value("left_flank_length")(input)?;
        let (input, right_flank_length) = parse_specific_equals_value("right_flank_length")(input)?;

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
        let (input, forward_anti_primary_gap_costs) =
            parse_named_cost_function("ForwardAntiPrimaryGap")(input)?;
        let (input, reverse_anti_primary_gap_costs) =
            parse_named_cost_function("ReverseAntiPrimaryGap")(input)?;

        trace!("Parsing primary edit costs");
        let (input, primary_edit_costs) = parse_named_cost_table("Primary Edit Costs")(input)?;
        trace!("Parsing secondary forward edit costs");
        let (input, secondary_forward_edit_costs) =
            parse_named_cost_table("Secondary Forward Edit Costs")(input)?;
        trace!("Parsing secondary reverse edit costs");
        let (input, secondary_reverse_edit_costs) =
            parse_named_cost_table("Secondary Reverse Edit Costs")(input)?;
        trace!("Parsing left flank edit costs");
        let (input, left_flank_edit_costs) =
            parse_named_cost_table("Left Flank Edit Costs")(input)?;
        trace!("Parsing right flank edit costs");
        let (input, right_flank_edit_costs) =
            parse_named_cost_table("Right Flank Edit Costs")(input)?;

        Ok((
            input,
            Self {
                left_flank_length,
                right_flank_length,
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
                left_flank_edit_costs,
                right_flank_edit_costs,

                rq_qr_offset_costs,
                rr_qq_offset_costs,
                length_costs,
                length_difference_costs,
                forward_anti_primary_gap_costs,
                reverse_anti_primary_gap_costs,
            },
        ))
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
        writeln!(f, "# Limits")?;
        writeln!(f, "left_flank_length = {}", self.left_flank_length)?;
        writeln!(f, "right_flank_length = {}", self.right_flank_length)?;

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
        writeln!(f, "ForwardAntiPrimaryGap")?;
        writeln!(f, "{}", self.forward_anti_primary_gap_costs)?;
        writeln!(f, "ReverseAntiPrimaryGap")?;
        writeln!(f, "{}", self.reverse_anti_primary_gap_costs)?;

        writeln!(f, "{}", self.primary_edit_costs)?;
        writeln!(f, "{}", self.secondary_forward_edit_costs)?;
        writeln!(f, "{}", self.secondary_reverse_edit_costs)?;
        writeln!(f, "{}", self.left_flank_edit_costs)?;
        writeln!(f, "{}", self.right_flank_edit_costs)?;

        Ok(())
    }
}

impl<AlphabetType: Alphabet, Cost: AStarCost> FromStr for TemplateSwitchConfig<AlphabetType, Cost> {
    type Err = crate::error::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let result = Self::parse_plain(s)
            .map(|(_, result)| result)
            .map_err(translate_nom_error)?;
        result.verify()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::implementation::alphabets::dna_alphabet_or_n::DnaAlphabetOrN;
    use generic_a_star::cost::U64Cost;

    use crate::config::TemplateSwitchConfig;

    #[test]
    fn display() {
        let config = TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default();
        let string = config.to_string();
        assert_eq!(
            config,
            TemplateSwitchConfig::read_plain(string.as_bytes()).unwrap()
        );
    }
}
