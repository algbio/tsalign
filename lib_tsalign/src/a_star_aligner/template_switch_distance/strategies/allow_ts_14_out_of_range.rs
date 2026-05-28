use thiserror::Error;

use crate::a_star_aligner::alignment_geometry::{AlignmentCoordinates, AlignmentRange};

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Ts14OutOfRangeStrategy {
    /// TSMs are not allowed to start or end outside of the specified alignment ranges.
    Disallow,
    /// TSMs are allowed to start or end outside of the specified alignment ranges.
    Allow,
}

impl From<bool> for Ts14OutOfRangeStrategy {
    fn from(value: bool) -> Self {
        if value { Self::Allow } else { Self::Disallow }
    }
}

#[derive(Default)]
pub struct AdditionalExplicitTSMStartsAndEnds {
    pub explicit_tsm_starts: Vec<AlignmentCoordinates>,
    pub explicit_tsm_ends: Vec<AlignmentCoordinates>,
}

#[derive(Debug, Error)]
pub enum AdditionalExplicitTSMStartsAndEndsError {}

impl AdditionalExplicitTSMStartsAndEnds {
    pub fn new(
        original_reference: &str,
        original_query: &str,
        range: &AlignmentRange,
        skip_characters: Vec<char>,
        has_embedded_rq_ranges: bool,
    ) -> Result<Self, AdditionalExplicitTSMStartsAndEndsError> {
        let mut explicit_tsm_starts = Vec::new();
        let mut explicit_tsm_ends = Vec::new();

        let original_reference: Vec<_> = original_reference.chars().collect();
        let original_query: Vec<_> = original_query.chars().collect();

        let do_skip = |c| skip_characters.contains(&c);

        // Find range start in original sequences.
        let original_reference_start = if has_embedded_rq_ranges {
            original_reference
                .iter()
                .copied()
                .enumerate()
                .find(|(_, c)| *c == '|')
                .map(|(i, _)| i)
                .expect(
                    "Using embedded RQ ranges, but reference sequence contains no '|' character.",
                )
        } else {
            original_reference
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| do_skip(*c))
                .nth(range.reference_offset())
                .map(|(i, _)| i)
                .expect("reference alignment range must be inside of reference sequence")
        };
        let original_query_start = if has_embedded_rq_ranges {
            original_query
                .iter()
                .copied()
                .enumerate()
                .find(|(_, c)| *c == '|')
                .map(|(i, _)| i)
                .expect("Using embedded RQ ranges, but query sequence contains no '|' character.")
        } else {
            original_query
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| do_skip(*c))
                .nth(range.query_offset())
                .map(|(i, _)| i)
                .expect("query alignment range must be inside of query sequence")
        };

        // Move backwards from range start in original sequences to collect explicit TSM starts.
        let mut original_reference_index = original_reference_start;
        let mut original_query_index = original_query_start;
        let mut cleaned_reference_index = range.reference_offset();
        let mut cleaned_query_index = range.query_offset();

        // Loop invariant: the pair of original characters has already been added to the explicit TSM starts.
        //                 Since the start of the range is not an "additional" TSM start,
        //                 this invariant is true at the beginning with an empty list of explicit TSM starts.
        while original_reference_index > 0 && original_query_index > 0 {
            original_reference_index -= 1;
            original_query_index -= 1;
            let skip_reference = do_skip(original_reference[original_reference_index]);
            let skip_query = do_skip(original_query[original_query_index]);

            if skip_reference && skip_query {
                // Skip complete gap in both sequences.
            } else if skip_reference {
                // Gap-char alignment.
                cleaned_query_index -= 1;
                explicit_tsm_starts.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
            } else if skip_query {
                // Char-gap alignment.
                cleaned_reference_index -= 1;
                explicit_tsm_starts.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
            } else {
                // Char-char alignment.
                cleaned_reference_index -= 1;
                cleaned_query_index -= 1;
                explicit_tsm_starts.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
            }
        }

        // Find range end in original sequences.
        let original_reference_end = if has_embedded_rq_ranges {
            original_reference
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| *c == '|')
                .nth(1)
                .map(|(i, _)| i + 1)
                .expect("Using embedded RQ ranges, but reference sequence contains no second '|' character.")
        } else {
            original_reference
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| do_skip(*c))
                .nth(range.reference_limit())
                .map(|(i, _)| i)
                .unwrap_or(original_reference.len())
        };
        let original_query_end = if has_embedded_rq_ranges {
            original_query
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| *c == '|')
                .nth(1)
                .map(|(i, _)| i + 1)
                .expect("Using embedded RQ ranges, but query sequence contains no second '|' character.")
        } else {
            original_query
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, c)| do_skip(*c))
                .nth(range.query_limit())
                .map(|(i, _)| i)
                .unwrap_or(original_query.len())
        };

        // Move forwards from range end (exclusive) in original sequences to collect explicit TSM ends.
        let mut original_reference_index = original_reference_end;
        let mut original_query_index = original_query_end;
        let mut cleaned_reference_index = range.reference_limit();
        let mut cleaned_query_index = range.query_limit();

        // Loop invariant: the pair of original characters has already been added to the explicit TSM ends.
        //                 Since the end of the range is not an "additional" TSM end,
        //                 this invariant is true at the beginning with an empty list of explicit TSM ends.
        while original_reference_index < original_reference.len()
            && original_query_index < original_query.len()
        {
            let skip_reference = do_skip(original_reference[original_reference_index]);
            let skip_query = do_skip(original_query[original_query_index]);

            if skip_reference && skip_query {
                // Skip complete gap in both sequences.
            } else if skip_reference {
                // Gap-char alignment.
                explicit_tsm_ends.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
                cleaned_query_index += 1;
            } else if skip_query {
                // Char-gap alignment.
                explicit_tsm_ends.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
                cleaned_reference_index += 1;
            } else {
                // Char-char alignment.
                explicit_tsm_ends.push(AlignmentCoordinates::new(
                    cleaned_reference_index,
                    cleaned_query_index,
                ));
                cleaned_reference_index += 1;
                cleaned_query_index += 1;
            }

            original_reference_index += 1;
            original_query_index += 1;
        }

        explicit_tsm_starts.reverse();
        Ok(Self {
            explicit_tsm_starts,
            explicit_tsm_ends,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::a_star_aligner::{
        alignment_geometry::{AlignmentCoordinates, AlignmentRange},
        template_switch_distance::strategies::allow_ts_14_out_of_range::AdditionalExplicitTSMStartsAndEnds,
    };

    #[test]
    fn new() {
        let original_reference = "AC--GT|ACG-|-ACGT-";
        let original_query = "-ACGT-|A-GT|-AC-GT-";
        let range = AlignmentRange::new(4, 4, 7, 7);
        let skip_characters = vec!['-'];

        let additional_explicit_tsm_starts_and_ends = AdditionalExplicitTSMStartsAndEnds::new(
            original_reference,
            original_query,
            &range,
            skip_characters,
            true,
        )
        .unwrap();

        let expected_starts = vec![
            AlignmentCoordinates::new(0, 0),
            AlignmentCoordinates::new(1, 0),
            AlignmentCoordinates::new(2, 1),
            AlignmentCoordinates::new(2, 2),
            AlignmentCoordinates::new(2, 3),
            AlignmentCoordinates::new(3, 4),
        ];
        let expected_ends = vec![
            AlignmentCoordinates::new(7, 7),
            AlignmentCoordinates::new(8, 8),
            AlignmentCoordinates::new(9, 9),
            AlignmentCoordinates::new(10, 9),
            AlignmentCoordinates::new(11, 10),
        ];

        assert_eq!(
            additional_explicit_tsm_starts_and_ends.explicit_tsm_starts,
            expected_starts
        );
        assert_eq!(
            additional_explicit_tsm_starts_and_ends.explicit_tsm_ends,
            expected_ends
        );
    }
}
