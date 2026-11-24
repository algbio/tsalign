//! Unit tests of the TS arrangement.
//!
//! The arrangements are shown as text with one line per row, in the order of the SVG
//! (reference complement inners, reference complement, reference inners, reference, query, query inners,
//! query complement, query complement inners). Characters are shown as follows:
//! - source, complement and inner characters: the character of the sequence (lower case if substituted);
//! - hidden characters: `.`; gaps: `-`; separators: `|`; spacers: `_`; blanks: ` `;
//! - optional (uncertainty range) characters: `*` if they increase the inner, `~` if they decrease it.

use lib_tsalign::{
    aligner::{TemplateSwitchAncestor, TemplateSwitchDescendant, TemplateSwitchDirection},
    alignment_result::{AlignmentType, TSMUncertaintyRange},
};

use super::{
    TsArrangement, complement::ComplementChar, index_types::SourceColumn, inner::InnerChar,
    source::SourceChar,
};
use crate::svg::UncertaintyRangeMode;

const A: TemplateSwitchAncestor = TemplateSwitchAncestor::Reference;
const AQ: TemplateSwitchAncestor = TemplateSwitchAncestor::Query;
const D: TemplateSwitchDescendant = TemplateSwitchDescendant::Reference;
const DQ: TemplateSwitchDescendant = TemplateSwitchDescendant::Query;

fn complement(c: char) -> char {
    match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        other => other,
    }
}

fn source_char(sequence: &[char], column: SourceColumn, lower_case: bool) -> char {
    let c = sequence[column.primitive()];
    if lower_case {
        c.to_ascii_lowercase()
    } else {
        c
    }
}

fn show_source(sequence: &[char], c: &SourceChar) -> char {
    match c {
        SourceChar::Source {
            column, lower_case, ..
        } => source_char(sequence, *column, *lower_case),
        SourceChar::OptionalSource {
            increases_inner, ..
        } => {
            if *increases_inner {
                '*'
            } else {
                '~'
            }
        }
        SourceChar::Hidden { .. } => '.',
        SourceChar::Gap { .. } => '-',
        SourceChar::Separator => '|',
        SourceChar::Spacer => '_',
        SourceChar::Blank => ' ',
    }
}

fn show_complement(complement_sequence: &[char], c: &ComplementChar) -> char {
    match c {
        ComplementChar::Complement {
            column, lower_case, ..
        } => source_char(complement_sequence, *column, *lower_case),
        ComplementChar::Hidden { .. } => '.',
        ComplementChar::Gap { .. } => '-',
        ComplementChar::Blank => ' ',
    }
}

fn show_inner(descendant_sequence: &[char], c: &InnerChar) -> char {
    match c {
        InnerChar::Inner {
            column, lower_case, ..
        } => source_char(descendant_sequence, *column, *lower_case),
        InnerChar::OptionalInner {
            increases_inner, ..
        } => {
            if *increases_inner {
                '*'
            } else {
                '~'
            }
        }
        InnerChar::Gap { .. } => '-',
        InnerChar::Blank => ' ',
    }
}

/// Renders the arrangement as text, see the module documentation. Trailing blanks are removed.
fn show(arrangement: &TsArrangement, reference: &str, query: &str) -> String {
    let reference: Vec<char> = reference.chars().collect();
    let query: Vec<char> = query.chars().collect();
    let reference_c: Vec<char> = reference.iter().copied().map(complement).collect();
    let query_c: Vec<char> = query.iter().copied().map(complement).collect();
    let descendant = |descendant| match descendant {
        TemplateSwitchDescendant::Reference => &reference,
        TemplateSwitchDescendant::Query => &query,
    };

    let mut lines = Vec::new();
    let inner_line = |label: String, inner: &super::inner::TsInner| {
        let sequence = descendant(inner.template_switch().descendant);
        let line: String = inner
            .sequence()
            .iter_values()
            .map(|c| show_inner(sequence, c))
            .collect();
        format!("{label:>6}: {line}")
    };

    for (_, inner) in arrangement.reference_complement_inners().rev() {
        lines.push(inner_line(
            format!("I{}", inner.template_switch().index),
            inner,
        ));
    }
    lines.push(format!(
        "{:>6}: {}",
        "RC",
        arrangement
            .reference_complement()
            .iter_values()
            .map(|c| show_complement(&reference_c, c))
            .collect::<String>()
    ));
    for (_, inner) in arrangement.reference_inners().rev() {
        lines.push(inner_line(
            format!("I{}", inner.template_switch().index),
            inner,
        ));
    }
    lines.push(format!(
        "{:>6}: {}",
        "R",
        arrangement
            .reference()
            .iter_values()
            .map(|c| show_source(&reference, c))
            .collect::<String>()
    ));
    lines.push(format!(
        "{:>6}: {}",
        "Q",
        arrangement
            .query()
            .iter_values()
            .map(|c| show_source(&query, c))
            .collect::<String>()
    ));
    for (_, inner) in arrangement.query_inners() {
        lines.push(inner_line(
            format!("I{}", inner.template_switch().index),
            inner,
        ));
    }
    lines.push(format!(
        "{:>6}: {}",
        "QC",
        arrangement
            .query_complement()
            .iter_values()
            .map(|c| show_complement(&query_c, c))
            .collect::<String>()
    ));
    for (_, inner) in arrangement.query_complement_inners() {
        lines.push(inner_line(
            format!("I{}", inner.template_switch().index),
            inner,
        ));
    }

    lines
        .into_iter()
        .map(|line| line.trim_end().to_owned())
        .collect::<Vec<_>>()
        .join("\n")
}

fn entrance(
    first_offset: isize,
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
    direction: TemplateSwitchDirection,
) -> AlignmentType {
    AlignmentType::TemplateSwitchEntrance {
        first_offset,
        uncertainty_range: Some(TSMUncertaintyRange::new()),
        descendant,
        ancestor,
        direction,
    }
}

fn reverse(
    first_offset: isize,
    descendant: TemplateSwitchDescendant,
    ancestor: TemplateSwitchAncestor,
) -> AlignmentType {
    entrance(
        first_offset,
        descendant,
        ancestor,
        TemplateSwitchDirection::Reverse,
    )
}

fn exit(anti_descendant_gap: isize) -> AlignmentType {
    AlignmentType::TemplateSwitchExit {
        anti_descendant_gap,
    }
}

fn flat(alignment: &[(usize, AlignmentType)]) -> Vec<AlignmentType> {
    alignment
        .iter()
        .flat_map(|(multiplicity, alignment_type)| {
            std::iter::repeat_n(*alignment_type, *multiplicity)
        })
        .collect()
}

/// Creates the arrangement of a complete alignment of the two sequences, as `create_ts_svg` does by default.
fn arrange(
    reference: &str,
    query: &str,
    alignment: &[(usize, AlignmentType)],
    uncertainty_range_mode: UncertaintyRangeMode,
) -> TsArrangement {
    arrange_with_offset(reference, query, (0, 0), alignment, uncertainty_range_mode)
}

fn arrange_with_offset(
    reference: &str,
    query: &str,
    (reference_offset, query_offset): (usize, usize),
    alignment: &[(usize, AlignmentType)],
    uncertainty_range_mode: UncertaintyRangeMode,
) -> TsArrangement {
    let mut arrangement = TsArrangement::new(
        reference_offset,
        query_offset,
        reference.len(),
        query.len(),
        flat(alignment),
        uncertainty_range_mode,
    )
    .unwrap();
    arrangement.remove_empty_columns();
    arrangement
}

fn assert_arrangement(arrangement: &TsArrangement, reference: &str, query: &str, expected: &str) {
    let actual = show(arrangement, reference, query);
    let expected = expected.trim_matches('\n');
    assert_eq!(
        actual, expected,
        "\nactual:\n{actual}\n\nexpected:\n{expected}\n"
    );
}

#[test]
fn primary_alignment_only() {
    let reference = "ACGTACGT";
    let query = "ACTTACCGT";
    let arrangement = arrange(
        reference,
        query,
        &[
            (2, AlignmentType::PrimaryMatch),
            (1, AlignmentType::PrimarySubstitution),
            (3, AlignmentType::PrimaryMatch),
            (1, AlignmentType::PrimaryInsertion),
            (1, AlignmentType::PrimaryDeletion),
            (1, AlignmentType::PrimaryInsertion),
            (1, AlignmentType::PrimaryMatch),
        ],
        UncertaintyRangeMode::None,
    );
    assert_eq!(arrangement.template_switches().count(), 0);
    assert_arrangement(
        &arrangement,
        reference,
        query,
        r"
    RC: ...... . .
     R: ACgTAC-G-T
     Q: ACtTACC-GT
    QC: ....... ..",
    );
}

type Case = (
    &'static str,
    &'static str,
    &'static str,
    Vec<(usize, AlignmentType)>,
    UncertaintyRangeMode,
    &'static str,
);

/// The arrangements of all TS kinds and shapes.
///
/// Each case is the name, the reference, the query, the alignment, the uncertainty range mode and the expected arrangement.
/// The `offsets` case aligns the sequences without their first characters.
fn cases() -> Vec<Case> {
    use AlignmentType::*;
    use UncertaintyRangeMode as U;
    vec![
        // 2 -> 1: descendant query, ancestor reference (an inversion in the query).
        (
            "query_reference",
            "AAAACCGTTT",
            "AAAACGGTTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I0:     GGC
    RC: ....GGC ...
     R: AAAACCG_TTT
     Q: AAAA... TTT
    QC: ....... ...
",
        ),
        // 1 -> 2: descendant reference, ancestor query.
        (
            "reference_query",
            "AAAACGGTTT",
            "AAAACCGTTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, D, AQ)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    RC: ....... ...
     R: AAAA... TTT
     Q: AAAACCG_TTT
    QC: ....GGC ...
    I0:     GGC
",
        ),
        // 1 -> 1: descendant and ancestor reference (a hairpin within the reference).
        (
            "reference_reference",
            "AAAACGTTTGG",
            "AAAAGGTTTGG",
            vec![
                (4, PrimaryMatch),
                (1, reverse(5, D, A)),
                (2, SecondaryMatch),
                (1, exit(2)),
                (5, PrimaryMatch),
            ],
            U::None,
            r"
    I0:          GC
    RC: ......  .AA..
     R: AAAA..  TTTGG
     Q: AAAAGG__TTTGG
    QC: ......  .....
",
        ),
        // 2 -> 2: descendant and ancestor query.
        (
            "query_query",
            "AAAAGGTTTGG",
            "AAAACGTTTGG",
            vec![
                (4, PrimaryMatch),
                (1, reverse(5, DQ, AQ)),
                (2, SecondaryMatch),
                (1, exit(2)),
                (5, PrimaryMatch),
            ],
            U::None,
            r"
    RC: ......  .....
     R: AAAAGG__TTTGG
     Q: AAAA..  TTTGG
    QC: ......  .AA..
    I0:          GC
",
        ),
        // Forward TS: a duplication in the query.
        (
            "forward",
            "AACGTTGCA",
            "AACGTTCGTGCA",
            vec![
                (6, PrimaryMatch),
                (1, entrance(-4, DQ, A, TemplateSwitchDirection::Forward)),
                (3, SecondaryMatch),
                (1, exit(0)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    RC: ......    ...
    I0:   CGT
     R: AACGTT____GCA
     Q: AACGTT... GCA
    QC: ......... ...
",
        ),
        // Negative anti-descendant gap: characters of the anti-descendant are repeated.
        (
            "negative_gap",
            "AAAACCGTTT",
            "AAAACGGAACCGTTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(-2)),
                (8, PrimaryMatch),
            ],
            U::None,
            r"
    I0:           GGC
    RC: ....      GGC...
     R: AAAA____AACCGTTT
     Q: AAAA... AACCGTTT
    QC: ....... ........
",
        ),
        // Zero anti-descendant gap: the inner is inserted into the descendant.
        (
            "zero_gap",
            "AAAACCGTTT",
            "AAAACCGCGGTTT",
            vec![
                (7, PrimaryMatch),
                (1, reverse(0, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(0)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I0:     GGC
    RC: ....GGC    ...
     R: AAAACCG____TTT
     Q: AAAACCG... TTT
    QC: .......... ...
",
        ),
        // Larger positive anti-descendant gap than the inner: the anti-descendant skips characters.
        (
            "large_gap",
            "AAAACCGTTTT",
            "AAAACGGTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(5)),
                (2, PrimaryMatch),
            ],
            U::None,
            r"
    I0:     GGC
    RC: ....GGC....
     R: AAAACCGTTTT
     Q: AAAA...  TT
    QC: .......  ..
",
        ),
        // Inner with a secondary insertion, deletion and substitution.
        (
            "inner_edits",
            "AAAACCGTATTT",
            "AAAATAACGGTTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(6, DQ, A)),
                (1, SecondaryMatch),
                (1, SecondarySubstitution),
                (1, SecondaryInsertion),
                (1, SecondaryDeletion),
                (2, SecondaryMatch),
                (1, SecondaryInsertion),
                (1, SecondaryMatch),
                (1, exit(5)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I0:     TGGC-AaT
    RC: ....G-GCA-tA..
     R: AAAAC CGT ATTT
     Q: AAAA. ... .TT
    QC: ..... ... ...
",
        ),
        // Two adjacent TSs.
        (
            "adjacent",
            "AAAACCGGTATTT",
            "AAAACGGATACTTT",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (1, reverse(4, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (1, PrimaryMatch),
                (1, PrimaryInsertion),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I1:          AT A
    I0:     GGC
    RC: ....GGC .AT A ..
     R: AAAACCG_GTA_T-TT
     Q: AAAA... ... CTTT
    QC: ....... ... ....
",
        ),
        // Two TSs with overlapping ancestors.
        (
            "overlapping",
            "AAAACCGTTTAAAA",
            "AAAACGGTTTCGGAAAA",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
                (1, reverse(-3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(0)),
                (4, PrimaryMatch),
            ],
            U::None,
            r"
    I1:     GGC
    I0:     GGC
    RC: ....GGC ...    ....
     R: AAAACCG_TTT____AAAA
     Q: AAAA... TTT... AAAA
    QC: ....... ...... ....
",
        ),
        // TS that spans the complete sequences.
        (
            "complete",
            "AACCG",
            "CGGTT",
            vec![(1, reverse(5, DQ, A)), (5, SecondaryMatch), (1, exit(5))],
            U::None,
            r"
    I0: TTGGC
    RC: TTGGC
     R: AACCG
     Q: .....
    QC: .....
",
        ),
        // Alignment of a substring: separators.
        (
            "offsets",
            "GAAAACCGTTTG",
            "CAAAACGGTTTC",
            vec![
                (4, PrimaryMatch),
                (1, reverse(3, DQ, A)),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I0:       GGC
    RC: . ....GGC ... .
     R: G|AAAACCG_TTT|G
     Q: C|AAAA... TTT|C
    QC: . ....... ... .
",
        ),
        // Uncertainty range.
        (
            "uncertainty_none",
            "AAAACCGTTT",
            "AAAACGGTTT",
            vec![
                (4, PrimaryMatch),
                (
                    1,
                    AlignmentType::TemplateSwitchEntrance {
                        first_offset: 3,
                        uncertainty_range: Some(TSMUncertaintyRange {
                            start_left_shift: 1,
                            start_right_shift: 1,
                            end_left_shift: 1,
                            end_right_shift: 1,
                        }),
                        descendant: DQ,
                        ancestor: A,
                        direction: TemplateSwitchDirection::Reverse,
                    },
                ),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::None,
            r"
    I0:     GGC
    RC: ....GGC ...
     R: AAAACCG_TTT
     Q: AAAA... TTT
    QC: ....... ...
",
        ),
        (
            "uncertainty_inner_only",
            "AAAACCGTTT",
            "AAAACGGTTT",
            vec![
                (4, PrimaryMatch),
                (
                    1,
                    AlignmentType::TemplateSwitchEntrance {
                        first_offset: 3,
                        uncertainty_range: Some(TSMUncertaintyRange {
                            start_left_shift: 1,
                            start_right_shift: 1,
                            end_left_shift: 1,
                            end_right_shift: 1,
                        }),
                        descendant: DQ,
                        ancestor: A,
                        direction: TemplateSwitchDirection::Reverse,
                    },
                ),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::InnerOnly,
            r"
    I0:    *~G~*
    RC: ....GGC ...
     R: AAAACCG_TTT
     Q: AAAA... TTT
    QC: ....... ...
",
        ),
        (
            "uncertainty_full",
            "AAAACCGTTT",
            "AAAACGGTTT",
            vec![
                (4, PrimaryMatch),
                (
                    1,
                    AlignmentType::TemplateSwitchEntrance {
                        first_offset: 3,
                        uncertainty_range: Some(TSMUncertaintyRange {
                            start_left_shift: 1,
                            start_right_shift: 1,
                            end_left_shift: 1,
                            end_right_shift: 1,
                        }),
                        descendant: DQ,
                        ancestor: A,
                        direction: TemplateSwitchDirection::Reverse,
                    },
                ),
                (3, SecondaryMatch),
                (1, exit(3)),
                (3, PrimaryMatch),
            ],
            U::Full,
            r"
    I0:    *~G~*
    RC: ....GGC ...
     R: AAAACCG_TTT
     Q: AAA*~.~ *TT
    QC: ....... ...
",
        ),
    ]
}

#[test]
fn arrangements() {
    let mut failures = Vec::new();
    for (name, reference, query, alignment, mode, expected) in cases() {
        let offsets = if name == "offsets" { (1, 1) } else { (0, 0) };
        let arrangement = arrange_with_offset(reference, query, offsets, &alignment, mode);
        let actual = show(&arrangement, reference, query);
        if actual != expected.trim_matches('\n') {
            failures.push(format!(
                "=== {name}\nactual:\n{actual}\nexpected:\n{}",
                expected.trim_matches('\n')
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n\n"));
}
