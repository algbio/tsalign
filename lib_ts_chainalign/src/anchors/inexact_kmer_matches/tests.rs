use generic_a_star::cost::U16Cost;
use num_traits::{Bounded, Zero};

use crate::{
    anchors::{
        inexact_kmer_matches::{
            compute_inexact_kmers, generate_kmer_deletions, generate_kmer_insertions,
            generate_kmer_substitutions,
        },
        kmers::Kmer8,
    },
    costs::GapAffineCosts,
};

fn generate_insertions(kmer: &[u8]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    for i in 0..=kmer.len() {
        for c in b"ACGT" {
            let mut new_kmer = kmer.to_vec();
            new_kmer.insert(i, *c);
            result.push(new_kmer);
        }
    }
    result
}

fn generate_insertion_tuples(
    kmer: &[u8],
    offset: usize,
    cost: u16,
) -> Vec<(Kmer8, usize, U16Cost)> {
    generate_insertions(kmer)
        .into_iter()
        .map(|kmer| (Kmer8::from(kmer.as_slice()), offset, U16Cost::from(cost)))
        .collect()
}

fn generate_double_insertion_tuples(
    kmer: &[u8],
    offset: usize,
    costs: &GapAffineCosts<U16Cost>,
) -> Vec<(Kmer8, usize, U16Cost)> {
    let mut result = Vec::new();
    for i in 0..=kmer.len() {
        for j in i..=kmer.len() {
            for c1 in b"ACGT" {
                for c2 in b"ACGT" {
                    let mut new_kmer = kmer.to_vec();
                    new_kmer.insert(i, *c1);
                    new_kmer.insert(j + 1, *c2);
                    let cost = if i == j {
                        costs.gap_open + costs.gap_extend
                    } else {
                        costs.gap_open + costs.gap_open
                    };
                    result.push((Kmer8::from(new_kmer.as_slice()), offset, cost));
                }
            }
        }
    }
    result
}

#[test]
fn test_generate_kmer_insertions_0() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"AAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_insertions(kmer, k, 0, &costs, &mut output);

    let mut expected_output = vec![(Kmer8::from(b"AAGA".as_slice()), U16Cost::from(0u16))];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_insertions_1() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"AAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_insertions(kmer, k, 1, &costs, &mut output);

    let mut expected_output = vec![
        (Kmer8::from(b"AAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"CAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"GAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"TAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"ACAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AGAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"ATAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAAGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AACGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AATGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGAA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGCA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGTA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGAA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGAC".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGAG".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAGAT".as_slice()), U16Cost::from(3u16)),
    ];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_insertions_2() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 2;
    let kmer = Kmer8::from(b"GA".as_slice());
    let mut output = Vec::new();
    generate_kmer_insertions(kmer, k, 2, &costs, &mut output);
    let kmer = kmer.to_vec(k);

    let double_insertions = [
        b"AA", b"AC", b"AG", b"AT", b"CA", b"CC", b"CG", b"CT", b"GA", b"GC", b"GG", b"GT", b"TA",
        b"TC", b"TG", b"TT",
    ];

    let mut new_kmer = Vec::new();
    let mut expected_output = Vec::new();

    // Generate insertions.
    for i in 0..=2 {
        for j in i..=2 {
            for insertion in &double_insertions {
                new_kmer.extend_from_slice(&kmer[..i]);
                new_kmer.extend_from_slice(&insertion[0..=0]);
                new_kmer.extend_from_slice(&kmer[i..j]);
                new_kmer.extend_from_slice(&insertion[1..=1]);
                new_kmer.extend_from_slice(&kmer[j..]);
                expected_output.push((
                    Kmer8::from(new_kmer.as_slice()),
                    if i == j {
                        U16Cost::from(4u16)
                    } else {
                        U16Cost::from(6u16)
                    },
                ));
                new_kmer.clear();
            }
        }
    }

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_deletions_0() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"AAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_deletions(kmer, k, 0, &costs, &mut output);

    let mut expected_output = vec![(Kmer8::from(b"AAGA".as_slice()), U16Cost::from(0u16))];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_deletions_1() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"AAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_deletions(kmer, k, 1, &costs, &mut output);

    let mut expected_output = vec![
        (Kmer8::from(b"AGA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAA".as_slice()), U16Cost::from(3u16)),
        (Kmer8::from(b"AAG".as_slice()), U16Cost::from(3u16)),
    ];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_deletions_2() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"CAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_deletions(kmer, k, 2, &costs, &mut output);

    let mut expected_output = vec![
        (Kmer8::from(b"GA".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"CA".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"AA".as_slice()), U16Cost::from(6u16)),
        (Kmer8::from(b"AG".as_slice()), U16Cost::from(6u16)),
        (Kmer8::from(b"CG".as_slice()), U16Cost::from(6u16)),
    ];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_deletions_empty() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"CAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_deletions(kmer, k, 4, &costs, &mut output);

    let mut expected_output = vec![(Kmer8::from(b"".as_slice()), U16Cost::from(6u16))];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_substitutions_0() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let kmer = Kmer8::from(b"AAGA".as_slice());
    let mut output = Vec::new();
    generate_kmer_substitutions(kmer, k, 0, &costs, &mut output);

    let mut expected_output = vec![(Kmer8::from(b"AAGA".as_slice()), U16Cost::from(0u16))];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_substitutions_1() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 2;
    let kmer = Kmer8::from(b"GA".as_slice());
    let mut output = Vec::new();
    generate_kmer_substitutions(kmer, k, 1, &costs, &mut output);

    let mut expected_output = vec![
        (Kmer8::from(b"AA".as_slice()), U16Cost::from(2u16)),
        (Kmer8::from(b"CA".as_slice()), U16Cost::from(2u16)),
        (Kmer8::from(b"TA".as_slice()), U16Cost::from(2u16)),
        (Kmer8::from(b"GC".as_slice()), U16Cost::from(2u16)),
        (Kmer8::from(b"GG".as_slice()), U16Cost::from(2u16)),
        (Kmer8::from(b"GT".as_slice()), U16Cost::from(2u16)),
    ];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_generate_kmer_substitutions_2() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 2;
    let kmer = Kmer8::from(b"GA".as_slice());
    let mut output = Vec::new();
    generate_kmer_substitutions(kmer, k, 2, &costs, &mut output);

    let mut expected_output = vec![
        (Kmer8::from(b"AC".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"AG".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"AT".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"CC".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"CG".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"CT".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"TC".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"TG".as_slice()), U16Cost::from(4u16)),
        (Kmer8::from(b"TT".as_slice()), U16Cost::from(4u16)),
    ];

    output.sort_unstable();
    output.dedup();
    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}

#[test]
fn test_compute_inexact_kmers_0() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let sequence = b"AAGACGTA";
    let output = compute_inexact_kmers::<u16, _>(sequence, k, 0, &costs);

    let mut expected_output = vec![vec![
        (Kmer8::from(b"AAGA".as_slice()), 0, U16Cost::zero()),
        (Kmer8::from(b"AGAC".as_slice()), 1, U16Cost::zero()),
        (Kmer8::from(b"GACG".as_slice()), 2, U16Cost::zero()),
        (Kmer8::from(b"ACGT".as_slice()), 3, U16Cost::zero()),
        (Kmer8::from(b"CGTA".as_slice()), 4, U16Cost::zero()),
    ]];

    for vec in &mut expected_output {
        vec.sort_unstable();
        vec.dedup();
    }

    assert_eq!(output, expected_output);
}

#[test]
fn test_compute_inexact_kmers_1() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 4;
    let sequence = b"AAGACGTA";
    let output = compute_inexact_kmers::<u16, _>(sequence, k, 1, &costs);

    fn generate_insertions(kmer: &[u8]) -> Vec<Vec<u8>> {
        let mut result = Vec::new();
        for i in 0..=kmer.len() {
            for c in b"ACGT" {
                let mut new_kmer = kmer.to_vec();
                new_kmer.insert(i, *c);
                result.push(new_kmer);
            }
        }
        result
    }

    let mut expected_output = vec![
        vec![
            // AAGA
            (Kmer8::from(b"AGA".as_slice()), 0, U16Cost::from(3u16)),
            (Kmer8::from(b"AAA".as_slice()), 0, U16Cost::from(3u16)),
            (Kmer8::from(b"AAG".as_slice()), 0, U16Cost::from(3u16)),
            // AGAC
            (Kmer8::from(b"GAC".as_slice()), 1, U16Cost::from(3u16)),
            (Kmer8::from(b"AAC".as_slice()), 1, U16Cost::from(3u16)),
            (Kmer8::from(b"AGC".as_slice()), 1, U16Cost::from(3u16)),
            (Kmer8::from(b"AGA".as_slice()), 1, U16Cost::from(3u16)),
            // GACG
            (Kmer8::from(b"ACG".as_slice()), 2, U16Cost::from(3u16)),
            (Kmer8::from(b"GCG".as_slice()), 2, U16Cost::from(3u16)),
            (Kmer8::from(b"GAG".as_slice()), 2, U16Cost::from(3u16)),
            (Kmer8::from(b"GAC".as_slice()), 2, U16Cost::from(3u16)),
            // ACGT
            (Kmer8::from(b"CGT".as_slice()), 3, U16Cost::from(3u16)),
            (Kmer8::from(b"AGT".as_slice()), 3, U16Cost::from(3u16)),
            (Kmer8::from(b"ACT".as_slice()), 3, U16Cost::from(3u16)),
            (Kmer8::from(b"ACG".as_slice()), 3, U16Cost::from(3u16)),
            // CGTA
            (Kmer8::from(b"GTA".as_slice()), 4, U16Cost::from(3u16)),
            (Kmer8::from(b"CTA".as_slice()), 4, U16Cost::from(3u16)),
            (Kmer8::from(b"CGA".as_slice()), 4, U16Cost::from(3u16)),
            (Kmer8::from(b"CGT".as_slice()), 4, U16Cost::from(3u16)),
        ],
        vec![
            (Kmer8::from(b"AAGA".as_slice()), 0, U16Cost::zero()),
            (Kmer8::from(b"CAGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"GAGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"TAGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ACGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AGGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ATGA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AAAA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AACA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AATA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AAGC".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AAGG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AAGT".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AGAC".as_slice()), 1, U16Cost::zero()),
            (Kmer8::from(b"CGAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"GGAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"TGAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AAAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"ACAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"ATAC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGCC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGGC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGTC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGAA".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGAG".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AGAT".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"GACG".as_slice()), 2, U16Cost::zero()),
            (Kmer8::from(b"AACG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"CACG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"TACG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GCCG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GGCG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GTCG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GAAG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GAGG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GATG".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GACA".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GACC".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"GACT".as_slice()), 2, U16Cost::from(2u16)),
            (Kmer8::from(b"ACGT".as_slice()), 3, U16Cost::zero()),
            (Kmer8::from(b"CCGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"GCGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"TCGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"AAGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"AGGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ATGT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACAT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACCT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACTT".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACGA".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACGC".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"ACGG".as_slice()), 3, U16Cost::from(2u16)),
            (Kmer8::from(b"CGTA".as_slice()), 4, U16Cost::zero()),
            (Kmer8::from(b"AGTA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"GGTA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"TGTA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CATA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CCTA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CTTA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGAA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGCA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGGA".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGTC".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGTG".as_slice()), 4, U16Cost::from(2u16)),
            (Kmer8::from(b"CGTT".as_slice()), 4, U16Cost::from(2u16)),
        ],
        generate_insertions(b"AAGA")
            .into_iter()
            .map(|kmer| (Kmer8::from(kmer.as_slice()), 0, U16Cost::from(3u16)))
            .chain(
                generate_insertions(b"AGAC")
                    .into_iter()
                    .map(|kmer| (Kmer8::from(kmer.as_slice()), 1, U16Cost::from(3u16))),
            )
            .chain(
                generate_insertions(b"GACG")
                    .into_iter()
                    .map(|kmer| (Kmer8::from(kmer.as_slice()), 2, U16Cost::from(3u16))),
            )
            .chain(
                generate_insertions(b"ACGT")
                    .into_iter()
                    .map(|kmer| (Kmer8::from(kmer.as_slice()), 3, U16Cost::from(3u16))),
            )
            .chain(
                generate_insertions(b"CGTA")
                    .into_iter()
                    .map(|kmer| (Kmer8::from(kmer.as_slice()), 4, U16Cost::from(3u16))),
            )
            .collect(),
    ];

    for vec in &mut expected_output {
        vec.sort_unstable();
        vec.dedup();
    }

    assert_eq!(output, expected_output, "{}", {
        use std::fmt::Write;
        let mut result = String::new();
        writeln!(result, "Output:").unwrap();
        for (i, vec) in output.iter().enumerate() {
            let k = k + i - 1;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }

        writeln!(result, "Expected output:").unwrap();
        for (i, vec) in expected_output.iter().enumerate() {
            let k = k + i - 1;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }
        result
    });
}

#[test]
fn test_compute_inexact_kmers_2_2() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 2;
    let sequence = b"AAG";
    let output = compute_inexact_kmers::<u16, _>(sequence, k, 2, &costs);

    let mut expected_output = vec![
        vec![
            // Del Del
            (Kmer8::from(b"".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"".as_slice()), 1, U16Cost::from(4u16)),
        ],
        vec![
            // Del
            (Kmer8::from(b"A".as_slice()), 0, U16Cost::from(3u16)),
            (Kmer8::from(b"A".as_slice()), 1, U16Cost::from(3u16)),
            (Kmer8::from(b"G".as_slice()), 1, U16Cost::from(3u16)),
            // Del Sub
            (Kmer8::from(b"C".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"G".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"T".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"C".as_slice()), 1, U16Cost::from(5u16)),
            (Kmer8::from(b"T".as_slice()), 1, U16Cost::from(5u16)),
        ],
        vec![
            // Nothing
            (Kmer8::from(b"AA".as_slice()), 0, U16Cost::from(0u16)),
            (Kmer8::from(b"AG".as_slice()), 1, U16Cost::from(0u16)),
            // Sub
            (Kmer8::from(b"CA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"GA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"TA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AC".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AT".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"CG".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"GG".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"TG".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AA".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AC".as_slice()), 1, U16Cost::from(2u16)),
            (Kmer8::from(b"AT".as_slice()), 1, U16Cost::from(2u16)),
            // Sub Sub
            (Kmer8::from(b"CC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CA".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"GA".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"TA".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"CC".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"GC".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"TC".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"CT".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"GT".as_slice()), 1, U16Cost::from(4u16)),
            (Kmer8::from(b"TT".as_slice()), 1, U16Cost::from(4u16)),
            // Del Ins are always more expensive than Sub Sub in this case.
        ],
        // Ins
        generate_insertion_tuples(b"AA", 0, 3)
            .into_iter()
            .chain(generate_insertion_tuples(b"AG", 1, 3))
            // Sub Ins
            .chain(generate_insertion_tuples(b"CA", 0, 5))
            .chain(generate_insertion_tuples(b"GA", 0, 5))
            .chain(generate_insertion_tuples(b"TA", 0, 5))
            .chain(generate_insertion_tuples(b"AC", 0, 5))
            .chain(generate_insertion_tuples(b"AG", 0, 5))
            .chain(generate_insertion_tuples(b"AT", 0, 5))
            .chain(generate_insertion_tuples(b"CG", 1, 5))
            .chain(generate_insertion_tuples(b"GG", 1, 5))
            .chain(generate_insertion_tuples(b"TG", 1, 5))
            .chain(generate_insertion_tuples(b"AA", 1, 5))
            .chain(generate_insertion_tuples(b"AC", 1, 5))
            .chain(generate_insertion_tuples(b"AT", 1, 5))
            .collect(),
        // Ins Ins
        generate_double_insertion_tuples(b"AA", 0, &costs)
            .into_iter()
            .chain(generate_double_insertion_tuples(b"AG", 1, &costs))
            .collect(),
    ];

    for vec in &mut expected_output {
        vec.sort_unstable();
        let mut previous_kmer = Kmer8::default();
        let mut previous_offset = usize::MAX;
        let mut previous_cost = U16Cost::max_value();
        vec.retain(|(kmer, offset, cost)| {
            if *kmer == previous_kmer && *offset == previous_offset {
                debug_assert!(previous_cost <= *cost);
                false
            } else {
                previous_kmer = *kmer;
                previous_offset = *offset;
                previous_cost = *cost;
                true
            }
        });
    }

    assert_eq!(output, expected_output, "{}", {
        use std::fmt::Write;
        let mut result = String::new();
        writeln!(result, "Output:").unwrap();
        for (i, vec) in output.iter().enumerate() {
            let k = k + i - 2;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }

        writeln!(result, "Expected output:").unwrap();
        for (i, vec) in expected_output.iter().enumerate() {
            let k = k + i - 2;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }
        result
    });
}

#[test]
fn test_compute_inexact_kmers_2_3() {
    let costs = GapAffineCosts::new(
        U16Cost::from(2u16),
        U16Cost::from(3u16),
        U16Cost::from(1u16),
    );
    let k = 3;
    let sequence = b"ATG";
    let output = compute_inexact_kmers::<u16, _>(sequence, k, 2, &costs);

    let mut expected_output = vec![
        vec![
            // Del Del
            (Kmer8::from(b"A".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"G".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"T".as_slice()), 0, U16Cost::from(6u16)),
        ],
        vec![
            // Del
            (Kmer8::from(b"AT".as_slice()), 0, U16Cost::from(3u16)),
            (Kmer8::from(b"AG".as_slice()), 0, U16Cost::from(3u16)),
            (Kmer8::from(b"TG".as_slice()), 0, U16Cost::from(3u16)),
            // Del Sub
            (Kmer8::from(b"CT".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"GT".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"TT".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AA".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AC".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"CG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"GG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"TG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AA".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AC".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AT".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"AG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"CG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"GG".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"TA".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"TC".as_slice()), 0, U16Cost::from(5u16)),
            (Kmer8::from(b"TT".as_slice()), 0, U16Cost::from(5u16)),
        ],
        vec![
            // Nothing
            (Kmer8::from(b"ATG".as_slice()), 0, U16Cost::from(0u16)),
            // Sub
            (Kmer8::from(b"CTG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"GTG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"TTG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AAG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ACG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"AGG".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ATA".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ATC".as_slice()), 0, U16Cost::from(2u16)),
            (Kmer8::from(b"ATT".as_slice()), 0, U16Cost::from(2u16)),
            // Sub Sub
            (Kmer8::from(b"CAG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CCG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CGG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GAG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GCG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GGG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TAG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TCG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TGG".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CTA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CTC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"CTT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GTA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GTC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"GTT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TTA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TTC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"TTT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AAA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AAC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AAT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"ACA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"ACC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"ACT".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AGA".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AGC".as_slice()), 0, U16Cost::from(4u16)),
            (Kmer8::from(b"AGT".as_slice()), 0, U16Cost::from(4u16)),
        ]
        .into_iter()
        // Del Ins
        .chain(generate_insertion_tuples(b"TG", 0, 6))
        .chain(generate_insertion_tuples(b"AG", 0, 6))
        .chain(generate_insertion_tuples(b"AT", 0, 6))
        .collect(),
        // Ins
        generate_insertion_tuples(b"ATG", 0, 3)
            .into_iter()
            // Sub Ins
            .chain(generate_insertion_tuples(b"CTG", 0, 5))
            .chain(generate_insertion_tuples(b"GTG", 0, 5))
            .chain(generate_insertion_tuples(b"TTG", 0, 5))
            .chain(generate_insertion_tuples(b"AAG", 0, 5))
            .chain(generate_insertion_tuples(b"ACG", 0, 5))
            .chain(generate_insertion_tuples(b"AGG", 0, 5))
            .chain(generate_insertion_tuples(b"ATA", 0, 5))
            .chain(generate_insertion_tuples(b"ATC", 0, 5))
            .chain(generate_insertion_tuples(b"ATT", 0, 5))
            .collect(),
        // Ins Ins
        generate_double_insertion_tuples(b"ATG", 0, &costs),
    ];

    for vec in &mut expected_output {
        vec.sort_unstable();
        let mut previous_kmer = Kmer8::default();
        let mut previous_offset = usize::MAX;
        let mut previous_cost = U16Cost::max_value();
        vec.retain(|(kmer, offset, cost)| {
            if *kmer == previous_kmer && *offset == previous_offset {
                debug_assert!(previous_cost <= *cost);
                false
            } else {
                previous_kmer = *kmer;
                previous_offset = *offset;
                previous_cost = *cost;
                true
            }
        });
    }

    assert_eq!(output, expected_output, "{}", {
        use std::fmt::Write;
        let mut result = String::new();
        writeln!(result, "Output:").unwrap();
        for (i, vec) in output.iter().enumerate() {
            let k = k + i - 2;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }

        writeln!(result, "Expected output:").unwrap();
        for (i, vec) in expected_output.iter().enumerate() {
            let k = k + i - 2;
            for kmer in vec {
                writeln!(result, "{}, {}, {}", kmer.0.to_string(k), kmer.1, kmer.2).unwrap()
            }
        }
        result
    });
}
