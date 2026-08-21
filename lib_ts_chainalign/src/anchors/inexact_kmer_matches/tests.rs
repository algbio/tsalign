use generic_a_star::cost::U16Cost;

use crate::{
    anchors::{
        inexact_kmer_matches::{
            generate_kmer_deletions, generate_kmer_insertions, generate_kmer_substitutions,
        },
        kmers::Kmer8,
    },
    costs::GapAffineCosts,
};

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
fn test_generate_inexact_kmers() {
    todo!()
}
