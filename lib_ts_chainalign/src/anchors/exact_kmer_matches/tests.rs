use crate::anchors::{exact_kmer_matches::compute_exact_kmers, kmers::Kmer8};

#[test]
fn test_compute_exact_kmers() {
    let k = 4;
    let sequence = b"AAGACGTA";
    let output = compute_exact_kmers(sequence, k);

    let mut expected_output = vec![
        (Kmer8::from(b"AAGA".as_slice()), 0),
        (Kmer8::from(b"AGAC".as_slice()), 1),
        (Kmer8::from(b"GACG".as_slice()), 2),
        (Kmer8::from(b"ACGT".as_slice()), 3),
        (Kmer8::from(b"CGTA".as_slice()), 4),
    ];

    expected_output.sort_unstable();
    expected_output.dedup();
    assert_eq!(output, expected_output);
}
