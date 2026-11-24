#!/usr/bin/env bash

cargo build

if [ ! -f test_files/STXBP5-AS1_99.toml ]; then
    cargo build --release
    target/release/tsalign align -p test_files/STXBP5-AS1_99.fa -o test_files/STXBP5-AS1_99.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy preprocess-price --rq-ranges R196..237Q196..204
    target/release/tsalign align -p test_files/STXBP5-AS1_99.fa -o test_files/STXBP5-AS1_99.no_ts.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --rq-ranges R196..237Q196..204 --no-ts
fi

target/debug/tsalign show -i test_files/STXBP5-AS1_99.toml -n test_files/STXBP5-AS1_99.no_ts.toml -ps test_files/STXBP5-AS1_99.svg --log-level trace
