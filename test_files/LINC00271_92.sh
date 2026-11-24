#!/usr/bin/env bash

cargo build

if [ ! -f test_files/LINC00271_92.toml ]; then
    cargo build --release
    target/release/tsalign align -p test_files/LINC00271_92.fa -o test_files/LINC00271_92.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy preprocess-price --rq-ranges R196..227Q196..202
    target/release/tsalign align -p test_files/LINC00271_92.fa -o test_files/LINC00271_92.no_ts.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --rq-ranges R196..227Q196..202 --no-ts
fi

target/debug/tsalign show -i test_files/LINC00271_92.toml -n test_files/LINC00271_92.no_ts.toml -ps test_files/LINC00271_92.svg --log-level trace
