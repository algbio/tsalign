#!/usr/bin/env bash

cargo build

target/debug/tsalign align -p test_files/chr11_39595337_39595347.fa -o test_files/chr11_39595337_39595347_no_ts.toml --alignment-method a-star-template-switch --skip-characters 'N-' --alphabet dna -c 'test_files/config/bench' --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy lookahead --allow-ts-14-out-of-range -k 0 --max-chaining-successors 0 --max-exact-cost-function-cost 0 --chaining-closed-list special --chaining-open-list linear-heap --no-ts --rq-ranges R490..527Q490..527
