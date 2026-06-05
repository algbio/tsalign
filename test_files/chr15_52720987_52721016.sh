#!/usr/bin/env bash

cargo run -- align -p test_files/chr15_52720987_52721016.fa -o test_files/chr15_52720987_52721016.toml --alignment-method a-star-template-switch --skip-characters N- --alphabet dna -c test_files/config/bench --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy preprocess-price --allow-ts-14-out-of-range -k 0 --max-chaining-successors 0 --max-exact-cost-function-cost 0 --chaining-closed-list special --chaining-open-list linear-heap --rq-ranges R490..555Q490..542
