#!/usr/bin/env bash

RUST_BACKTRACE=1 /usr/bin/time -v cargo run -- align -p test_files/ari-41.fa --alignment-method a-star-template-switch --skip-characters - --alphabet dna-n -c test_files/config/experiments --ts-node-ord-strategy anti-diagonal --ts-min-length-strategy lookahead --cost-limit 7 --rq-ranges R490..514Q490..521 --log-level trace