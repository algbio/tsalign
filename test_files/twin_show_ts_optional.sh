#!/usr/bin/env bash
set -euo pipefail

cargo run -- align -p test_files/twin_show_ts_optional.fa -c test_files/config/small --alphabet dna-n -o test_files/twin_show_ts_optional_no_ts.toml --rq-ranges R65..80Q65..75 --no-ts
cargo run -- align -p test_files/twin_show_ts_optional.fa -c test_files/config/small --alphabet dna-n -o test_files/twin_show_ts_optional.toml --rq-ranges R65..80Q65..75
cargo run -- show -i test_files/twin_show_ts_optional.toml -n test_files/twin_show_ts_optional_no_ts.toml -pas test_files/twin_show_ts_optional.svg -e full