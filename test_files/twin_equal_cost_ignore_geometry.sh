#!/usr/bin/env bash

cargo run --release -- align -p test_files/twin_equal_cost_ignore_geometry.fa -o test_files/twin_equal_cost_ignore_geometry.equal_cost.toml -c test_files/config/small --ts-uncertainty-range-extension-mode equal-cost
cargo run --release -- align -p test_files/twin_equal_cost_ignore_geometry.fa -o test_files/twin_equal_cost_ignore_geometry.equal_cost_ignore_geometry.toml -c test_files/config/small --ts-uncertainty-range-extension-mode equal-cost-ignore-geometry
cargo run -- show -i test_files/twin_equal_cost_ignore_geometry.equal_cost.toml -pas test_files/twin_equal_cost_ignore_geometry.equal_cost.inner-only.svg --uncertainty-range-mode inner-only
cargo run -- show -i test_files/twin_equal_cost_ignore_geometry.equal_cost_ignore_geometry.toml -pas test_files/twin_equal_cost_ignore_geometry.equal_cost_ignore_geometry.inner-only.svg --uncertainty-range-mode inner-only
cargo run -- show -i test_files/twin_equal_cost_ignore_geometry.equal_cost.toml -pas test_files/twin_equal_cost_ignore_geometry.equal_cost.full.svg --uncertainty-range-mode full
cargo run -- show -i test_files/twin_equal_cost_ignore_geometry.equal_cost_ignore_geometry.toml -pas test_files/twin_equal_cost_ignore_geometry.equal_cost_ignore_geometry.full.svg --uncertainty-range-mode full