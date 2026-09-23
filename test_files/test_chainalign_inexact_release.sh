#!/bin/bash

cargo run --release -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign --alphabet dna --max-anchor-mutations 0 -p "$@"
cargo run --release -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign --alphabet dna --max-anchor-mutations 1 -p "$@"
cargo run --release -- align --alignment-method a-star-chain-ts -c test_files/config/chainalign --alphabet dna --max-anchor-mutations 2 -p "$@"