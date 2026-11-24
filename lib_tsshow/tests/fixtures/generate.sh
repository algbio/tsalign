#!/usr/bin/env bash
# Regenerates the alignment fixtures of the lib_tsshow tests from their pair fasta inputs.
#
# The fixtures are committed, so that the snapshots of the visualisation do not change when the aligners change.
# Only rerun this deliberately (e.g. after a change of the toml format), then review and re-bless the snapshots:
#   lib_tsshow/tests/fixtures/generate.sh && TSSHOW_BLESS=1 cargo test -p lib_tsshow
#
# Inputs are the crafted sequences in `inputs/` and some sequences from `test_files/`.
# Configs: `test_files/config/small`, and variants of it in `configs/`:
# `forward` (cheap forward TSs), `forward_only` (cheap forward, expensive reverse TSs)
# and `inner` (expensive reverse TSs, so that the TS inners contain secondary edits).
set -euo pipefail

FIXTURES="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY="$(cd "$FIXTURES/../../.." && pwd)"
cd "$REPOSITORY"

cargo build --release -p tsalign
TSALIGN="${CARGO_TARGET_DIR:-$REPOSITORY/target}/release/tsalign"
CACHE="${CARGO_TARGET_DIR:-$REPOSITORY/target}/tsshow-fixture-cache"
mkdir -p "$CACHE"

# classic <name> <pair fasta> <config> <arguments for both runs> [arguments for the TS run...]
# Writes <name>.toml and <name>.no_ts.toml.
classic() {
    local name="$1" input="$2" config="$3"
    # Word splitting of the common arguments is intended.
    # shellcheck disable=SC2206
    local common=($4)
    shift 4
    "$TSALIGN" align -l warn -p "$input" -c "$config" -o "$FIXTURES/$name.toml" "${common[@]}" "$@"
    "$TSALIGN" align -l warn -p "$input" -c "$config" -o "$FIXTURES/$name.no_ts.toml" --no-ts "${common[@]}"
}

# chain <name> <pair fasta>
# Writes <name>.toml (tschainalign with its default config) and <name>.no_ts.toml (classic aligner without TSs).
chain() {
    local name="$1" input="$2"
    "$TSALIGN" align -l warn -p "$input" --alignment-method a-star-chain-ts --alphabet dna \
        -k 7 --max-anchor-mutations 1 --cache-directory "$CACHE" -o "$FIXTURES/$name.toml"
    "$TSALIGN" align -l warn -p "$input" -c test_files/config/chainalign --alphabet dna \
        -o "$FIXTURES/$name.no_ts.toml" --no-ts
}

SMALL=test_files/config/small
CONFIGS="$FIXTURES/configs"
INPUTS="$FIXTURES/inputs"

# Classic aligner (a-star-template-switch).
classic classic_qr_reverse test_files/twin_show_ts_small.fa "$SMALL" ""
classic classic_rq_reverse test_files/twin_show_ts_small2.fa "$SMALL" ""
classic classic_rr_qq_reverse test_files/twin_show_ts_small_rr_qq.fa "$SMALL" ""
classic classic_both_descendants test_files/twin_both_descendants.fa "$SMALL" ""
classic classic_negative_gap test_files/twin_in_place_inversion_1S.fa "$SMALL" ""
classic classic_negative_gap_2 test_files/twin_in_place_inversion_2S.fa "$SMALL" ""
classic classic_primary_indels test_files/twin_show_ts_small3.fa "$SMALL" ""
classic classic_three_ts test_files/twin_show_repeated_chars.fa "$SMALL" ""
classic classic_large_uncertainty test_files/twin_variable_length_ts.fa "$SMALL" ""
classic classic_no_uncertainty test_files/twin_show_ts_small_rr_qq.fa "$SMALL" "" \
    --ts-uncertainty-range-extension-mode none
classic classic_ranges test_files/twin_show_ts_optional.fa "$SMALL" "--rq-ranges R65..80Q65..75"
classic classic_adjacent_inversions "$INPUTS/adjacent_inversions.fa" "$SMALL" ""
classic classic_flank_snps "$INPUTS/inversion_with_flank_snps.fa" "$SMALL" ""
classic classic_inner_insertion "$INPUTS/inner_insertion.fa" "$CONFIGS/inner" ""
classic classic_inner_deletion "$INPUTS/inner_deletion.fa" "$CONFIGS/inner" "" --ts-total-length-strategy none
classic classic_inner_substitution "$INPUTS/inner_substitution.fa" "$CONFIGS/inner" ""
classic classic_forward_qr "$INPUTS/forward_duplication.fa" "$CONFIGS/forward" ""
classic classic_forward_qq "$INPUTS/forward_replacement.fa" "$CONFIGS/forward_only" ""

# tschainalign (a-star-chain-ts).
chain chain_adjacent_inversions "$INPUTS/adjacent_inversions.fa"
chain chain_inner_substitution "$INPUTS/inner_substitution.fa"
chain chain_flank_snps "$INPUTS/inversion_with_flank_snps.fa"
chain chain_rr_qq test_files/twin_show_ts_small_rr_qq.fa
chain chain_primary_indels test_files/twin_show_ts_small3.fa

# The runtime of the aligners differs between runs.
sed -i -E 's/^duration_seconds = .*/duration_seconds = 0.0/; s/^([a-z_]*_seconds) = [0-9.e+-]+$/\1 = 0.0/' "$FIXTURES"/*.toml
