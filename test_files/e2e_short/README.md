# Short inputs for the end-to-end tests

The end-to-end tests of tschainalign (`tsalign-tests/tests/chain_test_files.rs`) run every input with several values of
`k` and `max_anchor_mutations`. On some of the original files in `test_files`, single runs take minutes (issue #223).
Slow files that do not test a specific edge case are tested with fewer parameter combinations instead.
The files here are shorter copies of the files where that does not help:
they are built the same way (same blocks in the same order, only shorter homopolymer runs).
The original files are unchanged, since other scripts in `test_files` use them.

The classic A* aligner with template switches (`-c test_files/config/chainalign`) finds the same template switches with
the same cost on the original and on the short file:

| File | Why it is shortened | Original | Short | Optimal alignment, cost |
|---|---|---|---|---|
| `twin_both_descendants.fa` | Tests the option to force all TSM descendants to be the same; the classic aligner with template switches alone takes 2.5 minutes | blocks of 20 | blocks of 10 | RRR and QRR template switch over one block, cost 4 |
| `twin_equal_cost_ignore_geometry.fa` | Every run, also with exact anchors, takes more than a minute | A64 GGGG/CCCC T64 | A15 GGGG/CCCC T15 | one QRR template switch over the whole sequence, cost 2 |
