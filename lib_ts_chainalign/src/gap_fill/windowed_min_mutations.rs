//! Fill the gap between two chained inexact anchors with an optimal alignment.

/// The gap alignments of the windowed min mutations lower bound use node indices of type `u32`, and they do not
/// return errors. Their searches are limited to a single gap, so they need far fewer nodes than the arena can hold.
/// This lower bound is disabled in the CLI and the Python bindings (#224).
const ARENA_FULL: &str = "a windowed gap alignment does not fill the node arena";

pub mod gap_affine;
pub mod ts_12_jump;
pub mod ts_34_jump;
