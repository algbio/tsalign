#![deny(clippy::mod_module_files)] // Not `forbid`: the clap derive allows `clippy::restriction`.

pub mod aligner;
pub mod alignment_geometry;
pub mod alignment_result;
pub mod arena_index;
pub mod config;
pub mod costs;
pub mod error;
mod io;
#[cfg(test)]
mod tests;

pub use aligner::configurable::Aligner;
