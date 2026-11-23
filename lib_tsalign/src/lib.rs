#![forbid(clippy::mod_module_files)]

pub mod aligner;
pub mod alignment_geometry;
pub mod alignment_result;
pub mod config;
pub mod costs;
pub mod error;
mod io;
#[cfg(test)]
mod tests;

pub use aligner::configurable::Aligner;
