//! Bindings for tschainalign (`lib_ts_chainalign`).

use std::{
    collections::{BTreeSet, HashMap},
    fs::File,
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use compact_genome::implementation::alphabets::{
    dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
};
use generic_a_star::cost::U32Cost;
use lib_ts_chainalign::{
    ChainAlignError,
    cache::{cache_file_path, infer_k, infer_max_n, write_cache_file_atomically},
    chain_align::{
        limits::AlignmentLimits,
        performance_parameters::{
            AlignmentPerformanceParameters, AnchorPruning, ChainingClosedList, ChainingOpenList,
            InexactPruningParameters, InexactPruningRules,
        },
    },
    chaining_lower_bounds::ChainingLowerBounds,
    config::{ChainingLowerBoundConfig, InexactLowerBoundKind},
    costs::{AlignmentCosts, DEFAULT_CONFIG_TSA},
    dna_rc_fn,
};
use lib_tsalign::{
    alignment_geometry::{AlignmentCoordinates, AlignmentRange},
    alignment_result::AlignmentResult,
    config::TemplateSwitchConfig,
    costs::U64Cost,
};
use log::info;
use pyo3::{
    exceptions::{PyFileNotFoundError, PyIOError, PyValueError},
    prelude::*,
};
use pythonize::{depythonize, pythonize};

use crate::{
    TSPairwiseAlignment, gap_characters_to_ascii, parse_arena_index, py_to_str, search_error,
};

type LowerBounds = Arc<ChainingLowerBounds<U32Cost>>;

/// The parameters that determine the lower bounds besides the costs: `k` and `max_n`.
type LowerBoundsKey = (u32, usize);

#[pyclass(name = "ChainAligner")]
pub struct TSChainAligner {
    alignment_costs: AlignmentCosts<U32Cost>,
    /// The anchor size, or `None` to infer it from the sequence lengths.
    k: Option<u32>,
    max_anchor_mutations: u8,
    inexact_lower_bound: InexactLowerBoundKind,
    performance_parameters: AlignmentPerformanceParameters<U32Cost>,
    cache_directory: Option<PathBuf>,
    force_no_preprocessing: bool,
    /// The lower bounds that were computed or loaded so far.
    lower_bounds: Mutex<HashMap<LowerBoundsKey, LowerBounds>>,
    /// Held while loading or computing lower bounds, such that concurrent calls do not compute the same lower bounds twice.
    loading: Mutex<()>,
}

fn value_error(message: impl Into<String>) -> PyErr {
    PyValueError::new_err(message.into())
}

/// Converts a Python integer into the given type, raising a [`ValueError`](PyValueError) if it does not fit.
fn convert_integer<T: TryFrom<i64>>(value: i64, name: &str, expected: &str) -> PyResult<T> {
    T::try_from(value).map_err(|_| value_error(format!("{name} must be {expected}, got {value}")))
}

/// Parses an option the same way as the CLI parses it.
fn parse_value_enum<T: clap::ValueEnum>(value: &str, name: &str) -> PyResult<T> {
    T::from_str(value, false).map_err(|_| {
        let expected: Vec<_> = T::value_variants()
            .iter()
            .filter_map(|variant| variant.to_possible_value())
            .map(|value| format!("'{}'", value.get_name()))
            .collect();
        value_error(format!(
            "invalid {name} '{value}', expected one of {}",
            expected.join(", ")
        ))
    })
}

fn parse_inexact_lower_bound(value: &str) -> PyResult<InexactLowerBoundKind> {
    match value {
        "max-match-run" => Ok(InexactLowerBoundKind::MaxMatchRun),
        "windowed-min-mutations" => Err(value_error(
            "inexact_chaining_lower_bound 'windowed-min-mutations' is disabled because it is broken \
             (as in the CLI); use 'max-match-run'",
        )),
        value => Err(value_error(format!(
            "invalid inexact_chaining_lower_bound '{value}', expected 'max-match-run'"
        ))),
    }
}

/// Parses a `.tsa` cost configuration and converts it into the costs of tschainalign.
fn parse_costs(costs: &str) -> PyResult<AlignmentCosts<U32Cost>> {
    // Like the CLI, accept configs with and without `N` in the alphabet.
    // The costs of tschainalign do not depend on the alphabet of the config, so neither does the cache file name.
    let conversion_error = |error: lib_ts_chainalign::costs::CostConversionError| {
        value_error(format!(
            "The costs are not supported by tschainalign: {error}"
        ))
    };
    match TemplateSwitchConfig::<DnaAlphabet, U32Cost>::read_plain(costs.as_bytes()) {
        Ok(config) => AlignmentCosts::try_from(config).map_err(conversion_error),
        Err(dna_error) => {
            match TemplateSwitchConfig::<DnaAlphabetOrN, U32Cost>::read_plain(costs.as_bytes()) {
                Ok(config) => AlignmentCosts::try_from(config).map_err(conversion_error),
                Err(_) => Err(value_error(format!("Invalid costs: {dna_error}"))),
            }
        }
    }
}

impl TSChainAligner {
    fn config(&self, k: u32) -> ChainingLowerBoundConfig {
        ChainingLowerBoundConfig {
            anchor_k: k,
            max_anchor_mutations: self.max_anchor_mutations,
            inexact_lower_bound: self.inexact_lower_bound.clone(),
        }
    }

    /// Checks that `k` is supported, like the CLI.
    fn validate_k(&self, k: u32, k_given: bool) -> PyResult<()> {
        self.config(k).validate().map_err(|error| {
            let k_source = if k_given {
                "given with k"
            } else {
                "inferred from the sequence lengths and max_anchor_mutations; pass k to choose it explicitly"
            };
            value_error(format!(
                "Unsupported tschainalign parameters: {error} (k was {k_source})"
            ))
        })
    }

    fn cache_file(&self, k: u32, max_n: usize) -> Option<PathBuf> {
        self.cache_directory.as_ref().map(|cache_directory| {
            cache_file_path(
                &self.alignment_costs,
                cache_directory,
                k,
                self.max_anchor_mutations,
                &self.inexact_lower_bound,
                max_n,
            )
        })
    }

    fn create_cache_directory(&self) -> PyResult<()> {
        if let Some(cache_directory) = &self.cache_directory {
            std::fs::create_dir_all(cache_directory).map_err(|error| {
                PyIOError::new_err(format!(
                    "Creating cache directory {cache_directory:?}: {error}"
                ))
            })?;
        }
        Ok(())
    }

    fn read_cache_file(&self, cache_file: &Path, k: u32) -> PyResult<ChainingLowerBounds<U32Cost>> {
        info!("Loading preprocessed data from cache at {cache_file:?}");
        let lower_bounds = File::open(cache_file)
            .and_then(ChainingLowerBounds::read)
            .map_err(|error| {
                PyIOError::new_err(format!("Reading cache file {cache_file:?}: {error}"))
            })?;
        if lower_bounds.alignment_costs() != &self.alignment_costs
            || lower_bounds.anchor_k() != k
            || lower_bounds.max_anchor_mutations() != self.max_anchor_mutations
        {
            return Err(value_error(format!(
                "The cache file {cache_file:?} was preprocessed with different costs, k or max_anchor_mutations; \
                 remove it to preprocess again"
            )));
        }
        Ok(lower_bounds)
    }

    /// Computes the lower bounds and writes them into the cache directory, if any.
    fn compute_lower_bounds(&self, k: u32, max_n: usize) -> PyResult<ChainingLowerBounds<U32Cost>> {
        info!("Preprocessing for max_n = {max_n} and k = {k}...");
        let lower_bounds =
            lib_ts_chainalign::preprocess(max_n, self.config(k), self.alignment_costs.clone())
                .map_err(|error| {
                    value_error(format!("Unsupported tschainalign parameters: {error}"))
                })?;

        if let Some(cache_file) = self.cache_file(k, max_n) {
            self.create_cache_directory()?;
            info!("Storing preprocessed data into cache at {cache_file:?}");
            write_cache_file_atomically(&cache_file, |file| lower_bounds.write(file))
                .map_err(|error| PyIOError::new_err(error.to_string()))?;
        }
        Ok(lower_bounds)
    }

    fn cached_lower_bounds(&self, k: u32, max_n: usize) -> Option<LowerBounds> {
        self.lower_bounds
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .get(&(k, max_n))
            .cloned()
    }

    /// Returns the lower bounds for `k` and `max_n` from memory, from the cache directory, or by preprocessing.
    ///
    /// Must be called without holding the GIL.
    fn lower_bounds(&self, k: u32, max_n: usize) -> PyResult<LowerBounds> {
        if let Some(lower_bounds) = self.cached_lower_bounds(k, max_n) {
            return Ok(lower_bounds);
        }
        let _loading = self
            .loading
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        // Another thread may have loaded them while this one was waiting.
        if let Some(lower_bounds) = self.cached_lower_bounds(k, max_n) {
            return Ok(lower_bounds);
        }

        let cache_file = self.cache_file(k, max_n);
        let lower_bounds = match &cache_file {
            Some(cache_file) if cache_file.exists() => self.read_cache_file(cache_file, k)?,
            _ if self.force_no_preprocessing => {
                return Err(PyFileNotFoundError::new_err(match cache_file {
                    Some(cache_file) => format!(
                        "force_no_preprocessing is set, but no cache file found at {cache_file:?}"
                    ),
                    None => format!(
                        "force_no_preprocessing is set, but the lower bounds for k = {k} and max_n = {max_n} \
                         were not preprocessed, and there is no cache_directory"
                    ),
                }));
            }
            _ => self.compute_lower_bounds(k, max_n)?,
        };

        let lower_bounds = Arc::new(lower_bounds);
        self.lower_bounds
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .insert((k, max_n), lower_bounds.clone());
        Ok(lower_bounds)
    }

    /// Preprocesses the lower bounds for `k` and `max_n` if they are not in memory yet.
    ///
    /// Must be called without holding the GIL.
    fn preprocess_one(&self, k: u32, max_n: usize) -> PyResult<()> {
        let _loading = self
            .loading
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if self.cached_lower_bounds(k, max_n).is_some() {
            return Ok(());
        }
        let lower_bounds = Arc::new(self.compute_lower_bounds(k, max_n)?);
        self.lower_bounds
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .insert((k, max_n), lower_bounds);
        Ok(())
    }
}

#[pymethods]
impl TSChainAligner {
    #[new]
    #[pyo3(signature = (
        costs=None,
        k=None,
        max_anchor_mutations=0,
        inexact_chaining_lower_bound=None,
        max_chaining_successors=10,
        max_exact_cost_function_cost=0,
        chaining_closed_list="special",
        chaining_open_list="linear-heap",
        anchor_pruning=None,
        inexact_pruning_rules=None,
        exact_overlap_diagonal_distance=None,
        cache_directory=None,
        force_no_preprocessing=false,
        arena_index="u32",
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        costs: Option<&str>,
        k: Option<i64>,
        max_anchor_mutations: i64,
        inexact_chaining_lower_bound: Option<&str>,
        max_chaining_successors: i64,
        max_exact_cost_function_cost: i64,
        chaining_closed_list: &str,
        chaining_open_list: &str,
        anchor_pruning: Option<&str>,
        inexact_pruning_rules: Option<&str>,
        exact_overlap_diagonal_distance: Option<i64>,
        cache_directory: Option<PathBuf>,
        force_no_preprocessing: bool,
        arena_index: &str,
    ) -> PyResult<Self> {
        let arena_index = parse_arena_index(arena_index)?;
        let max_anchor_mutations: u8 = convert_integer(
            max_anchor_mutations,
            "max_anchor_mutations",
            "between 0 and 255",
        )?;
        let k: Option<u32> = k
            .map(|k| convert_integer(k, "k", "a non-negative integer"))
            .transpose()?;
        let max_successors: usize = convert_integer(
            max_chaining_successors,
            "max_chaining_successors",
            "a non-negative integer",
        )?;
        let max_exact_cost_function_cost: u32 = convert_integer(
            max_exact_cost_function_cost,
            "max_exact_cost_function_cost",
            "a non-negative 32-bit integer",
        )?;
        let closed_list: ChainingClosedList =
            parse_value_enum(chaining_closed_list, "chaining_closed_list")?;
        let open_list: ChainingOpenList =
            parse_value_enum(chaining_open_list, "chaining_open_list")?;

        // Like the CLI, reject options that only affect inexact anchors if the anchors are exact,
        // such that they are not silently ignored.
        if max_anchor_mutations == 0 {
            for (name, given) in [
                (
                    "inexact_chaining_lower_bound",
                    inexact_chaining_lower_bound.is_some(),
                ),
                ("anchor_pruning", anchor_pruning.is_some()),
                ("inexact_pruning_rules", inexact_pruning_rules.is_some()),
                (
                    "exact_overlap_diagonal_distance",
                    exact_overlap_diagonal_distance.is_some(),
                ),
            ] {
                if given {
                    return Err(value_error(format!(
                        "{name} has no effect with exact anchors (max_anchor_mutations=0); \
                         remove {name} or allow mutations in anchors"
                    )));
                }
            }
        }

        let inexact_lower_bound = inexact_chaining_lower_bound
            .map(parse_inexact_lower_bound)
            .transpose()?
            .unwrap_or(InexactLowerBoundKind::MaxMatchRun);
        let anchor_pruning: AnchorPruning = anchor_pruning
            .map(|value| parse_value_enum(value, "anchor_pruning"))
            .transpose()?
            .unwrap_or_default();
        if anchor_pruning != AnchorPruning::Inexact
            && (inexact_pruning_rules.is_some() || exact_overlap_diagonal_distance.is_some())
        {
            return Err(value_error(
                "inexact_pruning_rules and exact_overlap_diagonal_distance require anchor_pruning='inexact'",
            ));
        }
        let default_inexact_pruning = InexactPruningParameters::default();
        let inexact_pruning = InexactPruningParameters {
            rules: inexact_pruning_rules
                .map(|rules| rules.parse::<InexactPruningRules>())
                .transpose()
                .map_err(value_error)?
                .unwrap_or(default_inexact_pruning.rules),
            exact_overlap_diagonal_distance: exact_overlap_diagonal_distance
                .map(|distance| {
                    convert_integer(
                        distance,
                        "exact_overlap_diagonal_distance",
                        "a non-negative integer",
                    )
                })
                .transpose()?
                .unwrap_or(default_inexact_pruning.exact_overlap_diagonal_distance),
        };

        let alignment_costs = parse_costs(costs.unwrap_or(DEFAULT_CONFIG_TSA))?;

        let aligner = Self {
            alignment_costs,
            k,
            max_anchor_mutations,
            inexact_lower_bound,
            performance_parameters: AlignmentPerformanceParameters {
                max_successors,
                max_exact_cost_function_cost: max_exact_cost_function_cost.into(),
                closed_list,
                open_list,
                anchor_pruning,
                inexact_pruning,
                arena_index,
            },
            cache_directory,
            force_no_preprocessing,
            lower_bounds: Default::default(),
            loading: Default::default(),
        };
        if let Some(k) = k {
            aligner.validate_k(k, true)?;
        }
        Ok(aligner)
    }

    /// Preprocess the lower bounds for sequences of length up to `max_length`, like `tsalign preprocess`.
    fn preprocess(&self, py: Python<'_>, max_length: i64) -> PyResult<()> {
        let max_length: usize = convert_integer(max_length, "max_length", "at least 1")?;
        if max_length == 0 {
            return Err(value_error("max_length must be at least 1"));
        }

        // Infer and validate all ks before preprocessing anything.
        let mut ks_per_max_n = Vec::new();
        let mut current_max_n = infer_max_n(max_length, max_length);
        while current_max_n >= 1 {
            let next_max_n = current_max_n / 2;
            let ks: BTreeSet<_> = match self.k {
                Some(k) => [k].into(),
                None => (next_max_n + 1..=current_max_n)
                    .map(|max_n| infer_k(max_n, max_n, self.max_anchor_mutations))
                    .collect(),
            };
            for &k in &ks {
                self.validate_k(k, self.k.is_some())?;
            }
            ks_per_max_n.push((current_max_n, ks));
            current_max_n = next_max_n;
        }

        py.detach(|| {
            for (max_n, ks) in ks_per_max_n {
                for k in ks {
                    self.preprocess_one(k, max_n)?;
                }
            }
            Ok(())
        })
    }

    #[pyo3(signature = (reference, query, reference_name="reference", query_name="query", reference_start=None, reference_limit=None, query_start=None, query_limit=None, gap_characters="", cost_limit=None, iteration_limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn align(
        &self,
        py: Python<'_>,
        reference: Bound<'_, PyAny>,
        query: Bound<'_, PyAny>,
        reference_name: &str,
        query_name: &str,
        reference_start: Option<usize>,
        reference_limit: Option<usize>,
        query_start: Option<usize>,
        query_limit: Option<usize>,
        gap_characters: &str,
        cost_limit: Option<i64>,
        iteration_limit: Option<i64>,
    ) -> PyResult<Option<TSPairwiseAlignment>> {
        let gap_characters = gap_characters_to_ascii(gap_characters)?;
        let remove_gaps = |sequence: Vec<u8>| -> Vec<u8> {
            sequence
                .into_iter()
                .filter(|c| !gap_characters.contains(c))
                .collect()
        };
        let reference = remove_gaps(py_to_str(reference)?);
        let query = remove_gaps(py_to_str(query)?);

        let cost_limit: Option<u64> = cost_limit
            .map(|limit| convert_integer(limit, "cost_limit", "a non-negative integer"))
            .transpose()?;
        let iteration_limit: Option<u64> = iteration_limit
            .map(|limit| convert_integer(limit, "iteration_limit", "a non-negative integer"))
            .transpose()?;

        let range = AlignmentRange::new_offset_limit(
            AlignmentCoordinates::new(reference_start.unwrap_or(0), query_start.unwrap_or(0)),
            AlignmentCoordinates::new(
                reference_limit.unwrap_or(reference.len()),
                query_limit.unwrap_or(query.len()),
            ),
        );
        // Check the inputs before preprocessing.
        lib_ts_chainalign::validate_sequences(&reference, &query, &range)
            .map_err(|error| value_error(error.to_string()))?;
        if reference.is_empty() && query.is_empty() {
            return Err(value_error(
                "tschainalign needs at least one non-empty sequence",
            ));
        }

        let max_n = infer_max_n(reference.len(), query.len());
        let k = match self.k {
            Some(k) => k,
            None => {
                let k = infer_k(reference.len(), query.len(), self.max_anchor_mutations);
                self.validate_k(k, false)?;
                k
            }
        };

        let limits = AlignmentLimits {
            // Chain costs are at most `u32::MAX`, so a larger limit is no limit.
            cost_limit: cost_limit.map(|limit| u32::try_from(limit).unwrap_or(u32::MAX).into()),
            iteration_limit,
        };

        let reference_name = reference_name.to_string();
        let query_name = query_name.to_string();
        let result = py.detach(|| {
            let lower_bounds = self.lower_bounds(k, max_n)?;
            lib_ts_chainalign::align_with_limits(
                reference,
                query,
                range,
                &self.performance_parameters,
                &limits,
                &dna_rc_fn,
                &reference_name,
                &query_name,
                &lower_bounds,
            )
            .map_err(|error| match error {
                ChainAlignError::Search(error) => {
                    search_error(&error, self.performance_parameters.arena_index)
                }
                error => value_error(error.to_string()),
            })
        })?;

        match result {
            AlignmentResult::WithTarget { .. } => {
                // tschainalign computes 32-bit costs, while the result type of the bindings has 64-bit costs.
                // Both serialise the costs as plain integers.
                let result: AlignmentResult<U64Cost> = depythonize(&pythonize(py, &result)?)?;
                Ok(Some(TSPairwiseAlignment { result }))
            }
            AlignmentResult::WithoutTarget { .. } => Ok(None),
        }
    }
}
