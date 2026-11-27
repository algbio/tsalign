use std::io;

use compact_genome::implementation::alphabets::dna_alphabet_or_n::DnaAlphabetOrN;
use generic_a_star::AStarError;
use lib_tsalign::{
    Aligner,
    alignment_geometry::{AlignmentCoordinates, AlignmentRange},
    alignment_result::{AlignmentResult, TSMUncertaintyRangeExtensionMode},
    arena_index::ArenaIndexType,
    costs::U64Cost,
    error::AlignerInputError,
};
use lib_tsshow::plain_text::show_template_switches;
use pyo3::{
    create_exception,
    exceptions::{PyMemoryError, PyRuntimeError, PyValueError},
    prelude::*,
    types::PyDict,
};
use pythonize::{depythonize, pythonize};

mod chain;

create_exception!(
    _tsalign,
    ArenaFullError,
    PyMemoryError,
    "The A* search needs more nodes at the same time than its arena index type allows.\n\n\
     Raised by the ``align`` methods. Create the aligner with ``arena_index=\"u64\"`` to allow more nodes."
);

/// Parses the `arena_index` argument of the aligners.
fn parse_arena_index(arena_index: &str) -> PyResult<ArenaIndexType> {
    match arena_index {
        "u32" => Ok(ArenaIndexType::U32),
        "u64" => Ok(ArenaIndexType::U64),
        _ => Err(PyValueError::new_err(format!(
            "invalid arena_index '{arena_index}', expected one of 'u32', 'u64'"
        ))),
    }
}

/// Converts the error of a failed A* search into a Python exception.
///
/// A full node arena raises [`ArenaFullError`], with a hint to use a larger index type if there is one.
fn search_error(error: &AStarError, arena_index: ArenaIndexType) -> PyErr {
    match (error, arena_index) {
        (AStarError::ArenaFull(_), ArenaIndexType::U32) => ArenaFullError::new_err(format!(
            "{error}; create the aligner with arena_index=\"u64\" to allow more nodes \
             (this needs a bit more memory per node)"
        )),
        (AStarError::ArenaFull(_), _) => ArenaFullError::new_err(error.to_string()),
        _ => PyRuntimeError::new_err(error.to_string()),
    }
}

#[pyclass]
struct TSPairwiseAlignment {
    result: AlignmentResult<U64Cost>,
}

#[pymethods]
impl TSPairwiseAlignment {
    fn viz_template_switches(&self) -> PyResult<()> {
        show_template_switches(io::stdout(), &self.result, &None)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn stats<'a>(&'a self, py: Python<'a>) -> PyResult<Bound<'a, PyAny>> {
        Ok(pythonize(py, self.result.statistics())?)
    }

    fn cigar(&self) -> Option<String> {
        match &self.result {
            AlignmentResult::WithTarget { alignment, .. } => Some(alignment.cigar()),
            AlignmentResult::WithoutTarget { .. } => None,
        }
    }

    fn alignments<'a>(&'a self, py: Python<'a>) -> PyResult<Option<Bound<'a, PyAny>>> {
        match &self.result {
            AlignmentResult::WithTarget { alignment, .. } => {
                let mut container = Vec::new();
                alignment.iter_compact().for_each(|e| container.push(e));
                Ok(Some(pythonize(py, &container)?))
            }
            AlignmentResult::WithoutTarget { .. } => Ok(None),
        }
    }
}

fn py_to_str(o: Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    let str = o.str()?.to_str()?.as_bytes().to_vec();
    Ok(str)
}

fn gap_characters_to_ascii(gap_characters: &str) -> PyResult<Vec<u8>> {
    gap_characters
        .chars()
        .map(|c| {
            u8::try_from(c).ok().filter(u8::is_ascii).ok_or_else(|| {
                PyValueError::new_err("gap_characters must contain only ASCII characters")
            })
        })
        .collect()
}

#[pyclass(name = "Aligner")]
struct TSAligner {
    aligner: Aligner<DnaAlphabetOrN>,
    arena_index: ArenaIndexType,
}

#[pymethods]
impl TSAligner {
    #[new]
    #[pyo3(signature = (**kwargs))]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let Some(kwargs) = kwargs else {
            return Ok(Self {
                aligner: Aligner::new(),
                arena_index: ArenaIndexType::default(),
            });
        };

        let costs_str = if let Some(costs) = kwargs.get_item("costs")? {
            let py_str: &str = costs.extract()?;
            kwargs.del_item("costs")?;
            Some(py_str.to_string())
        } else {
            None
        };

        let arena_index = if let Some(arena_index) = kwargs.get_item("arena_index")? {
            let arena_index = parse_arena_index(arena_index.extract()?)?;
            kwargs.del_item("arena_index")?;
            arena_index
        } else {
            ArenaIndexType::default()
        };

        let mut aligner: Aligner<DnaAlphabetOrN> = depythonize(kwargs)?;
        aligner.set_arena_index(arena_index);

        if let Some(costs_str) = costs_str {
            aligner
                .set_costs_parse(&costs_str)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        }

        Ok(Self {
            aligner,
            arena_index,
        })
    }

    /// Align two sequences, accounting for template switches
    ///
    /// The function takes a reference and a query string, and performs a global alignment on both. The output alignment may contain (short-range) template switches.
    /// Optionally, settings can be specified on this aligner.
    #[pyo3(signature = (reference, query, reference_name="reference", query_name="query", reference_start=None, reference_limit=None, query_start=None, query_limit=None, gap_characters="", cost_limit=None, memory_limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn align(
        &self,
        reference: Bound<'_, PyAny>, // Accepting PyAny instead of PyString to allow using e.g. `Bio.Seq` types and alike. String representation will be used.
        query: Bound<'_, PyAny>,
        reference_name: &str,
        query_name: &str,
        reference_start: Option<usize>,
        reference_limit: Option<usize>,
        query_start: Option<usize>,
        query_limit: Option<usize>,
        gap_characters: &str,
        cost_limit: Option<u64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Option<TSPairwiseAlignment>> {
        let reference = py_to_str(reference)?;
        let query = py_to_str(query)?;

        let gap_characters = gap_characters_to_ascii(gap_characters)?;

        // The aligner removes the gap characters before applying the ranges,
        // so the default limits are the lengths without gap characters.
        let length_without_gaps = |sequence: &[u8]| {
            sequence
                .iter()
                .filter(|c| !gap_characters.contains(c))
                .count()
        };
        let reference_start = reference_start.unwrap_or(0);
        let reference_limit = reference_limit.unwrap_or_else(|| length_without_gaps(&reference));
        let query_start = query_start.unwrap_or(0);
        let query_limit = query_limit.unwrap_or_else(|| length_without_gaps(&query));
        let ranges = AlignmentRange::new_offset_limit(
            AlignmentCoordinates::new(reference_start, query_start),
            AlignmentCoordinates::new(reference_limit, query_limit),
        );

        let result = self
            .aligner
            .align(
                reference_name,
                &reference,
                query_name,
                &query,
                Some(ranges),
                &gap_characters,
                cost_limit,
                memory_limit,
                true,
                TSMUncertaintyRangeExtensionMode::EqualCostIgnoreGeometry,
            )
            .map_err(|error| match error {
                AlignerInputError::Search(error) => search_error(&error, self.arena_index),
                error => PyValueError::new_err(error.to_string()),
            })?;

        match result {
            result @ AlignmentResult::WithTarget { .. } => {
                let ts_alignment = TSPairwiseAlignment { result };
                Ok(Some(ts_alignment))
            }
            AlignmentResult::WithoutTarget { .. } => Ok(None),
        }
    }
}

/// Bindings for the `lib_tsalign` library.
#[pymodule(name = "_tsalign")]
fn tsalign(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_class::<TSPairwiseAlignment>()?;
    m.add_class::<TSAligner>()?;
    m.add_class::<chain::TSChainAligner>()?;
    m.add("ArenaFullError", m.py().get_type::<ArenaFullError>())?;
    m.add(
        "CHAIN_ALIGNER_DEFAULT_COSTS",
        lib_ts_chainalign::costs::DEFAULT_CONFIG_TSA,
    )?;
    Ok(())
}
