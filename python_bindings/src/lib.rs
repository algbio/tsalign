use std::io;

use compact_genome::implementation::alphabets::dna_alphabet_or_n::DnaAlphabetOrN;
use lib_tsalign::{
    a_star_aligner::{
        alignment_geometry::{AlignmentCoordinates, AlignmentRange},
        alignment_result::AlignmentResult,
        configurable_a_star_align::Aligner,
        template_switch_distance::AlignmentType,
    },
    costs::U64Cost,
};
use lib_tsshow::plain_text::show_template_switches;
use pyo3::{exceptions::PyRuntimeError, prelude::*, types::PyDict};
use pythonize::{depythonize, pythonize};

#[pyclass]
struct TSPairwiseAlignment {
    result: AlignmentResult<AlignmentType, U64Cost>,
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

#[pyclass(name = "Aligner")]
struct TSAligner {
    aligner: Aligner<DnaAlphabetOrN>,
}

#[pymethods]
impl TSAligner {
    #[new]
    #[pyo3(signature = (**kwargs))]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let Some(kwargs) = kwargs else {
            return Ok(Self {
                aligner: Aligner::new(),
            });
        };

        let costs_str = if let Some(costs) = kwargs.get_item("costs")? {
            let py_str: &str = costs.extract()?;
            kwargs.del_item("costs")?;
            Some(py_str.to_string())
        } else {
            None
        };

        let mut aligner: Aligner<DnaAlphabetOrN> = depythonize(kwargs)?;

        if let Some(costs_str) = costs_str {
            aligner
                .set_costs_parse(&costs_str)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        }

        Ok(Self { aligner })
    }

    /// Align two sequences, accounting for template switches
    ///
    /// The function takes a reference and a query string, and performs a global alignment on both. The output alignment may contain (short-range) template switches.
    /// Optionally, settings can be specified on this aligner.
    #[pyo3(signature = (reference, query, reference_name="reference", query_name="query", reference_start=None, reference_limit=None, query_start=None, query_limit=None, cost_limit=None, memory_limit=None))]
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
        cost_limit: Option<u64>,
        memory_limit: Option<usize>,
    ) -> PyResult<Option<TSPairwiseAlignment>> {
        let reference = py_to_str(reference)?;
        let query = py_to_str(query)?;

        let reference_start = reference_start.unwrap_or(0);
        let reference_limit = reference_limit.unwrap_or(reference.len());
        let query_start = query_start.unwrap_or(0);
        let query_limit = query_limit.unwrap_or(query.len());
        let ranges = AlignmentRange::new_offset_limit(
            AlignmentCoordinates::new(reference_start, query_start),
            AlignmentCoordinates::new(reference_limit, query_limit),
        );

        let result = self.aligner.align(
            reference_name,
            &reference,
            query_name,
            &query,
            Some(ranges),
            cost_limit,
            memory_limit,
            true,
        );

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
#[pymodule]
fn tsalign(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_class::<TSPairwiseAlignment>()?;
    m.add_class::<TSAligner>()?;
    Ok(())
}
