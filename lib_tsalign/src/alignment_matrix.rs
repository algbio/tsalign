use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::cost::AStarCost;
use index::{
    AlignmentMatrixIndex,
    iterators::{
        AlignmentMatrixInnerIterator, AlignmentMatrixQueryIterator,
        AlignmentMatrixReferenceIterator,
    },
};
use ndarray::Array2;
use num_traits::Bounded;

use crate::alignment_configuration::AlignmentConfiguration;

pub mod index;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentMatrix<Cost> {
    matrix: Array2<AlignmentMatrixEntry<Cost>>,
    configuration: AlignmentConfiguration<Cost>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentMatrixEntry<Cost> {
    pub cost: Cost,
    pub alignment_type: BaseAlignmentType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseAlignmentType {
    /// Marks the matrix origin at [0, 0].
    None,
    Insertion,
    Deletion,
    Match,
    Substitution,
}

impl<Cost: AStarCost> AlignmentMatrix<Cost> {
    pub fn new(
        configuration: AlignmentConfiguration<Cost>,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self {
            matrix: Array2::default((reference_length + 1, query_length + 1)),
            configuration,
        }
    }

    pub fn reference_index_iter(&self, query_index: usize) -> AlignmentMatrixReferenceIterator {
        AlignmentMatrixReferenceIterator::new(query_index, self.matrix.dim().0)
    }

    pub fn query_index_iter(&self, reference_index: usize) -> AlignmentMatrixQueryIterator {
        AlignmentMatrixQueryIterator::new(reference_index, self.matrix.dim().1)
    }

    pub fn inner_index_iter(&self) -> AlignmentMatrixInnerIterator {
        AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(
            self.matrix.dim().0,
            self.matrix.dim().1,
        ))
    }

    pub fn align<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SequenceType,
        query: &SequenceType,
    ) -> Cost {
        self.initialise();
        self.align_inner(reference, query);
        self.matrix[[self.matrix.dim().0 - 1, self.matrix.dim().1 - 1]].cost
    }

    fn initialise(&mut self) {
        // Initialise matrix origin.
        self.matrix[[0, 0]].cost = Cost::zero();
        self.matrix[[0, 0]].alignment_type = BaseAlignmentType::None;

        // Initialise matrix edges.
        for index in self.reference_index_iter(0).skip(1) {
            self.set_deletion_cost(index);
        }
        for index in self.query_index_iter(0).skip(1) {
            self.set_insertion_cost(index);
        }
    }

    fn align_inner<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SequenceType,
        query: &SequenceType,
    ) {
        for index in self.inner_index_iter() {
            self.set_min_cost(index, reference, query);
        }
    }

    fn set_insertion_cost(&mut self, index: AlignmentMatrixIndex) {
        self.matrix[index] = self.compute_insertion_entry(index);
    }

    fn set_deletion_cost(&mut self, index: AlignmentMatrixIndex) {
        self.matrix[index] = self.compute_deletion_entry(index);
    }

    fn set_min_cost<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        index: AlignmentMatrixIndex,
        reference: &SequenceType,
        query: &SequenceType,
    ) {
        // Handle matches and substitutions.
        let mut entry = self.compute_match_or_substitution_entry(index, reference, query);

        // Handle insertions.
        let insertion_entry = self.compute_insertion_entry(index);
        if insertion_entry < entry {
            entry = insertion_entry;
        }

        // Handle deletions.
        let deletion_entry = self.compute_deletion_entry(index);
        if deletion_entry < entry {
            entry = deletion_entry;
        }

        self.matrix[index] = entry;
    }

    fn compute_insertion_entry(&self, index: AlignmentMatrixIndex) -> AlignmentMatrixEntry<Cost> {
        let alignment_type = BaseAlignmentType::Insertion;
        let predecessor_cost = self.matrix[index.predecessor(alignment_type)].cost;

        AlignmentMatrixEntry {
            cost: predecessor_cost + self.configuration.cost(alignment_type),
            alignment_type,
        }
    }

    fn compute_deletion_entry(&self, index: AlignmentMatrixIndex) -> AlignmentMatrixEntry<Cost> {
        let alignment_type = BaseAlignmentType::Deletion;
        let predecessor_cost = self.matrix[index.predecessor(alignment_type)].cost;

        AlignmentMatrixEntry {
            cost: predecessor_cost + self.configuration.cost(alignment_type),
            alignment_type,
        }
    }

    fn compute_match_or_substitution_entry<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        index: AlignmentMatrixIndex,
        reference: &SequenceType,
        query: &SequenceType,
    ) -> AlignmentMatrixEntry<Cost> {
        let alignment_type = if reference[index.reference_index - 1] == query[index.query_index - 1]
        {
            BaseAlignmentType::Match
        } else {
            BaseAlignmentType::Substitution
        };
        let precedessor_cost = self.matrix[index.match_or_substitution_predecessor()].cost;

        AlignmentMatrixEntry {
            cost: precedessor_cost + self.configuration.cost(alignment_type),
            alignment_type,
        }
    }

    #[cfg(test)]
    fn manual_debug_fill(&mut self, entries: impl IntoIterator<Item = AlignmentMatrixEntry<Cost>>) {
        let mut entries = entries.into_iter();
        for index in self.inner_index_iter() {
            self.matrix[index] = entries.next().unwrap();
        }
        assert!(entries.next().is_none());
    }
}

impl<Cost: PartialOrd> PartialOrd for AlignmentMatrixEntry<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.cost.partial_cmp(&other.cost) {
            ord @ (Some(core::cmp::Ordering::Greater | core::cmp::Ordering::Less) | None) => ord,
            ord @ Some(core::cmp::Ordering::Equal) => {
                if self.alignment_type == other.alignment_type {
                    ord
                } else {
                    None
                }
            }
        }
    }
}

impl<Cost: Bounded> Default for AlignmentMatrixEntry<Cost> {
    fn default() -> Self {
        Self {
            cost: Cost::max_value(),
            alignment_type: BaseAlignmentType::None,
        }
    }
}

impl<Cost: AStarCost> core::fmt::Display for AlignmentMatrix<Cost> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut cost_column_widths = vec![0; self.matrix.dim().1];
        for reference_index in 0..self.matrix.dim().0 {
            for (query_index, cost_column_width) in cost_column_widths.iter_mut().enumerate() {
                let cost = self.matrix[[reference_index, query_index]].cost.as_u64();
                let mut local_cost_column_width = 1;

                let mut limit = 1;
                for _ in 1..=10 {
                    limit *= 10;
                    if cost >= limit {
                        local_cost_column_width += 1;
                    }
                }

                *cost_column_width = local_cost_column_width.max(*cost_column_width);
            }
        }

        for reference_index in 0..self.matrix.dim().0 {
            write!(f, "[ ")?;
            #[allow(clippy::needless_range_loop)]
            for query_index in 0..self.matrix.dim().1 {
                write!(
                    f,
                    "{: >width$}",
                    self.matrix[[reference_index, query_index]].cost,
                    width = cost_column_widths[query_index],
                )?;
                write!(
                    f,
                    "{} ",
                    match self.matrix[[reference_index, query_index]].alignment_type {
                        BaseAlignmentType::None => "N",
                        BaseAlignmentType::Insertion => "I",
                        BaseAlignmentType::Deletion => "D",
                        BaseAlignmentType::Match => "M",
                        BaseAlignmentType::Substitution => "S",
                    }
                )?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::{
            alphabets::dna_alphabet::DnaAlphabet, bit_vec_sequence_store::BitVectorSequenceStore,
        },
        interface::sequence_store::SequenceStore,
    };
    use generic_a_star::cost::U64Cost;
    use traitsequence::interface::Sequence;

    use crate::{
        alignment_configuration::AlignmentConfiguration,
        alignment_matrix::{AlignmentMatrixEntry, BaseAlignmentType},
    };

    use super::AlignmentMatrix;

    #[test]
    fn simple_alignments() {
        let mut sequence_store = BitVectorSequenceStore::<DnaAlphabet>::new();
        let reference = sequence_store.add_from_slice_u8(b"ACG").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCG").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::<U64Cost>::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 3u64.into());

        let mut manual_matrix = matrix.clone();
        manual_matrix.manual_debug_fill(
            [
                (0u64, BaseAlignmentType::Match),
                (3, BaseAlignmentType::Deletion),
                (6, BaseAlignmentType::Deletion),
                (3, BaseAlignmentType::Insertion),
                (0, BaseAlignmentType::Match),
                (3, BaseAlignmentType::Deletion),
                (6, BaseAlignmentType::Insertion),
                (3, BaseAlignmentType::Match),
                (2, BaseAlignmentType::Substitution),
                (9, BaseAlignmentType::Insertion),
                (6, BaseAlignmentType::Insertion),
                (3, BaseAlignmentType::Match),
            ]
            .into_iter()
            .map(|(cost, alignment_type)| AlignmentMatrixEntry {
                cost: cost.into(),
                alignment_type,
            }),
        );
        assert_eq!(
            matrix, manual_matrix,
            "matrix:\n{matrix}\nmanual_matrix:\n{manual_matrix}"
        );

        let reference = sequence_store.add_from_slice_u8(b"ACCCGT").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCGT").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::<U64Cost>::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 3u64.into());

        let reference = sequence_store.add_from_slice_u8(b"ACGCCCCCT").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCCCCGCT").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::<U64Cost>::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 4u64.into());
    }
}
