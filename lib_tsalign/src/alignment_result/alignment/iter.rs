use std::ops::{Deref, DerefMut};

use crate::a_star_aligner::alignment_result::IAlignmentType;

pub struct CompactAlignmentIter<'alignment, AlignmentType> {
    alignment: &'alignment [(usize, AlignmentType)],
    front_multiplicity: usize,
    back_multiplicity: usize,
}

pub struct CompactAlignmentIterCloned<'alignment, AlignmentType: IAlignmentType + Clone> {
    iter: CompactAlignmentIter<'alignment, AlignmentType>,
}

pub struct FlatAlignmentIter<'alignment, AlignmentType> {
    iter: CompactAlignmentIter<'alignment, AlignmentType>,
    length: usize,
}

pub struct FlatAlignmentIterCloned<'alignment, AlignmentType: IAlignmentType + Clone> {
    iter: FlatAlignmentIter<'alignment, AlignmentType>,
}

pub struct PeekFrontMultiplicityMut<
    'iter,
    'alignment: 'iter,
    AlignmentType: IAlignmentType,
    LengthDecreaser: PeekLengthDecreaser,
> {
    iter: &'iter mut CompactAlignmentIter<'alignment, AlignmentType>,
    initial_multiplicity: usize,
    length_decreaser: LengthDecreaser,
}

pub struct PeekBackMultiplicityMut<
    'iter,
    'alignment: 'iter,
    AlignmentType: IAlignmentType,
    LengthDecreaser: PeekLengthDecreaser,
> {
    iter: &'iter mut CompactAlignmentIter<'alignment, AlignmentType>,
    initial_multiplicity: usize,
    length_decreaser: LengthDecreaser,
}

pub trait PeekLengthDecreaser {
    fn decrease_length(&mut self, decrement: usize);
}

pub struct NoLengthDecreaser;

pub struct RefMutLengthDecreaser<'length> {
    length: &'length mut usize,
}

//////////////////////////////////////////////
// CompactAlignmentIter //////////////////////
//////////////////////////////////////////////

impl<'alignment, AlignmentType: IAlignmentType> CompactAlignmentIter<'alignment, AlignmentType> {
    /// Create a new iterator over the given alignment.
    pub(crate) fn new(alignment: &'alignment [(usize, AlignmentType)]) -> Self {
        let mut result = Self {
            alignment,
            front_multiplicity: alignment
                .first()
                .map(|(multiplicity, alignment_type)| {
                    if alignment_type.is_repeatable() {
                        *multiplicity
                    } else {
                        1.min(*multiplicity)
                    }
                })
                .unwrap_or(0),
            back_multiplicity: alignment
                .last()
                .map(|(multiplicity, alignment_type)| {
                    if alignment_type.is_repeatable() {
                        *multiplicity
                    } else {
                        1.min(*multiplicity)
                    }
                })
                .unwrap_or(0),
        };
        result.advance_zeroes();
        result
    }

    /// Peek at the first alignment in the iterator.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front(&self) -> Option<(usize, &'alignment AlignmentType)> {
        if self.front_multiplicity == 0 {
            return None;
        }

        self.alignment
            .first()
            .map(|(_, alignment_type)| (self.front_multiplicity, alignment_type))
    }

    /// Peek at the last alignment in the iterator.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back(&self) -> Option<(usize, &'alignment AlignmentType)> {
        if self.back_multiplicity == 0 {
            return None;
        }

        self.alignment
            .last()
            .map(|(_, alignment_type)| (self.back_multiplicity, alignment_type))
    }

    /// Peek at the first alignment in the iterator and clone it.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_cloned(&self) -> Option<(usize, AlignmentType)>
    where
        AlignmentType: IAlignmentType + Clone,
    {
        self.peek_front()
            .map(|(multiplicity, alignment_type)| (multiplicity, alignment_type.clone()))
    }

    /// Peek at the last alignment in the iterator and clone it.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_cloned(&self) -> Option<(usize, AlignmentType)>
    where
        AlignmentType: IAlignmentType + Clone,
    {
        self.peek_back()
            .map(|(multiplicity, alignment_type)| (multiplicity, alignment_type.clone()))
    }

    /// Mutably peek at the multiplicity of the first alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, NoLengthDecreaser>> {
        (self.front_multiplicity > 0)
            .then(|| PeekFrontMultiplicityMut::new(self, NoLengthDecreaser))
    }

    /// Mutably peek at the multiplicity of the last alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, NoLengthDecreaser>> {
        (self.back_multiplicity > 0).then(|| PeekBackMultiplicityMut::new(self, NoLengthDecreaser))
    }

    fn advance_zeroes(&mut self) {
        self.advance_zeroes_front();
        self.advance_zeroes_back();
    }

    fn advance_zeroes_front(&mut self) {
        while self.front_multiplicity == 0 && !self.alignment.is_empty() {
            self.pop_front();
        }
    }

    fn advance_zeroes_back(&mut self) {
        while self.back_multiplicity == 0 && !self.alignment.is_empty() {
            self.pop_back();
        }
    }

    fn pop_front(&mut self) {
        debug_assert_eq!(self.front_multiplicity, 0);

        if !self.alignment.is_empty() {
            self.alignment = &self.alignment[1..];
            self.front_multiplicity = self
                .alignment
                .first()
                .map(|(multiplicity, alignment_type)| {
                    if alignment_type.is_repeatable() {
                        *multiplicity
                    } else {
                        1.min(*multiplicity)
                    }
                })
                .unwrap_or(0);
        }
    }

    fn pop_back(&mut self) {
        debug_assert_eq!(self.back_multiplicity, 0);

        if !self.alignment.is_empty() {
            self.alignment = &self.alignment[..self.alignment.len() - 1];
            self.back_multiplicity = self
                .alignment
                .last()
                .map(|(multiplicity, alignment_type)| {
                    if alignment_type.is_repeatable() {
                        *multiplicity
                    } else {
                        1.min(*multiplicity)
                    }
                })
                .unwrap_or(0);
        }
    }
}

impl<'alignment, AlignmentType: IAlignmentType> Iterator
    for CompactAlignmentIter<'alignment, AlignmentType>
{
    type Item = (usize, &'alignment AlignmentType);

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.peek_front();
        self.front_multiplicity = 0;
        self.advance_zeroes_front();
        result
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.alignment.len(), Some(self.alignment.len()))
    }
}

impl<AlignmentType: IAlignmentType> DoubleEndedIterator
    for CompactAlignmentIter<'_, AlignmentType>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        let result = self.peek_back();
        self.back_multiplicity = 0;
        self.advance_zeroes_back();
        result
    }
}

impl<AlignmentType: IAlignmentType> ExactSizeIterator for CompactAlignmentIter<'_, AlignmentType> {}

//////////////////////////////////////////////
// CompactAlignmentIterCloned ////////////////
//////////////////////////////////////////////

impl<'alignment, AlignmentType: IAlignmentType + Clone>
    CompactAlignmentIterCloned<'alignment, AlignmentType>
{
    pub fn new(alignment: &'alignment [(usize, AlignmentType)]) -> Self {
        Self {
            iter: CompactAlignmentIter::new(alignment),
        }
    }

    /// Peek at the first alignment in the iterator.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front(&self) -> Option<(usize, &'alignment AlignmentType)> {
        self.iter.peek_front()
    }

    /// Peek at the last alignment in the iterator.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back(&self) -> Option<(usize, &'alignment AlignmentType)> {
        self.iter.peek_back()
    }

    /// Peek at the first alignment in the iterator and clone it.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_cloned(&self) -> Option<(usize, AlignmentType)> {
        self.iter.peek_front_cloned()
    }

    /// Peek at the last alignment in the iterator and clone it.
    ///
    /// Returns the remaining multiplicity along with the alignment type.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_cloned(&self) -> Option<(usize, AlignmentType)> {
        self.iter.peek_back_cloned()
    }

    /// Mutably peek at the multiplicity of the first alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, NoLengthDecreaser>> {
        self.iter.peek_front_multiplicity_mut()
    }

    /// Mutably peek at the multiplicity of the last alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, NoLengthDecreaser>> {
        self.iter.peek_back_multiplicity_mut()
    }
}

impl<AlignmentType: IAlignmentType + Clone> Iterator
    for CompactAlignmentIterCloned<'_, AlignmentType>
{
    type Item = (usize, AlignmentType);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter
            .next()
            .map(|(multiplicity, alignment_type)| (multiplicity, alignment_type.clone()))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<AlignmentType: IAlignmentType + Clone> DoubleEndedIterator
    for CompactAlignmentIterCloned<'_, AlignmentType>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter
            .next_back()
            .map(|(multiplicity, alignment_type)| (multiplicity, alignment_type.clone()))
    }
}

impl<AlignmentType: IAlignmentType + Clone> ExactSizeIterator
    for CompactAlignmentIterCloned<'_, AlignmentType>
{
}

//////////////////////////////////////////////
// FlatAlignmentIter /////////////////////////
//////////////////////////////////////////////

impl<'alignment, AlignmentType: IAlignmentType> FlatAlignmentIter<'alignment, AlignmentType> {
    pub fn new(alignment: &'alignment [(usize, AlignmentType)]) -> Self {
        Self {
            length: alignment
                .iter()
                .map(|(multiplicity, alignment_type)| {
                    if alignment_type.is_repeatable() {
                        *multiplicity
                    } else {
                        1
                    }
                })
                .sum(),
            iter: CompactAlignmentIter::new(alignment),
        }
    }

    /// Peek at the first alignment in the iterator.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front(&self) -> Option<&'alignment AlignmentType> {
        self.iter
            .peek_front()
            .map(|(_, alignment_type)| alignment_type)
    }

    /// Peek at the last alignment in the iterator.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back(&self) -> Option<&'alignment AlignmentType> {
        self.iter
            .peek_back()
            .map(|(_, alignment_type)| alignment_type)
    }

    /// Peek at the first alignment in the iterator and clone it.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_cloned(&self) -> Option<AlignmentType>
    where
        AlignmentType: IAlignmentType + Clone,
    {
        self.iter
            .peek_front_cloned()
            .map(|(_, alignment_type)| alignment_type)
    }

    /// Peek at the last alignment in the iterator and clone it.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_cloned(&self) -> Option<AlignmentType>
    where
        AlignmentType: IAlignmentType + Clone,
    {
        self.iter
            .peek_back_cloned()
            .map(|(_, alignment_type)| alignment_type)
    }

    /// Mutably peek at the multiplicity of the first alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<
        PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, RefMutLengthDecreaser<'iter>>,
    > {
        self.peek_front().is_some().then(|| {
            PeekFrontMultiplicityMut::new(
                &mut self.iter,
                RefMutLengthDecreaser {
                    length: &mut self.length,
                },
            )
        })
    }

    /// Mutably peek at the multiplicity of the last alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<
        PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, RefMutLengthDecreaser<'iter>>,
    > {
        self.peek_back().is_some().then(|| {
            PeekBackMultiplicityMut::new(
                &mut self.iter,
                RefMutLengthDecreaser {
                    length: &mut self.length,
                },
            )
        })
    }
}

impl<'alignment, AlignmentType: IAlignmentType> Iterator
    for FlatAlignmentIter<'alignment, AlignmentType>
{
    type Item = &'alignment AlignmentType;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(result) = self.peek_front() {
            *self.peek_front_multiplicity_mut().unwrap() -= 1;
            Some(result)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.length, Some(self.length))
    }
}

impl<AlignmentType: IAlignmentType> DoubleEndedIterator for FlatAlignmentIter<'_, AlignmentType> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some(result) = self.peek_back() {
            *self.peek_back_multiplicity_mut().unwrap() -= 1;
            Some(result)
        } else {
            None
        }
    }
}

impl<AlignmentType: IAlignmentType> ExactSizeIterator for FlatAlignmentIter<'_, AlignmentType> {}

//////////////////////////////////////////////
// FlatAlignmentIterCloned ///////////////////
//////////////////////////////////////////////

impl<'alignment, AlignmentType: IAlignmentType + Clone>
    FlatAlignmentIterCloned<'alignment, AlignmentType>
{
    pub fn new(alignment: &'alignment [(usize, AlignmentType)]) -> Self {
        Self {
            iter: FlatAlignmentIter::new(alignment),
        }
    }

    /// Peek at the first alignment in the iterator.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front(&self) -> Option<&'alignment AlignmentType> {
        self.iter.peek_front()
    }

    /// Peek at the last alignment in the iterator.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back(&self) -> Option<&'alignment AlignmentType> {
        self.iter.peek_back()
    }

    /// Peek at the first alignment in the iterator and clone it.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_cloned(&self) -> Option<AlignmentType> {
        self.iter.peek_front_cloned()
    }

    /// Peek at the last alignment in the iterator and clone it.
    ///
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_cloned(&self) -> Option<AlignmentType> {
        self.iter.peek_back_cloned()
    }

    /// Mutably peek at the multiplicity of the first alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_front_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<
        PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, RefMutLengthDecreaser<'iter>>,
    > {
        self.iter.peek_front_multiplicity_mut()
    }

    /// Mutably peek at the multiplicity of the last alignment in the iterator.
    ///
    /// It is allowed to decrease the multiplicity, which is interpreted as taking that many elements out of the iterator.
    /// If the iterator is empty, `None` is returned.
    pub fn peek_back_multiplicity_mut<'iter>(
        &'iter mut self,
    ) -> Option<
        PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, RefMutLengthDecreaser<'iter>>,
    > {
        self.iter.peek_back_multiplicity_mut()
    }
}

impl<AlignmentType: IAlignmentType + Clone> Iterator
    for FlatAlignmentIterCloned<'_, AlignmentType>
{
    type Item = AlignmentType;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().cloned()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<AlignmentType: IAlignmentType + Clone> DoubleEndedIterator
    for FlatAlignmentIterCloned<'_, AlignmentType>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().cloned()
    }
}

impl<AlignmentType: IAlignmentType + Clone> ExactSizeIterator
    for FlatAlignmentIterCloned<'_, AlignmentType>
{
}

//////////////////////////////////////////////
// Peek multiplicities ///////////////////////
//////////////////////////////////////////////

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn new(
        iter: &'iter mut CompactAlignmentIter<'alignment, AlignmentType>,
        length_decreaser: LengthDecreaser,
    ) -> Self {
        Self {
            initial_multiplicity: iter.front_multiplicity,
            iter,
            length_decreaser,
        }
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    Deref for PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    type Target = usize;

    fn deref(&self) -> &Self::Target {
        &self.iter.front_multiplicity
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    DerefMut for PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.iter.front_multiplicity
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    Drop for PeekFrontMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn drop(&mut self) {
        assert!(
            self.initial_multiplicity >= self.iter.front_multiplicity,
            "Increasing the multiplicity is forbidden while peeking mutably. Only decreasing is allowed."
        );
        self.length_decreaser
            .decrease_length(self.initial_multiplicity - self.iter.front_multiplicity);
        self.iter.advance_zeroes_front();
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn new(
        iter: &'iter mut CompactAlignmentIter<'alignment, AlignmentType>,
        length_decreaser: LengthDecreaser,
    ) -> Self {
        Self {
            initial_multiplicity: iter.back_multiplicity,
            iter,
            length_decreaser,
        }
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    Deref for PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    type Target = usize;

    fn deref(&self) -> &Self::Target {
        &self.iter.back_multiplicity
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    DerefMut for PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.iter.back_multiplicity
    }
}

impl<'iter, 'alignment: 'iter, AlignmentType: IAlignmentType, LengthDecreaser: PeekLengthDecreaser>
    Drop for PeekBackMultiplicityMut<'iter, 'alignment, AlignmentType, LengthDecreaser>
{
    fn drop(&mut self) {
        assert!(
            self.initial_multiplicity >= self.iter.back_multiplicity,
            "Increasing the multiplicity is forbidden while peeking mutably. Only decreasing is allowed."
        );
        self.length_decreaser
            .decrease_length(self.initial_multiplicity - self.iter.back_multiplicity);
        self.iter.advance_zeroes_back();
    }
}

//////////////////////////////////////////////
// LengthDecreasers //////////////////////////
//////////////////////////////////////////////

impl PeekLengthDecreaser for NoLengthDecreaser {
    fn decrease_length(&mut self, _decrement: usize) {
        // Do nothing.
    }
}

impl PeekLengthDecreaser for RefMutLengthDecreaser<'_> {
    fn decrease_length(&mut self, decrement: usize) {
        *self.length -= decrement;
    }
}
