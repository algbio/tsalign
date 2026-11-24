use std::{
    io::{Read, Write},
    mem,
    ops::{Index, IndexMut},
    slice,
};

use log::trace;

pub struct LowerBoundCostArray<const DIMENSION: usize, Cost> {
    dim: [usize; DIMENSION],
    data: Vec<Cost>,
}

impl<const DIMENSION: usize, Cost> LowerBoundCostArray<DIMENSION, Cost> {
    /// The size of the array in each dimension.
    pub fn dim(&self) -> [usize; DIMENSION] {
        self.dim
    }

    pub fn new_from_cost(dim: [usize; DIMENSION], cost: Cost) -> Self
    where
        Cost: Clone,
    {
        Self {
            dim,
            data: vec![cost; dim.into_iter().product()],
        }
    }

    pub fn write(&self, mut write: impl Write) -> std::io::Result<()>
    where
        Cost: Copy,
    {
        for dimension in self.dim {
            write.write_all(&dimension.to_ne_bytes())?;
        }

        let cost_size = mem::size_of::<Cost>();
        let data: &[u8] = unsafe {
            slice::from_raw_parts(self.data.as_ptr() as *const u8, self.data.len() * cost_size)
        };
        write.write_all(data)
    }

    pub fn read(mut read: impl Read) -> std::io::Result<Self>
    where
        Cost: Copy,
    {
        let mut buffer = [0; mem::size_of::<usize>()];
        let mut dim = [0usize; DIMENSION];
        for dimension in &mut dim {
            read.read_exact(&mut buffer)?;
            *dimension = usize::from_ne_bytes(buffer);
        }
        let dim = dim;
        trace!("Read dimensions: {dim:?}");

        let cost_size = mem::size_of::<Cost>();
        let data_len = dim.into_iter().product();
        let data_bytes_len = cost_size * data_len;

        let mut data = Vec::<Cost>::with_capacity(data_bytes_len);
        let mut data_bytes = unsafe {
            Vec::from_raw_parts(data.as_mut_ptr() as *mut u8, 0, data.capacity() * cost_size)
        };
        read.by_ref()
            .take(data_bytes_len.try_into().unwrap())
            .read_to_end(&mut data_bytes)?;
        unsafe {
            data.set_len(data_len);
        };
        data_bytes.leak();

        Ok(Self { dim, data })
    }
}

impl<const DIMENSION: usize, Cost> Index<[usize; DIMENSION]>
    for LowerBoundCostArray<DIMENSION, Cost>
{
    type Output = Cost;

    fn index(&self, index: [usize; DIMENSION]) -> &Self::Output {
        &self.data[coordinates_to_index(index, self.dim)]
    }
}

impl<const DIMENSION: usize, Cost> IndexMut<[usize; DIMENSION]>
    for LowerBoundCostArray<DIMENSION, Cost>
{
    fn index_mut(&mut self, index: [usize; DIMENSION]) -> &mut Self::Output {
        &mut self.data[coordinates_to_index(index, self.dim)]
    }
}

impl<Cost> FromIterator<Cost> for LowerBoundCostArray<1, Cost> {
    fn from_iter<T: IntoIterator<Item = Cost>>(iter: T) -> Self {
        let data: Vec<_> = iter.into_iter().collect();
        let dim = [data.len()];
        Self { dim, data }
    }
}

#[inline]
fn coordinates_to_index<const DIMENSION: usize>(
    coordinates: [usize; DIMENSION],
    dim: [usize; DIMENSION],
) -> usize {
    let mut result = 0;
    for (index, ordinate) in coordinates.into_iter().enumerate() {
        result += dim.into_iter().skip(index + 1).product::<usize>() * ordinate;
    }
    result
}
