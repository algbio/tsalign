use std::{
    fmt::Debug,
    mem,
    ops::{BitOrAssign, ShlAssign, ShrAssign},
};

use num_traits::{PrimInt, bounds::LowerBounded};

pub trait KmerStore:
    PrimInt + ShlAssign<i32> + ShrAssign<i32> + BitOrAssign + From<u8> + Into<u128>
{
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct Kmer<Store> {
    kmer: Store,
}

pub type Kmer64 = Kmer<u128>;
pub type Kmer32 = Kmer<u64>;
pub type Kmer16 = Kmer<u32>;
pub type Kmer8 = Kmer<u16>;

impl KmerStore for u128 {}
impl KmerStore for u64 {}
impl KmerStore for u32 {}
impl KmerStore for u16 {}

fn char_to_bits<Store: From<u8>>(char: u8) -> Store {
    match char {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        char => panic!("Unsupported character: {}", char as char),
    }
    .into()
}

fn bits_to_char<Store: Into<u128>>(bits: Store) -> u8 {
    match bits.into() & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!(),
    }
}

impl<Store: KmerStore> From<&[u8]> for Kmer<Store> {
    fn from(value: &[u8]) -> Self {
        assert!(value.len() <= mem::size_of::<Store>() * 4);

        let mut store = Store::zero();
        for c in value {
            store <<= 2;
            store |= char_to_bits(*c);
        }
        Self { kmer: store }
    }
}

impl<Store: KmerStore> FromIterator<u8> for Kmer<Store> {
    fn from_iter<T: IntoIterator<Item = u8>>(iter: T) -> Self {
        let mut store = Store::zero();
        for (i, c) in iter.into_iter().enumerate() {
            assert!(i < mem::size_of::<Store>() * 4);
            store <<= 2;
            store |= char_to_bits(c);
        }
        Self { kmer: store }
    }
}

impl<Store: KmerStore> LowerBounded for Kmer<Store> {
    fn min_value() -> Self {
        Self {
            kmer: Store::min_value(),
        }
    }
}

impl<Store: KmerStore> Kmer<Store> {
    pub fn push(&mut self, char: u8) {
        self.kmer <<= 2;
        self.kmer |= char_to_bits(char);
    }

    pub fn to_vec(self, k: usize) -> Vec<u8> {
        assert!(k <= mem::size_of::<Store>() * 4);
        let mut result = Vec::new();
        let mut kmer = self.kmer;
        for _ in 0..k {
            result.push(bits_to_char(kmer));
            kmer >>= 2;
        }
        result.reverse();
        result
    }

    pub fn to_string(self, k: usize) -> String {
        String::from_utf8_lossy(&self.to_vec(k)).to_string()
    }
}

impl<Store: KmerStore> Default for Kmer<Store> {
    fn default() -> Self {
        Self {
            kmer: Store::zero(),
        }
    }
}

impl<Store: KmerStore> Debug for Kmer<Store> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string(mem::size_of::<Store>() * 4))
    }
}
