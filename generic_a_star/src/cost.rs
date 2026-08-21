use std::{
    fmt::{Debug, Display},
    hash::Hash,
    iter::Sum,
    ops::{Add, AddAssign, Mul, Sub, SubAssign},
    str::FromStr,
};

use num_traits::{Bounded, CheckedAdd, CheckedSub, SaturatingAdd, SaturatingSub, Zero};

/// The cost of an A* node.
pub trait AStarCost:
    From<u8>
    + Add<Output = Self>
    + Sum
    + Sub<Output = Self>
    + Mul<Output = Self>
    + SaturatingAdd
    + SaturatingSub
    + CheckedAdd
    + CheckedSub
    + AddAssign
    + SubAssign
    + FromStr
    + Zero
    + Bounded
    + Display
    + Debug
    + Ord
    + Eq
    + Hash
    + Copy
{
    type CostType;

    fn as_f64(&self) -> f64;

    fn as_u64(&self) -> u64;

    fn as_primitive(&self) -> Self::CostType;

    fn from_primitive(value: Self::CostType) -> Self;

    fn as_usize(&self) -> usize;

    fn from_usize(value: usize) -> Self;
}

macro_rules! primitive_cost {
    ($name:ident, $primitive:ident) => {
        #[doc = concat!("The cost of an A* node.\n\nThis cost type uses [`", stringify!($primitive), "`] as the internal representation of cost.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name($primitive);

        impl $crate::cost::AStarCost for $name {
            type CostType = $primitive;

            fn as_f64(&self) -> f64 {
                self.0 as f64
            }

            fn as_u64(&self)->u64{
                self.0 as u64
            }

            fn as_primitive(&self) -> Self::CostType {
                self.0
            }

            fn from_primitive(value: Self::CostType) -> Self {
                Self(value)
            }

            fn as_usize(&self) -> usize {
                self.as_primitive().try_into().unwrap()
            }

            fn from_usize(value: usize) -> Self {
                Self::from_primitive(value.try_into().unwrap())
            }
        }

        impl From<$primitive> for $name {
            fn from(value: $primitive) -> Self {
                Self(value)
            }
        }

        impl From<u8> for $name {
            fn from(value: u8) -> Self {
                Self(value.into())
            }
        }

        impl std::ops::Add for $name {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                Self(self.0.checked_add(rhs.0).unwrap_or_else(|| panic!("Overflow when adding costs {self} + {rhs}")))
            }
        }

        impl std::iter::Sum for $name {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                iter.fold(Self::zero(), |acc, x| acc + x)
            }
        }

        impl std::ops::Sub for $name {
            type Output = Self;

            fn sub(self, rhs: Self) -> Self::Output {
                Self(self.0.checked_sub(rhs.0).unwrap_or_else(|| panic!("Overflow when subtracting costs {self} - {rhs}")))
            }
        }

        impl std::ops::Mul for $name {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                Self(self.0.checked_mul(rhs.0).unwrap_or_else(|| panic!("Overflow when multiplying costs {self} * {rhs}")))
            }
        }

        impl num_traits::SaturatingAdd for $name {
            fn saturating_add(&self, rhs: &Self) -> Self {
                Self(self.0.saturating_add(rhs.0))
            }
        }

        impl num_traits::SaturatingSub for $name {
            fn saturating_sub(&self, rhs: &Self) -> Self {
                Self(self.0.saturating_sub(rhs.0))
            }
        }

        impl num_traits::CheckedAdd for $name {
            fn checked_add(&self, rhs: &Self) -> Option<Self> {
                self.0.checked_add(rhs.0).map($name)
            }
        }

        impl num_traits::CheckedSub for $name {
            fn checked_sub(&self, rhs: &Self) -> Option<Self> {
                self.0.checked_sub(rhs.0).map($name)
            }
        }

        impl std::ops::AddAssign for $name {
            fn add_assign(&mut self, rhs: Self) {
                *self = *self + rhs;
            }
        }

        impl std::ops::SubAssign for $name {
            fn sub_assign(&mut self, rhs: Self) {
                *self = *self - rhs;
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                std::fmt::Display::fmt(&self.0, f)
            }
        }

        impl std::str::FromStr for $name {
            type Err = <$primitive as std::str::FromStr>::Err;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                $primitive::from_str(s).map(Self)
            }
        }

        impl num_traits::Zero for $name {
            fn zero() -> Self {
                $primitive::zero().into()
            }

            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }

        impl num_traits::Bounded for $name {
            fn min_value() -> Self {
                <$primitive as num_traits::Bounded>::min_value().into()
            }

            fn max_value() -> Self {
                <$primitive as num_traits::Bounded>::max_value().into()
            }
        }
    };
}

primitive_cost!(U16Cost, u16);
primitive_cost!(I16Cost, i16);
primitive_cost!(U32Cost, u32);
primitive_cost!(I32Cost, i32);
primitive_cost!(U64Cost, u64);
primitive_cost!(I64Cost, i64);

/// A pair of cost types where the first cost is prioritised over the second cost.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OrderedPairCost<A, B>(pub A, pub B);

impl<A: AStarCost, B: AStarCost> AStarCost for OrderedPairCost<A, B> {
    type CostType = A::CostType;

    fn as_f64(&self) -> f64 {
        self.0.as_f64()
    }

    fn as_u64(&self) -> u64 {
        self.0.as_u64()
    }

    fn as_primitive(&self) -> Self::CostType {
        self.0.as_primitive()
    }

    fn from_primitive(value: Self::CostType) -> Self {
        Self(A::from_primitive(value), B::zero())
    }

    fn as_usize(&self) -> usize {
        self.0.as_usize()
    }

    fn from_usize(value: usize) -> Self {
        Self(A::from_usize(value), B::zero())
    }
}

impl<A: Display, B: Display> Display for OrderedPairCost<A, B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

impl<A: Zero, B: Zero> Zero for OrderedPairCost<A, B> {
    fn zero() -> Self {
        Self(A::zero(), B::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero()
    }
}

impl<A: Bounded, B: Bounded> Bounded for OrderedPairCost<A, B> {
    fn min_value() -> Self {
        Self(A::min_value(), B::min_value())
    }

    fn max_value() -> Self {
        Self(A::max_value(), B::max_value())
    }
}

impl<A: From<u8>, B: From<u8>> From<u8> for OrderedPairCost<A, B> {
    fn from(value: u8) -> Self {
        Self(A::from(value), B::from(0))
    }
}

impl<A: Add<Output = A>, B: Add<Output = B>> Add for OrderedPairCost<A, B> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl<A, B> Sum for OrderedPairCost<A, B>
where
    OrderedPairCost<A, B>: Add<Output = OrderedPairCost<A, B>> + Zero,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<A: Sub<Output = A>, B: Sub<Output = B>> Sub for OrderedPairCost<A, B> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl<A: Mul<Output = A>, B: Mul<Output = B>> Mul for OrderedPairCost<A, B> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0, self.1 * rhs.1)
    }
}

impl<A: FromStr, B: FromStr> FromStr for OrderedPairCost<A, B> {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let Some(s) = s.strip_prefix('(') else {
            return Err(());
        };
        let Some(s) = s.strip_suffix(')') else {
            return Err(());
        };
        let mut parts = s.splitn(2, ',');
        let s = parts.next().ok_or(())?.trim();
        let a = s.parse::<A>().map_err(|_| ())?;
        let s = parts.next().ok_or(())?.trim();
        let b = s.parse::<B>().map_err(|_| ())?;
        Ok(Self(a, b))
    }
}

impl<A: AddAssign, B: AddAssign> AddAssign for OrderedPairCost<A, B> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
        self.1 += rhs.1;
    }
}

impl<A: SubAssign, B: SubAssign> SubAssign for OrderedPairCost<A, B> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
        self.1 -= rhs.1;
    }
}

impl<A: CheckedAdd, B: CheckedAdd> CheckedAdd for OrderedPairCost<A, B> {
    fn checked_add(&self, rhs: &Self) -> Option<Self> {
        Some(Self(
            self.0.checked_add(&rhs.0)?,
            self.1.checked_add(&rhs.1)?,
        ))
    }
}

impl<A: CheckedSub, B: CheckedSub> CheckedSub for OrderedPairCost<A, B> {
    fn checked_sub(&self, rhs: &Self) -> Option<Self> {
        Some(Self(
            self.0.checked_sub(&rhs.0)?,
            self.1.checked_sub(&rhs.1)?,
        ))
    }
}

impl<A: SaturatingAdd, B: SaturatingAdd> SaturatingAdd for OrderedPairCost<A, B> {
    fn saturating_add(&self, rhs: &Self) -> Self {
        Self(self.0.saturating_add(&rhs.0), self.1.saturating_add(&rhs.1))
    }
}

impl<A: SaturatingSub, B: SaturatingSub> SaturatingSub for OrderedPairCost<A, B> {
    fn saturating_sub(&self, rhs: &Self) -> Self {
        Self(self.0.saturating_sub(&rhs.0), self.1.saturating_sub(&rhs.1))
    }
}
