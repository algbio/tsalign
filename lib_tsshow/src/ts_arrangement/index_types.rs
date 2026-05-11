use strong_type::StrongType;

/// A column index in a source string.
#[derive(StrongType)]
#[strong_type(conversion, addable)]
pub struct SourceColumn(usize);

/// A column index over the characters in an arrangement.
///
/// This includes copies and hidden characters, but not gaps and blanks.
#[derive(StrongType)]
#[strong_type(conversion, addable)]
pub struct ArrangementCharColumn(usize);

/// A column index in an arrangement.
///
/// This is the index that represents where the character will actually be rendered.
#[derive(StrongType)]
#[strong_type(conversion, addable)]
pub struct ArrangementColumn(usize);

/// Identifier of a template switch inner sequence.
#[derive(StrongType)]
#[strong_type(conversion, addable)]
pub struct TsInnerIdentifier(usize);

macro_rules! index_type_ops {
    ($name:ty, $inner:ty, $signed_inner:ty) => {
        impl $name {
            pub fn saturating_add(&self, other: $inner) -> Self {
                Self(self.0.saturating_add(other))
            }

            pub fn checked_sub(&self, other: $inner) -> Option<Self> {
                self.0.checked_sub(other).map(Self)
            }

            pub fn saturating_sub(&self, other: $inner) -> Self {
                Self(self.0.saturating_sub(other))
            }
        }

        impl std::ops::Add<$inner> for $name {
            type Output = Self;

            fn add(self, other: $inner) -> Self::Output {
                Self(self.0 + other)
            }
        }

        impl std::ops::Sub<$inner> for $name {
            type Output = Self;

            fn sub(self, other: $inner) -> Self::Output {
                Self(self.0 - other)
            }
        }

        impl std::ops::Add<$name> for $inner {
            type Output = $name;

            fn add(self, other: $name) -> Self::Output {
                (self + other.0).into()
            }
        }

        impl std::ops::Sub<$name> for $inner {
            type Output = $name;

            fn sub(self, other: $name) -> Self::Output {
                (self - other.0).into()
            }
        }

        impl std::ops::Add<$signed_inner> for $name {
            type Output = Self;

            fn add(self, other: $signed_inner) -> Self::Output {
                Self(
                    <$inner>::try_from(<$signed_inner>::try_from(self.0).unwrap() + other).unwrap(),
                )
            }
        }

        impl std::ops::Sub<$signed_inner> for $name {
            type Output = Self;

            fn sub(self, other: $signed_inner) -> Self::Output {
                Self(
                    <$inner>::try_from(<$signed_inner>::try_from(self.0).unwrap() - other).unwrap(),
                )
            }
        }

        impl std::ops::AddAssign<$inner> for $name {
            fn add_assign(&mut self, other: $inner) {
                self.0 += other
            }
        }

        impl std::ops::SubAssign<$inner> for $name {
            fn sub_assign(&mut self, other: $inner) {
                self.0 -= other
            }
        }

        impl<'a> From<&'a $name> for $inner {
            fn from(value: &'a $name) -> Self {
                value.0
            }
        }
    };
}

index_type_ops!(SourceColumn, usize, isize);
index_type_ops!(ArrangementCharColumn, usize, isize);
index_type_ops!(ArrangementColumn, usize, isize);
index_type_ops!(TsInnerIdentifier, usize, isize);
