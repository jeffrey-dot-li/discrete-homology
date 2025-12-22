use std::convert::TryInto;

// https://github.com/coreylowman/dfdx/blob/5e0c3ddda5de3258d52f9d728426adc191ae1f5e/dfdx-core/src/shapes/shape.rs#L9

/// Represents a single dimension of a multi dimensional [Shape]
pub trait Dim<T: TryInto<usize> + Copy = u32>:
    'static + Copy + Clone + std::fmt::Debug + Send + Sync + Eq + PartialEq
{
    fn size(&self) -> T;
    fn from_size(size: T) -> Option<Self>;
}

/// Represents a single dimension where all
/// instances are guaranteed to be the same size at compile time.
pub trait ConstDim<T: TryInto<usize> + Copy = u32>: Default + Dim {
    const SIZE: T;
}

impl Dim<u32> for u32 {
    #[inline(always)]
    fn size(&self) -> u32 {
        *self
    }
    #[inline(always)]
    fn from_size(size: u32) -> Option<Self> {
        Some(size)
    }
}

/// Represents a [Dim] with size known at compile time
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Const<const M: u32>;
impl<const M: u32> Dim<u32> for Const<M> {
    #[inline(always)]
    fn size(&self) -> u32 {
        M
    }
    #[inline(always)]
    fn from_size(size: u32) -> Option<Self> {
        if size == M {
            Some(Const)
        } else {
            None
        }
    }
}

impl<const M: u32> ConstDim for Const<M> {
    const SIZE: u32 = M;
}

impl<const N: u32> core::ops::Add<Const<N>> for u32 {
    type Output = u32;
    fn add(self, _: Const<N>) -> Self::Output {
        self.size() + N
    }
}
impl<const N: u32> core::ops::Add<u32> for Const<N> {
    type Output = u32;
    fn add(self, rhs: u32) -> Self::Output {
        N + rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: u32, const M: u32> core::ops::Add<Const<N>> for Const<M>
where
    Const<{ M + N }>: Sized,
{
    type Output = Const<{ M + N }>;
    fn add(self, _: Const<N>) -> Self::Output {
        Const
    }
}

impl<const N: u32> core::ops::Mul<Const<N>> for u32 {
    type Output = u32;
    fn mul(self, _: Const<N>) -> Self::Output {
        self.size() * N
    }
}
impl<const N: u32> core::ops::Mul<u32> for Const<N> {
    type Output = u32;
    fn mul(self, rhs: u32) -> Self::Output {
        N * rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: u32, const M: u32> core::ops::Mul<Const<N>> for Const<M>
where
    Const<{ M * N }>: Sized,
{
    type Output = Const<{ M * N }>;
    fn mul(self, _: Const<N>) -> Self::Output {
        Const
    }
}

impl<const N: u32> core::ops::Div<Const<N>> for u32 {
    type Output = u32;
    fn div(self, _: Const<N>) -> Self::Output {
        self.size() / N
    }
}
impl<const N: u32> core::ops::Div<u32> for Const<N> {
    type Output = u32;
    fn div(self, rhs: u32) -> Self::Output {
        N / rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: u32, const M: u32> core::ops::Div<Const<N>> for Const<M>
where
    Const<{ M / N }>: Sized,
{
    type Output = Const<{ M / N }>;
    fn div(self, _: Const<N>) -> Self::Output {
        Const
    }
}
