/// https://github.com/coreylowman/dfdx/blob/5e0c3ddda5de3258d52f9d728426adc191ae1f5e/dfdx-core/src/shapes/shape.rs#L9

/// Represents a single dimension of a multi dimensional [Shape]
pub trait Dim: 'static + Copy + Clone + std::fmt::Debug + Send + Sync + Eq + PartialEq {
    fn size(&self) -> usize;
    fn from_size(size: usize) -> Option<Self>;
}

/// Represents a single dimension where all
/// instances are guaranteed to be the same size at compile time.
pub trait ConstDim: Default + Dim {
    const SIZE: usize;
}

impl Dim for usize {
    #[inline(always)]
    fn size(&self) -> usize {
        *self
    }
    #[inline(always)]
    fn from_size(size: usize) -> Option<Self> {
        Some(size)
    }
}

/// Represents a [Dim] with size known at compile time
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Const<const M: usize>;
impl<const M: usize> Dim for Const<M> {
    #[inline(always)]
    fn size(&self) -> usize {
        M
    }
    #[inline(always)]
    fn from_size(size: usize) -> Option<Self> {
        if size == M {
            Some(Const)
        } else {
            None
        }
    }
}

impl<const M: usize> ConstDim for Const<M> {
    const SIZE: usize = M;
}

impl<const N: usize> core::ops::Add<Const<N>> for usize {
    type Output = usize;
    fn add(self, _: Const<N>) -> Self::Output {
        self.size() + N
    }
}
impl<const N: usize> core::ops::Add<usize> for Const<N> {
    type Output = usize;
    fn add(self, rhs: usize) -> Self::Output {
        N + rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: usize, const M: usize> core::ops::Add<Const<N>> for Const<M>
where
    Const<{ M + N }>: Sized,
{
    type Output = Const<{ M + N }>;
    fn add(self, _: Const<N>) -> Self::Output {
        Const
    }
}

impl<const N: usize> core::ops::Mul<Const<N>> for usize {
    type Output = usize;
    fn mul(self, _: Const<N>) -> Self::Output {
        self.size() * N
    }
}
impl<const N: usize> core::ops::Mul<usize> for Const<N> {
    type Output = usize;
    fn mul(self, rhs: usize) -> Self::Output {
        N * rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: usize, const M: usize> core::ops::Mul<Const<N>> for Const<M>
where
    Const<{ M * N }>: Sized,
{
    type Output = Const<{ M * N }>;
    fn mul(self, _: Const<N>) -> Self::Output {
        Const
    }
}

impl<const N: usize> core::ops::Div<Const<N>> for usize {
    type Output = usize;
    fn div(self, _: Const<N>) -> Self::Output {
        self.size() / N
    }
}
impl<const N: usize> core::ops::Div<usize> for Const<N> {
    type Output = usize;
    fn div(self, rhs: usize) -> Self::Output {
        N / rhs.size()
    }
}

#[cfg(feature = "nightly")]
impl<const N: usize, const M: usize> core::ops::Div<Const<N>> for Const<M>
where
    Const<{ M / N }>: Sized,
{
    type Output = Const<{ M / N }>;
    fn div(self, _: Const<N>) -> Self::Output {
        Const
    }
}

/// Represents either `[T; N]` or `Vec<T>`
pub trait Array<T>: IntoIterator<Item = T> {
    type Dim: Dim;
    fn dim(&self) -> Self::Dim;
}
impl<T, const N: usize> Array<T> for [T; N] {
    type Dim = Const<N>;
    fn dim(&self) -> Self::Dim {
        Const
    }
}
impl<T> Array<T> for std::vec::Vec<T> {
    type Dim = usize;
    fn dim(&self) -> Self::Dim {
        self.len()
    }
}
