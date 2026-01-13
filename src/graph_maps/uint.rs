use num_traits::{PrimInt, Unsigned};
use std::fmt::Debug;
use std::hash::Hash;
pub trait UINT: Unsigned + PrimInt + Default + Debug + Hash {}
impl<T> UINT for T where T: Unsigned + PrimInt + Default + Debug + Hash {}
