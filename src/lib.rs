// #![feature(generic_const_exprs)]
// #![allow(incomplete_features)]

pub mod computations;
pub mod graph_maps;
pub mod graphs;
pub mod shape;

pub mod prelude {
    pub use super::graph_maps::uint::UINT;
    pub use super::graphs::cube::Newable;
    pub use super::graphs::*;
    pub use super::shape::*;
    pub use std::convert::*;
}
