pub mod graph_maps;
pub mod graphs;
pub mod shape;

pub mod prelude {
    pub use super::graphs::cube::Newable;
    pub use super::graphs::*;
    pub use super::shape::*;
    pub use std::convert::*;
}
