pub mod graph_maps;
pub mod graphs;
pub mod shape;

pub mod prelude {
    pub use crate::graphs::*;
    pub use crate::shape::*;
    pub use std::convert::*;
}

fn main() {
    println!("Hello, world!");
}
