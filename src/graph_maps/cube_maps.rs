use crate::graph_maps::{GraphMapError, VertGraphMap};
use crate::graphs::cube::CubeGraph;
use crate::graphs::UGraph;
use crate::prelude::*;

#[derive(Debug)]
struct CubeMap<'u, 'v, D: Dim, V: UGraph> {
    map: VertGraphMap<'u, 'v, CubeGraph<D>, V>,
    // degenerate_indices: Vec<bool>,
}

// Inclusion / Forgetful functor
impl<'u, 'v, D: Dim, V: UGraph> From<VertGraphMap<'u, 'v, CubeGraph<D>, V>>
    for CubeMap<'u, 'v, D, V>
{
    fn from(value: VertGraphMap<'u, 'v, CubeGraph<D>, V>) -> Self {
        // let dim = value.domain.dim().size() as usize;
        Self {
            map: value,
            // degenerate_indices: vec![false; dim],
        }
    }
}

impl<'u, 'v, D: Dim, V: UGraph> From<CubeMap<'u, 'v, D, V>>
    for VertGraphMap<'u, 'v, CubeGraph<D>, V>
{
    fn from(value: CubeMap<'u, 'v, D, V>) -> Self {
        value.map
    }
}

impl<'u, 'v, V: UGraph> CubeMap<'u, 'v, u32, V> {
    pub fn dim(&self) -> u32 {
        self.map.domain.dim()
    }

    pub fn try_combine<'w>(
        &self,
        other: &CubeMap<'w, 'v, u32, V>,
        combined: &'w CubeGraph<u32>,
    ) -> Result<CubeMap<'w, 'v, u32, V>, GraphMapError> {
        panic!("Not implemented")
    }
}
