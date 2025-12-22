use crate::graph_maps::VertGraphMap;
use crate::graphs::cube::CubeGraph;
use crate::graphs::UGraph;
use crate::prelude::*;

#[derive(Debug)]
struct CubeMap<'u, 'v, D: Dim, V: UGraph> {
    map: VertGraphMap<'u, 'v, CubeGraph<D>, V>,
}

// Inclusion / Forgetful functor
impl<'u, 'v, D: Dim, V: UGraph> From<VertGraphMap<'u, 'v, CubeGraph<D>, V>>
    for CubeMap<'u, 'v, D, V>
{
    fn from(value: VertGraphMap<'u, 'v, CubeGraph<D>, V>) -> Self {
        Self { map: value }
    }
}

impl<'u, 'v, D: Dim, V: UGraph> From<CubeMap<'u, 'v, D, V>>
    for VertGraphMap<'u, 'v, CubeGraph<D>, V>
{
    fn from(value: CubeMap<'u, 'v, D, V>) -> Self {
        value.map
    }
}
