use std::borrow::Cow;

use crate::graph_maps::{GraphMap, GraphMapError, VertGraphMap};
use crate::graphs::cube::CubeGraph;
use crate::graphs::UGraph;
use crate::prelude::*;

#[derive(Debug)]
pub struct CubeMap<'u, 'v, D: Dim, V: UGraph> {
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
    ) -> Result<CubeMap<'w, 'v, u32, V>, GraphMapError>
    where
        'u: 'w,
    {
        use GraphMapError as E;
        assert!(self.map.codomain == other.map.codomain);
        let n = self.map.domain.n();
        for (i, (x, y)) in self
            .map
            .vert_maps
            .iter()
            .zip(other.map.vert_maps.iter())
            .enumerate()
        {
            if !self.map.codomain.is_edge(*x, *y) {
                return Err(E::BadEdge(i as u32, i as u32 + n, *x, *y));
            }
        }

        let mut combined_verts =
            Vec::with_capacity(self.map.vert_maps.len() + other.map.vert_maps.len());

        combined_verts.extend_from_slice(&self.map.vert_maps);
        combined_verts.extend_from_slice(&other.map.vert_maps);
        let combined = Cow::Owned(CubeGraph::new(n + 1));

        let map = unsafe {
            VertGraphMap::new_unchecked(
                combined,
                // Codomain is almost certainly borrowed so this shouldn't be an issue.:where
                self.map.codomain.clone(),
                Cow::Owned(combined_verts),
            )
        };
        Ok(Self { map })
    }

    pub fn d(i: u32, sign: bool) -> CubeMap<'u, 'v, u32, V> {
        panic!()
    }
}
