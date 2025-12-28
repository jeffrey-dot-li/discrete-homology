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
    ) -> Result<(CubeMap<'w, 'v, u32, V>, CubeMap<'w, 'v, u32, V>), GraphMapError>
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
        let combined_len = self.map.vert_maps.len() + other.map.vert_maps.len();

        let mut combined_verts = Vec::with_capacity(combined_len);
        let mut combined_verts_2 = Vec::with_capacity(combined_len);

        combined_verts.extend_from_slice(&self.map.vert_maps);
        combined_verts.extend_from_slice(&other.map.vert_maps);
        combined_verts_2.extend_from_slice(&other.map.vert_maps);
        combined_verts_2.extend_from_slice(&self.map.vert_maps);

        let combined = Cow::Owned(CubeGraph::new(n + 1));
        let combined_2 = Cow::Owned(CubeGraph::new(n + 1));
        let map = unsafe {
            VertGraphMap::new_unchecked(
                combined,
                // Codomain is almost certainly borrowed so this shouldn't be an issue.:where
                self.map.codomain.clone(),
                Cow::Owned(combined_verts),
            )
        };
        let map2 = unsafe {
            VertGraphMap::new_unchecked(
                combined_2,
                // Codomain is almost certainly borrowed so this shouldn't be an issue.:where
                self.map.codomain.clone(),
                Cow::Owned(combined_verts_2),
            )
        };
        Ok((Self { map }, Self { map: map2 }))
    }

    pub fn d(i: u32, sign: bool) -> CubeMap<'u, 'v, u32, V> {
        panic!()
    }
}

pub fn combined_cube_maps<'u, 'v, V: UGraph>(
    maps: &[CubeMap<'u, 'v, u32, V>],
) -> Vec<CubeMap<'u, 'v, u32, V>> {
    let mut combined_maps = Vec::new();
    let len = maps.len();
    for i in 0..len {
        // if let Ok(combined) = maps[i].try_combine(&maps[i]) {
        //     combined_maps.push(combined.0);
        //     // combined_maps.push(combined.1);
        // } else {
        //     panic!("Could not combine map with itself");
        // }

        for j in (0..len) {
            if let Ok(combined) = maps[i].try_combine(&maps[j]) {
                combined_maps.push(combined.0);
                // combined_maps.push(combined.1);
            }
        }
    }
    combined_maps
}

#[cfg(test)]
mod tests {
    use crate::graph_maps::generate_maps_naive;

    use super::*;

    #[test]
    fn test_2cube_gsphere_combined() {
        let n = 2;

        use cube::CubeGraph;
        let source = CubeGraph::new(n);
        let cube_prev = CubeGraph::new(n - 1);
        // let target = extras::greene_sphere();
        let target = extras::C_N_graph(5);
        let (cube_prev_maps, cube_prev_numchecked) = generate_maps_naive(&cube_prev, &target);
        let cube_prev_maps = cube_prev_maps
            .into_iter()
            .map(|m| CubeMap::from(m))
            .collect::<Vec<_>>();
        // black_box prevents the compiler from optimizing away inputs/outputs
        let cube_n_combined_maps = combined_cube_maps(&cube_prev_maps);

        let cube_n_naive_maps = generate_maps_naive(&source, &target).0;

        assert!(
            cube_n_combined_maps.len() == cube_n_naive_maps.len(),
            "num maps combined was {}, but naive was {}",
            cube_n_combined_maps.len(),
            cube_n_naive_maps.len()
        );
    }
}
