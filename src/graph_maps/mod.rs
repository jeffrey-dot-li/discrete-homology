pub mod cube_maps;
use crate::prelude::*;
use std::fmt::Debug;

pub trait GraphMap<'a, U, V>
where
    U: UGraph,
    V: UGraph,
{
    fn domain(&self) -> &'a U;
    fn codomain(&self) -> &'a V;
    fn map(&self, u: u32) -> u32;
}
#[derive(Debug)]
pub struct VertGraphMap<'a, U, V>
where
    U: UGraph,
    V: UGraph,
{
    domain: &'a U,
    codomain: &'a V,
    vert_maps: Vec<u32>,
}
#[derive(Debug)]
pub enum GraphMapError {
    BadEdge(u32, u32, u32, u32),
    InvalidMap(String),
}

impl<'a, U, V> VertGraphMap<'a, U, V>
where
    U: UGraph,
    V: UGraph,
{
    pub fn try_new(
        domain: &'a U,
        codomain: &'a V,
        vert_maps: Vec<u32>,
    ) -> Result<Self, GraphMapError> {
        use GraphMapError as E;
        if vert_maps.len() != domain.n() as usize {
            return Err(E::InvalidMap(format!(
                "invalid map: len {:?} != {:?}",
                vert_maps.len(),
                domain.n()
            )));
        }

        // This has O(domain.edges) runtime

        for i in 0..domain.n() {
            let mapped = vert_maps[i as usize];

            if mapped >= codomain.n() {
                return Err(E::InvalidMap(format!(
                    "invalid map: codomain node {mapped} out of range"
                )));
            }
            let neighbors = domain.neighbors(i);
            for neighbor in neighbors {
                // We only need to check for neighbor < i
                // neighbor == i is guarenteed because codomain is a valid UGraph
                // We guarentee checking for neighbor > i because is_edge(i, neighbor) == is_edge(neighbor, i)
                // This way we don't need to check out of bounds twice
                if neighbor < i {
                    let mapped_neighbor = vert_maps[neighbor as usize];
                    if !codomain.is_edge(mapped_neighbor, mapped) {
                        return Err(E::BadEdge(i, neighbor, mapped, mapped_neighbor));
                    }
                }
            }
        }

        Ok(Self {
            domain,
            codomain,
            vert_maps,
        })
    }
}

impl<'a, U, V> GraphMap<'a, U, V> for VertGraphMap<'a, U, V>
where
    U: UGraph,
    V: UGraph,
{
    fn domain(&self) -> &'a U {
        self.domain
    }
    fn codomain(&self) -> &'a V {
        self.codomain
    }
    fn map(&self, u: u32) -> u32 {
        self.vert_maps[u as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphs::cube::CubeGraph;
    use crate::graphs::extras::greene_sphere;

    #[test]
    fn test_graph_map() {
        let cube = CubeGraph::<Const<2>>::default();
        let gsphere_graph = greene_sphere();

        let id_2cube = VertGraphMap::try_new(&cube, &cube, vec![0, 1, 2, 3]);
        assert!(id_2cube.is_ok(), "{:?}", id_2cube);

        let id_2gsphere = VertGraphMap::try_new(
            &gsphere_graph,
            &gsphere_graph,
            (0..gsphere_graph.n()).collect(),
        );
        assert!(id_2gsphere.is_ok());

        let neg_gsphere = VertGraphMap::try_new(
            &gsphere_graph,
            &gsphere_graph,
            (0..gsphere_graph.n()).rev().collect(),
        );
        assert!(neg_gsphere.is_ok());
    }

    #[test]
    fn test_graph_map_bad_edge() {
        let cube = CubeGraph::<Const<2>>::default();

        // 2-cube vertices and their neighbors:
        // 0: [0,1,2] (connected to 1 and 2, differs by one bit)
        // 1: [0,1,3]
        // 2: [0,2,3]
        // 3: [1,2,3]
        //
        // Map [0,3,1,2]: This should fail because:
        // - Vertex 0 has edge to vertex 1 in domain
        // - But f(0)=0 and f(1)=3, and (0,3) is NOT an edge in 2-cube
        //   (0=00b and 3=11b differ by 2 bits, not 1)
        let mapping = vec![0, 3, 1, 2];
        let invalid_map = VertGraphMap::try_new(&cube, &cube, mapping.clone());

        match invalid_map {
            Err(GraphMapError::BadEdge(_v1, _v2, _m1, _m2)) => {
                // Expected error
            }
            Err(GraphMapError::InvalidMap(msg)) => {
                panic!("Expected BadEdge, got InvalidMap: {msg}");
            }
            Ok(_) => {
                panic!(
                    "Expected BadEdge error, but GraphMap was created successfully.\n\
                     Mapping: {mapping:?}"
                );
            }
        }
    }

    #[test]
    fn test_graph_map_invalid_length() {
        let cube = CubeGraph::<Const<2>>::default();

        // 2-cube has 4 vertices, but we provide only 3 mappings
        let mapping = vec![0, 1, 2];
        let invalid_map = VertGraphMap::try_new(&cube, &cube, mapping.clone());

        match invalid_map {
            Err(GraphMapError::InvalidMap(_msg)) => {
                // Expected error
            }
            Err(GraphMapError::BadEdge(v1, v2, m1, m2)) => {
                panic!(
                    "Expected InvalidMap (wrong length), got BadEdge({v1}, {v2}).\n\
                     Mapping length: {}, Domain vertices: {}",
                    mapping.len(),
                    cube.n()
                );
            }
            Ok(_) => {
                panic!(
                    "Expected InvalidMap error, but GraphMap was created successfully.\n\
                     Mapping: {mapping:?}, Domain vertices: {}",
                    cube.n()
                );
            }
        }
    }

    #[test]
    fn test_graph_map_out_of_range() {
        let cube2 = CubeGraph::<Const<2>>::default();
        let cube3 = CubeGraph::<Const<3>>::default();

        // Try to map 2-cube (4 vertices) to 3-cube (8 vertices) with a vertex out of range
        // Using vertex 10 which doesn't exist in the 3-cube (only has vertices 0-7)
        let mapping = vec![0, 1, 2, 10];
        let invalid_map = VertGraphMap::try_new(&cube2, &cube3, mapping.clone());

        match invalid_map {
            Err(GraphMapError::InvalidMap(_msg)) => {
                // Expected error
            }
            Err(GraphMapError::BadEdge(v1, v2, m1, m2)) => {
                panic!(
                    "Expected InvalidMap (out of range), got {invalid_map:?}.\n\
                     Mapping: {mapping:?}, Codomain vertices: 0-{}",
                    cube3.n() - 1
                );
            }
            Ok(_) => {
                panic!(
                    "Expected InvalidMap error, but GraphMap was created successfully.\n\
                     Mapping: {mapping:?}, Codomain vertices: 0-{}",
                    cube3.n() - 1
                );
            }
        }
    }
}
