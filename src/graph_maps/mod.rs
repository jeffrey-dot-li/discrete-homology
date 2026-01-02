pub mod cube_isomorphism;
pub mod cube_maps;
pub mod permutation_generator;
pub mod stack_map;
use crate::prelude::*;
use std::borrow::Cow;
use std::fmt::Debug;

pub trait GraphMap<U, V>
where
    U: UGraph,
    V: UGraph,
{
    fn domain(&self) -> &U;
    fn codomain(&self) -> &V;
    fn map(&self, u: u32) -> u32;
    fn mapped_vertices(&self) -> impl Iterator<Item = u32>;

    unsafe fn change_domain(
        &self,
        new_domain: U,
        mapped_vertices: impl IntoIterator<Item = u32>,
    ) -> Self;
}

#[derive(Debug, Clone)]
pub struct VertGraphMap<'u, 'v, U, V>
where
    U: UGraph,
    V: UGraph,
{
    domain: Cow<'u, U>,
    codomain: Cow<'v, V>,
    vert_maps: Vec<u32>,
}

#[derive(Debug)]
pub enum GraphMapError {
    BadEdge(u32, u32, u32, u32),
}

impl<'u, 'v, U, V> VertGraphMap<'u, 'v, U, V>
where
    U: UGraph,
    V: UGraph,
{
    /// # Safety
    ///
    /// Yeah this doesn't actually check if vert_maps are legit
    pub unsafe fn new_unchecked(
        domain: impl Into<Cow<'u, U>>,
        codomain: impl Into<Cow<'v, V>>,
        vert_maps: Cow<'_, Vec<u32>>,
    ) -> Self {
        Self {
            domain: domain.into(),
            codomain: codomain.into(),
            vert_maps: vert_maps.into_owned(),
        }
    }

    pub fn try_from(
        domain: impl Into<Cow<'u, U>>,
        codomain: impl Into<Cow<'v, V>>,
        mapped_verts: impl IntoIterator<Item = u32>,
        // Workspace is a length n mutable array
        workspace: &mut [u32],
    ) -> Result<Self, GraphMapError> {
        let domain = domain.into();
        let codomain = codomain.into();
        let n = domain.n() as usize;

        assert!(workspace.len() >= n);

        use GraphMapError as E;
        let mut max = 0;

        for (i, mapped_i) in mapped_verts.into_iter().enumerate() {
            assert!(
                mapped_i < codomain.n(),
                "invalid map: codomain node {mapped_i} out of range {:?}",
                codomain.n()
            );
            let neighbors = domain.neighbors(i as u32);
            for neighbor in neighbors {
                let neighbor = neighbor as usize;
                if neighbor < i {
                    let mapped_neighbor = workspace[neighbor];
                    if !codomain.is_edge(mapped_neighbor, mapped_i) {
                        return Err(E::BadEdge(
                            i as u32,
                            neighbor as u32,
                            mapped_i,
                            mapped_neighbor,
                        ));
                    }
                }
            }
            // We have confirmed that edges j - i where j <= i are good
            workspace[i] = mapped_i;
            max = i;
        }
        assert!(max == (n - 1), "max should be n - 1, got {max}");
        Ok(Self {
            codomain,
            domain,
            vert_maps: workspace[0..n].to_vec(),
        })
    }
}

pub fn generate_maps_naive<'u, 'v, U: UGraph, V: UGraph>(
    source: &'u U,
    target: &'v V,
) -> (Vec<VertGraphMap<'u, 'v, U, V>>, u64) {
    let n = source.n() as usize;
    let m = target.n() as usize;
    let total_checks = (m as u64).pow(n as u32);
    let mut generator = permutation_generator::PermutationGenerator::new(n as u32, m as u32, 0);

    let mut maps: Vec<VertGraphMap<'u, 'v, U, V>> = Vec::new();
    let mut workspace: Vec<u32> = vec![0; n];

    for _ in 0..total_checks {
        let next = generator.next();

        let next_iter = next.unwrap();
        let valid_map: Result<VertGraphMap<'u, 'v, U, V>, GraphMapError> = VertGraphMap::try_from(
            Cow::Borrowed(source),
            Cow::Borrowed(target),
            next_iter,
            &mut workspace,
        );
        if let Ok(map) = valid_map {
            maps.push(map);
        }
    }
    let last = generator.next();
    debug_assert!(
        last.is_none(),
        "Permutation generator not exhausted {last:?} {generator:?}"
    );
    (maps, total_checks)
}






impl<U, V> GraphMap<U, V> for VertGraphMap<'_, '_, U, V>
where
    U: UGraph,
    V: UGraph,
{
    fn domain(&self) -> &U {
        self.domain.as_ref()
    }
    fn codomain(&self) -> &V {
        self.codomain.as_ref()
    }
    fn map(&self, u: u32) -> u32 {
        self.vert_maps[u as usize]
    }
    fn mapped_vertices(&self) -> impl Iterator<Item = u32> {
        self.vert_maps.iter().copied()
    }

    unsafe fn change_domain(
        &self,
        new_domain: U,
        mapped_vertices: impl IntoIterator<Item = u32>,
    ) -> Self {
        let mapped_verts: Vec<u32> = mapped_vertices.into_iter().collect();
        Self {
            domain: Cow::Owned(new_domain),
            codomain: self.codomain.clone(),
            vert_maps: mapped_verts,
        }
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
        let mut workspace_len2: [u32; 4] = [0; 4];
        let mut gsphere_workspace: Vec<u32> = vec![0; gsphere_graph.n() as usize];

        let id_2cube = VertGraphMap::try_from(
            Cow::Borrowed(&cube),
            Cow::Borrowed(&cube),
            0..cube.n(),
            &mut workspace_len2,
        );
        assert!(id_2cube.is_ok(), "{:?}", id_2cube);

        let id_2gsphere = VertGraphMap::try_from(
            Cow::Borrowed(&gsphere_graph),
            Cow::Borrowed(&gsphere_graph),
            0..gsphere_graph.n(),
            &mut gsphere_workspace,
        );
        assert!(id_2gsphere.is_ok());

        let neg_gsphere = VertGraphMap::try_from(
            Cow::Borrowed(&gsphere_graph),
            Cow::Borrowed(&gsphere_graph),
            0..gsphere_graph.n(),
            &mut gsphere_workspace,
        );
        assert!(neg_gsphere.is_ok());
    }

    #[test]
    fn test_graph_map_bad_edge() {
        let cube = CubeGraph::<Const<2>>::default();
        let mut workspace = [0u32; 4];

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
        let mapping = [0, 3, 1, 2];
        let invalid_map = VertGraphMap::try_from(
            Cow::Borrowed(&cube),
            Cow::Borrowed(&cube),
            mapping,
            &mut workspace,
        );

        match invalid_map {
            Err(GraphMapError::BadEdge(_v1, _v2, _m1, _m2)) => {
                // Expected error
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
    #[should_panic(expected = "max should be n - 1")]
    fn test_graph_map_invalid_length() {
        let cube = CubeGraph::<Const<2>>::default();
        let mut workspace = [0u32; 4];

        // 2-cube has 4 vertices, but we provide only 3 mappings
        // This should panic because the iterator doesn't have enough elements
        let mapping = [0, 1, 2];
        let _invalid_map = VertGraphMap::try_from(
            Cow::Borrowed(&cube),
            Cow::Borrowed(&cube),
            mapping,
            &mut workspace,
        );
    }

    #[test]
    #[should_panic(expected = "invalid map: codomain node 10 out of range")]
    fn test_graph_map_out_of_range() {
        let cube2 = CubeGraph::<Const<2>>::default();
        let cube3 = CubeGraph::<Const<3>>::default();
        let mut workspace = [0u32; 4];

        // Try to map 2-cube (4 vertices) to 3-cube (8 vertices) with a vertex out of range
        // Using vertex 10 which doesn't exist in the 3-cube (only has vertices 0-7)
        // This should panic because vertex 10 is out of range
        let mapping = [0, 1, 2, 10];
        let _invalid_map = VertGraphMap::try_from(
            Cow::Borrowed(&cube2),
            Cow::Borrowed(&cube3),
            mapping,
            &mut workspace,
        );
    }

    #[test]
    fn test_cube2_graphs() {
        let cube2 = CubeGraph::<Const<2>>::default();
        let (maps, total_checks) = generate_maps_naive(&cube2, &cube2);
        assert_eq!(total_checks, 256); // 2^2 = 4 vertices, so 4^4 = 256 total maps checked
        assert_eq!(maps.len(), 84); // There are 24 valid graph maps from 2-cube to 2
    }
}
