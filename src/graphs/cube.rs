use crate::graphs::{AdjMatrix, UGraph};
use crate::prelude::*;
use crate::shape::{Const, Dim};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CubeGraph<D: Dim> {
    n: D,
}

impl<const N: u32> Default for CubeGraph<Const<N>> {
    fn default() -> Self {
        Self { n: Const::<N> }
    }
}

impl CubeGraph<u32> {
    pub fn new(n: u32) -> Self {
        Self { n }
    }
}

impl<D: Dim> From<CubeGraph<D>> for AdjMatrix {
    fn from(value: CubeGraph<D>) -> Self {
        let n = value.n.size();
        let verts = 2_usize.pow(n);
        let mut adj: AdjMatrix = vec![vec![false; verts]; verts];
        for i in 0..verts {
            adj[i][i] = true;
            for j in (i + 1)..verts {
                if cube_share_edge(i, j) {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }
        adj
    }
}

impl<D: Dim> GraphNeighbors for CubeGraph<D> {
    fn neighbors(&self, v: u32) -> impl Iterator<Item = u32> {
        let mut items = ((-1i32)..(self.n.size() as i32))
            .map(move |i| {
                if i < 0 {
                    return v;
                }

                // Change one bit at the ith position
                2u32.pow(i as u32) ^ v
            })
            .collect::<Vec<_>>();
        items.sort();
        items.into_iter()
    }
}
impl<D: Dim> UGraph for CubeGraph<D> {
    fn n(&self) -> u32 {
        2_u32.pow(self.n.size())
    }

    fn degree(&self, _v: u32) -> u32 {
        self.n.size()
    }

    fn is_edge<V: Into<u32>>(&self, a: V, b: V) -> bool {
        let a: u32 = a.into();
        let b: u32 = b.into();
        a == b || cube_share_edge(a as usize, b as usize)
    }
}

pub fn cube_share_edge(a: usize, b: usize) -> bool {
    if a == b {
        return true;
    }
    // let (a, b) = ordered(a, b);

    let x: usize = a ^ b;
    x.is_power_of_two()
}

#[cfg(test)]
mod tests {
    use crate::graphs::UGraph;

    use super::*;

    #[test]
    fn test_differ_by_one_bit() {
        // assert!(differ_by_one_bit(0, 0));
        assert!(cube_share_edge(0, 1));
        assert!(cube_share_edge(1, 0));
        assert!(cube_share_edge(3, 1));
        assert!(!cube_share_edge(1, 2));
    }

    #[test]
    fn test_n_cube() {
        let cube = CubeGraph::<Const<2>>::default();

        assert_eq!(cube.n(), 4);

        // Define expected neighbors for each vertex
        let expected = [
            (0u32, vec![0, 1, 2]),
            (1u32, vec![0, 1, 3]),
            (2u32, vec![0, 2, 3]),
            (3u32, vec![1, 2, 3]),
        ];
        let expected_neighbors = expected
            .iter()
            .map(|(_v, neighbors)| neighbors.clone())
            .collect::<Vec<_>>();

        let cube_neighbors = (0..4)
            .map(|v| cube.neighbors(v).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        assert_eq!(cube_neighbors, expected_neighbors);
        // Check that all edges match expectations using is_edge
        assert!(
            expected.iter().all(|(v, neighbors)| {
                (0..4).all(|n| cube.is_edge(*v, n) == neighbors.contains(&n))
            }),
            "Edge structure doesn't match expected 2-cube connectivity"
        );
    }

    #[test]
    fn test_into_adj_matrix() {
        let cube = CubeGraph::<Const<3>>::default();
        let n = cube.n() as usize;

        // Convert CSRGraph to AdjMatrix
        let adj_matrix: AdjMatrix = cube.into();

        // Check matrix dimensions
        assert_eq!(adj_matrix.len(), n, "Matrix should have {n} rows");
        for (i, row) in adj_matrix.iter().enumerate() {
            assert_eq!(row.len(), n, "Row {i} should have {n} columns");
        }

        // Check all self-loops are present (diagonal is true)
        for (i, row) in adj_matrix.iter().enumerate() {
            assert!(row[i], "Self-loop must exist at ({i}, {i})");
        }

        // Check symmetry
        for i in 0..n {
            for j in 0..n {
                assert_eq!(
                    adj_matrix[i][j], adj_matrix[j][i],
                    "Matrix must be symmetric: ({i}, {j}) != ({j}, {i})"
                );
            }
        }

        // Check edges are correct (vertices differ by one bit should be connected)
        for (i, row) in adj_matrix.iter().enumerate() {
            for (j, is_connected) in row.iter().enumerate() {
                let should_be_connected = i == j || cube_share_edge(i, j);
                assert_eq!(
                    *is_connected, should_be_connected,
                    "Edge ({i}, {j}): expected {should_be_connected}, got {}",
                    *is_connected
                );
            }
        }
    }

    #[test]
    fn test_roundtrip_conversion() {
        // Test that AdjMatrix -> CSRGraph -> AdjMatrix preserves the structure
        let original_cube = CubeGraph::<Const<2>>::default();
        let adj1: AdjMatrix = original_cube.into();

        let csr = CSRGraph::try_from(adj1.clone()).expect("Should convert to CSRGraph");
        let adj2: AdjMatrix = csr.into();

        // Check that both matrices are identical
        assert_eq!(
            adj1, adj2,
            "Roundtrip conversion should preserve matrix structure"
        );
    }
}
