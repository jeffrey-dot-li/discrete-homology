use crate::graphs::{AdjMatrix, CSRGraph, UGraph};

pub struct CubeGraph {
    g: CSRGraph,
}
pub fn differ_by_one_bit(a: usize, b: usize) -> bool {
    let x: usize = a ^ b;
    x.is_power_of_two()
}

pub fn n_cube(n: usize) -> impl UGraph {
    let verts = 2_usize.pow(n as u32);

    let mut adj: AdjMatrix = vec![vec![false; verts]; verts];

    for i in 0..verts {
        adj[i][i] = true;
        for j in (i + 1)..verts {
            if differ_by_one_bit(i, j) {
                adj[i][j] = true;
                adj[j][i] = true;
            }
        }
    }
    CSRGraph::try_from(adj).unwrap()
}

#[cfg(test)]
mod tests {
    use crate::graphs::UGraph;

    use super::*;

    #[test]
    fn test_n_cube() {
        let cube = n_cube(2);
        assert_eq!(cube.neighbors(0u32), &[0, 1, 2]);
        assert_eq!(cube.neighbors(1u32), &[0, 1, 3]);
        assert_eq!(cube.neighbors(2u32), &[0, 2, 3]);
        assert_eq!(cube.neighbors(3u32), &[1, 2, 3]);
    }

    #[test]
    fn test_into_adj_matrix() {
        let cube = n_cube(3);
        let n = cube.n() as usize;

        // Convert CSRGraph to AdjMatrix
        let adj_matrix: AdjMatrix = cube.into();

        // Check matrix dimensions
        assert_eq!(adj_matrix.len(), n, "Matrix should have {n} rows");
        for (i, row) in adj_matrix.iter().enumerate() {
            assert_eq!(row.len(), n, "Row {i} should have {n} columns");
        }

        // Check all self-loops are present (diagonal is true)
        for i in 0..n {
            assert!(adj_matrix[i][i], "Self-loop must exist at ({i}, {i})");
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
        for i in 0..n {
            for j in 0..n {
                let should_be_connected = i == j || differ_by_one_bit(i, j);
                assert_eq!(
                    adj_matrix[i][j], should_be_connected,
                    "Edge ({i}, {j}): expected {should_be_connected}, got {}",
                    adj_matrix[i][j]
                );
            }
        }
    }

    #[test]
    fn test_roundtrip_conversion() {
        // Test that AdjMatrix -> CSRGraph -> AdjMatrix preserves the structure
        let original_cube = n_cube(2);
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
