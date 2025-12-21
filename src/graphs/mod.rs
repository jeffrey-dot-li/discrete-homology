pub mod cube;
pub mod extras;

use std::convert::TryFrom;
pub type AdjMatrix = Vec<Vec<bool>>;

// Simple Undirected Graph
pub trait UGraph: TryFrom<AdjMatrix> + Into<AdjMatrix> {
    fn neighbors<V>(&self, v: V) -> &[u32]
    where
        V: Into<u32>;
    fn n(&self) -> u32;

    fn degree(&self, v: u32) -> u32 {
        self.neighbors(v).len() as u32
    }
    fn is_edge<V>(&self, a: V, b: V) -> bool
    where
        V: Into<u32>,
    {
        // This is actually suboptimal because we know that neighbors is sorted
        self.neighbors(a).contains(&b.into())
    }
}
#[derive(Debug, Clone)]
pub struct CSRGraph {
    offsets: Vec<u32>,
    neighbor_list: Vec<u32>,
}

impl TryFrom<AdjMatrix> for CSRGraph {
    type Error = String;
    fn try_from(mat: AdjMatrix) -> Result<Self, Self::Error> {
        let n: usize = mat.len();
        if n == 0 {
            return Ok(Self {
                offsets: vec![0],
                neighbor_list: vec![],
            });
        }
        for (i, row) in mat.iter().enumerate() {
            if row.len() != n {
                return Err(format!(
                    "adjacency matrix must be square: row {i} has len {}, expected {n}",
                    row.len()
                ));
            }
        }

        // Symmetry + diagonal check
        for i in 0..n {
            if !mat[i][i] {
                return Err(format!("Must have identity relation on node mat[{i}][{i}]"));
            }
            for j in (i + 1)..n {
                if mat[i][j] != mat[j][i] {
                    return Err(format!(
                        "matrix must be symmetric mat[{i}][{j}] = {} but mat[{j}][{i}] = {}",
                        mat[i][j], mat[j][i]
                    ));
                }
            }
        }
        // Start with 1 degree for self-loop
        let mut deg: Vec<u32> = vec![1; n];
        for i in 0..n {
            for j in (i + 1)..n {
                if mat[i][j] {
                    deg[i] = deg[i]
                        .checked_add(1)
                        .ok_or_else(|| format!("degree overflow at node {i}"))?;
                    deg[j] = deg[j]
                        .checked_add(1)
                        .ok_or_else(|| format!("degree overflow at node {j}"))?;
                }
            }
        }

        // Prefix sum -> offsets

        let mut offsets: Vec<u32> = vec![0; n + 1];
        for i in 0..n {
            offsets[i + 1] = offsets[i]
                .checked_add(deg[i])
                .ok_or_else(|| "offset overflow (too many edges for u32)".to_string())?;
        }
        let total: u32 = offsets[n];
        let mut neighbors: Vec<u32> = vec![0; total as usize];
        let mut cursor: Vec<u32> = offsets[..n].to_vec();

        for i in 0..n {
            // Add self-loop
            neighbors[cursor[i] as usize] = i as u32;
            cursor[i] += 1;

            for j in (i + 1)..n {
                if mat[i][j] {
                    let iu = i as u32;
                    let ju = j as u32;

                    let pi = cursor[i] as usize;
                    neighbors[pi] = ju;
                    cursor[i] += 1;

                    let pj = cursor[j] as usize;
                    neighbors[pj] = iu;
                    cursor[j] += 1;
                }
            }
        }

        // Sort each adjacency list so neighbor iteration is deterministic - least to greatest

        for u in 0..n {
            let start = offsets[u] as usize;
            let end = offsets[u + 1] as usize;
            neighbors[start..end].sort_unstable();
        }
        Ok(Self {
            offsets,
            neighbor_list: neighbors,
        })
    }
}
impl From<CSRGraph> for AdjMatrix {
    fn from(value: CSRGraph) -> Self {
        let n = value.offsets.len() - 1;
        let mut mat = vec![vec![false; n]; n];

        for (i, row) in mat.iter_mut().enumerate() {
            let start = value.offsets[i] as usize;
            let end = value.offsets[i + 1] as usize;
            for &neighbor in &value.neighbor_list[start..end] {
                row[neighbor as usize] = true;
            }
        }

        mat
    }
}

impl UGraph for CSRGraph {
    fn neighbors<V>(&self, v: V) -> &[u32]
    where
        V: Into<u32>,
    {
        let v: u32 = v.into();
        let start = self.offsets[v as usize] as usize;
        let end = self.offsets[(v + 1) as usize] as usize;
        &self.neighbor_list[start..end]
    }

    fn n(&self) -> u32 {
        (self.offsets.len() - 1) as u32
    }
}

impl CSRGraph {
    pub fn new(offsets: Vec<u32>, neighbors: Vec<u32>) -> Self {
        Self {
            offsets,
            neighbor_list: neighbors,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tryfrom_non_square_matrix() {
        // Create a non-square matrix (3 rows but row 1 has only 2 columns)
        let mat: AdjMatrix = vec![
            vec![true, false, true],
            vec![false, true], // Wrong length!
            vec![true, false, true],
        ];

        let result = CSRGraph::try_from(mat.clone());

        match result {
            Err(err) => {
                assert!(
                    err.contains("square"),
                    "Error should mention 'square', got: {err}"
                );
            }
            Ok(_) => {
                panic!(
                    "Expected error for non-square matrix, but CSRGraph was created successfully.\n\
                     Matrix dimensions: {} rows\n\
                     Row 0 length: {}\n\
                     Row 1 length: {} (should be {})",
                    mat.len(),
                    mat[0].len(),
                    mat[1].len(),
                    mat.len()
                );
            }
        }
    }

    #[test]
    fn test_tryfrom_missing_self_loop() {
        // Create a matrix missing a self-loop at position [1][1]
        let mat: AdjMatrix = vec![
            vec![true, false, false],
            vec![false, false, false], // Missing self-loop at [1][1]!
            vec![false, false, true],
        ];

        let result = CSRGraph::try_from(mat.clone());

        match result {
            Err(err) => {
                assert!(
                    err.contains("identity") || err.contains("mat[1][1]"),
                    "Error should mention identity relation or mat[1][1], got: {err}"
                );
            }
            Ok(_) => {
                let diagonal_info: Vec<String> = mat
                    .iter()
                    .enumerate()
                    .map(|(i, row)| format!("mat[{i}][{i}] = {}", row[i]))
                    .collect();
                panic!(
                    "Expected error for missing self-loop, but CSRGraph was created successfully.\n\
                     Diagonal values:\n  {}",
                    diagonal_info.join("\n  ")
                );
            }
        }
    }

    #[test]
    fn test_tryfrom_non_symmetric_matrix() {
        // Create a non-symmetric matrix where mat[0][1] != mat[1][0]
        let mat: AdjMatrix = vec![
            vec![true, true, false],  // mat[0][1] = true
            vec![false, true, false], // mat[1][0] = false (asymmetric!)
            vec![false, false, true],
        ];

        let result = CSRGraph::try_from(mat.clone());

        match result {
            Err(err) => {
                assert!(
                    err.contains("symmetric"),
                    "Error should mention 'symmetric', got: {err}"
                );
            }
            Ok(_) => {
                panic!(
                    "Expected error for non-symmetric matrix, but CSRGraph was created successfully.\n\
                     Asymmetry at: mat[0][1] = {}, mat[1][0] = {}",
                    mat[0][1], mat[1][0]
                );
            }
        }
    }

    #[test]
    fn test_tryfrom_valid_matrix() {
        // Test that a valid matrix successfully converts
        let mat: AdjMatrix = vec![
            vec![true, true, false],
            vec![true, true, true],
            vec![false, true, true],
        ];

        let result = CSRGraph::try_from(mat.clone());

        match result {
            Ok(graph) => {
                assert_eq!(graph.n(), 3, "Graph should have 3 vertices");

                // Verify the graph structure matches the matrix
                for i in 0..3u32 {
                    let neighbors = graph.neighbors(i);

                    for j in 0..3u32 {
                        let should_be_neighbor = mat[i as usize][j as usize];
                        let is_neighbor = neighbors.contains(&j);
                        assert_eq!(
                            is_neighbor, should_be_neighbor,
                            "Vertex {i} neighbor {j}: expected {should_be_neighbor}, got {is_neighbor}. Neighbors: {neighbors:?}"
                        );
                    }
                }
            }
            Err(err) => {
                panic!("Expected successful conversion for valid matrix, got error: {err}");
            }
        }
    }

    #[test]
    fn test_tryfrom_empty_matrix() {
        // Test edge case: empty matrix
        let mat: AdjMatrix = vec![];

        let result = CSRGraph::try_from(mat);

        match result {
            Ok(graph) => {
                assert_eq!(graph.n(), 0, "Empty graph should have 0 vertices");
            }
            Err(err) => {
                panic!("Expected successful conversion for empty matrix, got error: {err}");
            }
        }
    }
}
