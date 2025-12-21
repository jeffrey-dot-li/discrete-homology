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
    fn is_neighbor<V>(&self, a: V, b: V) -> bool
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

        // Fixed
        // Add self-loops and edges from neighbors
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = true; // Self-loop (required but not explicitly stored in CSR)

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
