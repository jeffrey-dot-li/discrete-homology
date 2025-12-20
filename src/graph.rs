use core::panic;
use std::{convert::TryFrom, usize};
pub type AdjMatrix = Vec<Vec<bool>>;

pub trait UGraph: TryFrom<AdjMatrix> + Into<AdjMatrix> {
    fn neighbors(&self, v: u32) -> &[u32];
    fn n(&self) -> u32;
}

pub struct CSRGraph {
    offsets: Vec<u32>,
    neighbors: Vec<u32>,
}

impl TryFrom<AdjMatrix> for CSRGraph {
    type Error = String;
    fn try_from(mat: AdjMatrix) -> Result<Self, Self::Error> {
        let n: usize = mat.len();
        if n == 0 {
            return Ok(Self {
                offsets: vec![0],
                neighbors: vec![],
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
        let mut deg: Vec<u32> = vec![0; n];
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
        Ok(Self { offsets, neighbors })
    }
}
impl From<CSRGraph> for AdjMatrix {
    fn from(value: CSRGraph) -> Self {
        panic!()
    }
}

impl UGraph for CSRGraph {
    fn neighbors(&self, v: u32) -> &[u32] {
        let start = self.offsets[v as usize] as usize;
        let end = self.offsets[(v + 1) as usize] as usize;
        &self.neighbors[start..end]
    }

    fn n(&self) -> u32 {
        (self.offsets.len() - 1) as u32
    }
}

impl CSRGraph {
    pub fn new(offsets: Vec<u32>, neighbors: Vec<u32>) -> Self {
        Self { offsets, neighbors }
    }
}
