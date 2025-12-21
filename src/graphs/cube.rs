use crate::graphs::{AdjMatrix, CSRGraph};

pub fn differ_by_one_bit(a: usize, b: usize) -> bool {
    let x: usize = a ^ b;
    x.is_power_of_two()
}

pub fn n_cube(n: usize) -> CSRGraph {
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
        assert_eq!(cube.neighbors(0), &[0, 1, 2]);
        assert_eq!(cube.neighbors(1), &[0, 1, 3]);
        assert_eq!(cube.neighbors(2), &[0, 2, 3]);
        assert_eq!(cube.neighbors(3), &[1, 2, 3]);
    }
}
