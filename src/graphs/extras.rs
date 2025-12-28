use crate::prelude::*;
pub fn greene_sphere() -> CSRGraph {
    let t = true;
    let f = false;
    let adj: AdjMatrix = vec![
        vec![t, f, t, f, t, f, t, f, t, f],
        vec![f, t, t, f, f, f, f, f, t, t],
        vec![t, t, t, t, f, f, f, f, f, f],
        vec![f, f, t, t, t, f, f, f, f, t],
        vec![t, f, f, t, t, t, f, f, f, f],
        vec![f, f, f, f, t, t, t, f, f, t],
        vec![t, f, f, f, f, t, t, t, f, f],
        vec![f, f, f, f, f, f, t, t, t, t],
        vec![t, t, f, f, f, f, f, t, t, f],
        vec![f, t, f, t, f, t, f, t, f, t],
    ];
    CSRGraph::try_from(adj).unwrap()
}

pub fn C_N_graph(n: u32) -> CSRGraph {
    let mut adj: AdjMatrix = vec![vec![false; n as usize]; n as usize];
    for i in 0..n {
        adj[i as usize][((i + 1) % n) as usize] = true;
        adj[i as usize][((i + n - 1) % n) as usize] = true;
        adj[i as usize][i as usize] = true;
    }
    CSRGraph::try_from(adj).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_greene_sphere() {
        let greene_sphere = greene_sphere();
        let n = greene_sphere.n();
        assert_eq!(
            (greene_sphere.neighbors(0u32)).collect::<Vec<_>>(),
            vec![0, 2, 4, 6, 8]
        );
        assert_eq!(
            greene_sphere.neighbors(n - 1).collect::<Vec<_>>(),
            vec![1, 3, 5, 7, 9]
        );

        for i in 1..n - 1 {
            let prev = (i + n - 3) % ((n - 2) as i32) as u32 + 1;
            let next = ((i) % (n - 2)) as u32 + 1;
            let s = if i % 2 == 0 { 0 } else { n - 1 };
            assert_eq!(greene_sphere.degree(i), 4);
            assert!(
                greene_sphere.is_edge(i, prev),
                "{i}, {prev} were not neighbors"
            );
            assert!(greene_sphere.is_edge(i, next));
            assert!(greene_sphere.is_edge(i, s));
        }
    }

    #[test]
    fn test_c_n_graph() {
        // Test C_3 (triangle cycle)
        let c3 = C_N_graph(3);
        assert_eq!(c3.n(), 3, "C_3 should have 3 vertices");

        // Each vertex in C_3 should have degree 3 (self-loop + 2 neighbors)
        for i in 0..3 {
            assert_eq!(c3.degree(i), 3, "Vertex {i} in C_3 should have degree 3");

            // Check self-loop
            assert!(c3.is_edge(i, i), "Vertex {i} should have a self-loop");

            // Check neighbors
            let next = (i + 1) % 3;
            let prev = (i + 2) % 3; // (i - 1 + 3) % 3
            assert!(
                c3.is_edge(i, next),
                "Vertex {i} should be connected to {next}"
            );
            assert!(
                c3.is_edge(i, prev),
                "Vertex {i} should be connected to {prev}"
            );
        }

        // Test C_4 (square cycle)
        let c4 = C_N_graph(4);
        assert_eq!(c4.n(), 4, "C_4 should have 4 vertices");

        for i in 0..4 {
            assert_eq!(c4.degree(i), 3, "Vertex {i} in C_4 should have degree 3");

            // Check that vertex is only connected to self, next, and previous
            let next = (i + 1) % 4;
            let prev = (i + 3) % 4; // (i - 1 + 4) % 4

            for j in 0..4 {
                let should_be_connected = j == i || j == next || j == prev;
                assert_eq!(
                    c4.is_edge(i, j),
                    should_be_connected,
                    "C_4: Edge ({i}, {j}) should be {should_be_connected}"
                );
            }
        }

        // Test C_5 (pentagon cycle)
        let c5 = C_N_graph(5);
        assert_eq!(c5.n(), 5, "C_5 should have 5 vertices");

        // Verify cycle structure
        let expected_neighbors: Vec<Vec<u32>> = vec![
            vec![0, 1, 4], // Vertex 0: self, next (1), prev (4)
            vec![0, 1, 2], // Vertex 1: prev (0), self, next (2)
            vec![1, 2, 3], // Vertex 2: prev (1), self, next (3)
            vec![2, 3, 4], // Vertex 3: prev (2), self, next (4)
            vec![0, 3, 4], // Vertex 4: prev (3), self, next (0)
        ];

        for (i, expected) in expected_neighbors.iter().enumerate() {
            let actual: Vec<u32> = c5.neighbors(i as u32).collect();
            assert_eq!(
                actual, *expected,
                "C_5: Vertex {i} neighbors mismatch. Expected: {expected:?}, Got: {actual:?}"
            );
        }

        // Verify symmetry
        for i in 0u32..5 {
            for j in 0u32..5 {
                assert_eq!(
                    c5.is_edge(i, j),
                    c5.is_edge(j, i),
                    "C_5 should be symmetric: edge ({i}, {j}) != edge ({j}, {i})"
                );
            }
        }
    }
}
