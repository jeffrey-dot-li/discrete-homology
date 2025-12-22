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
}
