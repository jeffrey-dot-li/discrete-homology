//! Discrete homology groups from the sparse boundary operators.
//!
//! For a graph `G`, let `C_n` be the free module on the non-degenerate
//! `n`-cubes (`Q_n -> G`) and `d_n : C_n -> C_{n-1}` the alternating face sum.
//! Writing `rank` for the rank over a field `F_p`,
//!
//! ```text
//! beta_n = |C_n| - rank(d_n) - rank(d_{n+1})
//! ```
//!
//! because `ker d_n` has dimension `|C_n| - rank(d_n)` and `im d_{n+1}` has
//! dimension `rank(d_{n+1})`.
//!
//! This is the chain complex of
//! Ebel, Kapulkin, Kershaw, *Discrete homology computations by reduction to
//! zero differentials*, arXiv:2608.04262, Section 1.
//!
//! NOTE ON CONVENTIONS: [`crate::computations::boundary::partial`] realises
//! `sum_i (-1)^i (A o d_i^+ - A o d_i^-)` with `i` zero-indexed, whereas the
//! paper writes `sum_{i=1}^{n} (-1)^i (A o delta_i^0 - A o delta_i^1)`. The two
//! differ by an overall sign (`d_repo = -d_paper`), which does not affect
//! kernels, images, or ranks, hence not the Betti numbers below.

use std::collections::HashMap;

use super::fp::{mod_inv, mod_norm, normalize, sub_scaled};
use crate::computations::boundary::{build_hashmap, partial};
use crate::graph_maps::cube_maps::{combined_cube_maps, CubeMap};
use crate::graph_maps::hashset::UIntOpenHashSet;
use crate::graph_maps::polynomial::Poly;
use crate::graph_maps::stack_map::{generate_maps_naive_stack, StackGraphMap};
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;
use crate::prelude::*;

/// A cube map of any degree into `V`, with the vertex tuple packed as one
/// base-`|V|` digit per domain vertex (see [`StackGraphMap`]).
pub type CubeMapN<'u, 'v, V> = CubeMap<u32, V, StackGraphMap<'u, 'v, CubeGraph<u32>, V, u64>>;

/// One chain level `C_n`: the basis of non-degenerate `n`-cubes.
pub struct Level {
    /// Packed values of the non-degenerate `n`-cubes, sorted ascending. The
    /// position in this vector is the row index of that cube.
    pub basis: Vec<u64>,
    /// Packed value -> row index.
    pub index: HashMap<u64, u32>,
    /// Membership set consumed by [`partial`].
    pub set: UIntOpenHashSet<u64>,
}

impl Level {
    /// Build the level from the (possibly degenerate) `n`-cube maps.
    pub fn from_maps<'u, 'v, V: UGraph>(maps: &[CubeMapN<'u, 'v, V>]) -> Self {
        let mut basis: Vec<u64> = maps
            .iter()
            .filter(|m| !m.is_degenerate())
            .map(|m| m.inner().value())
            .collect();
        basis.sort_unstable();
        basis.dedup();
        let index: HashMap<u64, u32> = basis
            .iter()
            .enumerate()
            .map(|(i, &v)| (v, i as u32))
            .collect();
        let set = build_hashmap(maps.iter().filter(|m| !m.is_degenerate()).map(|m| m.inner()));
        Self { basis, index, set }
    }

    /// `|C_n|`, the number of non-degenerate `n`-cubes.
    pub fn dim(&self) -> usize {
        self.basis.len()
    }
}

/// Generate every `n`-cube map naively (enumerating all `|V|^|Q_n|` tuples).
/// Only usable for small `n`; higher degrees should be built with
/// [`combined_cube_maps`] instead.
pub fn gen_level<'u, 'v, V: UGraph>(
    source: &'u CubeGraph<u32>,
    target: &'v V,
) -> Vec<CubeMapN<'u, 'v, V>> {
    generate_maps_naive_stack(source, target)
        .0
        .into_iter()
        .map(CubeMap::from)
        .collect()
}

/// Rank of a sparse matrix over `F_p`, given as a list of columns. Each column
/// is a list of `(row, coefficient)`; coefficients need not be reduced, ordered,
/// or free of duplicate rows on input (see [`normalize`]). Rows are arbitrary
/// `u32` labels, so the number of rows is not needed here.
///
/// Sparse Gaussian elimination with a pivot per leading row; pivot columns are
/// normalised to leading coefficient 1.
pub fn rank_mod_p(columns: &[Vec<(u32, i64)>], p: i64) -> usize {
    let mut pivots: HashMap<u32, Vec<(u32, i64)>> = HashMap::new();
    for col in columns {
        let mut work = normalize(col, p);

        // Reduce `work` against existing pivots until its leading entry hits an
        // unused row (then it becomes a new pivot) or it vanishes.
        while let Some(&(lead, c)) = work.first() {
            match pivots.get(&lead) {
                Some(pivot) => work = sub_scaled(&work, pivot, c, p),
                None => {
                    let inv = mod_inv(c, p);
                    let normalised: Vec<(u32, i64)> = work
                        .iter()
                        .map(|&(r, cc)| (r, mod_norm(cc * inv, p)))
                        .collect();
                    pivots.insert(lead, normalised);
                    break;
                }
            }
        }
    }
    pivots.len()
}

/// `rank(d_k)`, where `d_k : C_k -> C_{k-1}` is represented by the boundary
/// polynomials of the non-degenerate `k`-cubes, with rows indexed by `lower`.
fn rank_dk<'u, 'v, V: UGraph, const CAP: usize>(
    degree_maps: &[CubeMapN<'u, 'v, V>],
    lower: &Level,
    p: i64,
) -> usize {
    let columns: Vec<Vec<(u32, i64)>> = degree_maps
        .iter()
        .filter(|m| !m.is_degenerate())
        .map(|m| {
            let poly: Poly<u64, CAP> = partial(m.inner(), &lower.set);
            poly.terms()
                .map(|(idx, c)| (lower.index[&idx], c as i64))
                .collect()
        })
        .collect();
    rank_mod_p(&columns, p)
}

/// The chain-level data behind the Betti numbers: level sizes `|C_n|` and the
/// ranks of `d_1, d_2, d_3` over `F_p`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Report {
    /// `|C_0|, |C_1|, |C_2|` -- the number of non-degenerate `n`-cubes.
    pub c: [usize; 3],
    /// `rank(d_1), rank(d_2), rank(d_3)`.
    pub rank_d: [usize; 3],
    /// `beta_0, beta_1, beta_2`.
    pub betti: [i64; 3],
}

/// Compute [`Report`] for `target` over `F_p`.
///
/// Cubes of degree 0..2 are enumerated naively; degree 3 is built by pairing
/// degree-2 cubes ([`combined_cube_maps`]), which is what makes `beta_2`
/// reachable at all (`beta_2` needs all of `d_3`).
///
/// `p` must be prime. A large prime gives the rational Betti numbers whenever
/// the groups are torsion-free in the degrees computed; small primes can only
/// ever shrink a rank if `p`-torsion is present.
pub fn report<V: UGraph>(target: &V, p: i64) -> Report {
    let s0 = CubeGraph::new(0u32);
    let s1 = CubeGraph::new(1u32);
    let s2 = CubeGraph::new(2u32);

    let deg0: Vec<CubeMapN<'_, '_, V>> = gen_level(&s0, target);
    let deg1: Vec<CubeMapN<'_, '_, V>> = gen_level(&s1, target);
    let deg2: Vec<CubeMapN<'_, '_, V>> = gen_level(&s2, target);
    let deg3: Vec<CubeMapN<'_, '_, V>> = combined_cube_maps(&deg2).collect();

    let l0 = Level::from_maps(&deg0);
    let l1 = Level::from_maps(&deg1);
    let l2 = Level::from_maps(&deg2);

    let r1 = rank_dk::<V, 2>(&deg1, &l0, p);
    let r2 = rank_dk::<V, 4>(&deg2, &l1, p);
    let r3 = rank_dk::<V, 6>(&deg3, &l2, p);

    let c = [l0.dim(), l1.dim(), l2.dim()];
    let rank_d = [r1, r2, r3];
    let betti = [
        c[0] as i64 - rank_d[0] as i64,
        c[1] as i64 - rank_d[0] as i64 - rank_d[1] as i64,
        c[2] as i64 - rank_d[1] as i64 - rank_d[2] as i64,
    ];
    Report { c, rank_d, betti }
}

/// Betti numbers `[beta_0, beta_1, beta_2]` of `target` over `F_p`.
pub fn betti_numbers<V: UGraph>(target: &V, p: i64) -> Vec<i64> {
    report(target, p).betti.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphs::extras::{c_n_graph, greene_sphere};

    /// Large prime: avoids `p`-torsion masking ranks for these small graphs.
    const P: i64 = 1_000_003;

    #[test]
    fn test_rank_mod_p_basic() {
        // [[1, 1], [1, 1]] over F_p has rank 1.
        let cols = vec![vec![(0u32, 1i64), (1, 1)], vec![(0, 1), (1, 1)]];
        assert_eq!(rank_mod_p(&cols, P), 1);

        // Identity 2x2.
        let cols = vec![vec![(0u32, 1i64)], vec![(1, 1)]];
        assert_eq!(rank_mod_p(&cols, P), 2);

        // Zero column contributes nothing.
        let cols = vec![vec![], vec![(3u32, 5)]];
        assert_eq!(rank_mod_p(&cols, P), 1);

        // Coefficients reduce mod p: 2 is invertible, and 2*c cancels.
        let cols = vec![vec![(0u32, 2i64), (1, 3)], vec![(0, 1), (1, 3)]];
        assert_eq!(rank_mod_p(&cols, P), 2);

        // Dependent after scaling: col2 = 2 * col1.
        let cols = vec![vec![(0u32, 3i64), (2, 5)], vec![(0, 6), (2, 10)]];
        assert_eq!(rank_mod_p(&cols, P), 1);

        // Duplicate rows within a column must be merged, not treated as two
        // separate entries: [1,0] + [1,0] == 2*[1,0], still rank 1, not 2.
        let cols = vec![vec![(0u32, 1i64), (0, 1)], vec![(1u32, 1i64)]];
        assert_eq!(rank_mod_p(&cols, P), 2);

        // And a pair of duplicate rows that cancels leaves a zero column.
        let cols = vec![vec![(0u32, 1i64), (0, -1)], vec![(1u32, 1i64)]];
        assert_eq!(rank_mod_p(&cols, P), 1);
    }

    /// `d^2 = 0` with degenerate faces dropped: for every non-degenerate
    /// 3-cube, `d_2(d_3(A)) == 0`. This is the theorem that makes the
    /// non-degenerate chain complex a chain complex at all.
    #[test]
    fn test_d_squared_zero_greene_sphere() {
        let target = greene_sphere();
        let s1 = CubeGraph::new(1u32);
        let s2 = CubeGraph::new(2u32);

        let deg1: Vec<CubeMapN<'_, '_, _>> = gen_level(&s1, &target);
        let deg2: Vec<CubeMapN<'_, '_, _>> = gen_level(&s2, &target);
        let deg3: Vec<CubeMapN<'_, '_, _>> = combined_cube_maps(&deg2).collect();

        let l1 = Level::from_maps(&deg1);
        let l2 = Level::from_maps(&deg2);

        let by_val2: HashMap<u64, &CubeMapN<'_, '_, _>> = deg2
            .iter()
            .filter(|m| !m.is_degenerate())
            .map(|m| (m.inner().value(), m))
            .collect();

        let mut checked = 0usize;
        for m3 in deg3.iter().filter(|m| !m.is_degenerate()) {
            let b3: Poly<u64, 6> = partial(m3.inner(), &l2.set);
            let mut acc: HashMap<u64, i64> = HashMap::new();
            for (v2, c3) in b3.terms() {
                let m2 = by_val2[&v2];
                let b2: Poly<u64, 4> = partial(m2.inner(), &l1.set);
                for (v1, c2) in b2.terms() {
                    *acc.entry(v1).or_insert(0) += (c3 as i64) * (c2 as i64);
                }
            }
            assert!(
                acc.values().all(|&c| c == 0),
                "d^2 != 0 on 3-cube {:?}: {:?}",
                m3.inner().value(),
                acc
            );
            checked += 1;
        }
        assert!(checked > 0, "no non-degenerate 3-cubes were checked");
    }

    /// The paper's headline low-degree fact: `dim H_2(Gsph) = 1`, Gsph being the
    /// smallest graph with non-trivial second homology.
    #[test]
    fn test_betti_greene_sphere() {
        let r = report(&greene_sphere(), P);
        // 10 vertices, all on one component.
        assert_eq!(r.c[0], 10, "|C_0|, got {r:?}");
        assert_eq!(r.rank_d[0], 9, "rank(d_1) = |V| - beta_0, got {r:?}");
        // beta_0 = 1 (connected), beta_1 = 0 (Gsph is a discrete 2-sphere).
        assert_eq!(r.betti[0], 1, "beta_0(Gsph), got {r:?}");
        assert_eq!(r.betti[1], 0, "beta_1(Gsph), got {r:?}");
        assert_eq!(
            r.betti[2], 1,
            "dim H_2(Gsph) should be 1 (Ebel-Kapulkin-Kershaw, sec. 1/4), got {r:?}"
        );
    }

    /// Prints chain level sizes, boundary ranks, and Betti numbers for the
    /// graphs shipped with the crate:
    /// `cargo test --lib dump_report -- --nocapture`.
    #[test]
    fn dump_report() {
        let cases = [
            ("C3", c_n_graph(3)),
            ("C4", c_n_graph(4)),
            ("C5", c_n_graph(5)),
            ("C6", c_n_graph(6)),
            ("C7", c_n_graph(7)),
            ("Gsph", greene_sphere()),
        ];
        eprintln!(
            "{:<6} {:>4} {:>7} {:>7} | {:>6} {:>6} {:>6} | {:>4} {:>4} {:>4}",
            "graph", "|C0|", "|C1|", "|C2|", "r(d1)", "r(d2)", "r(d3)", "b0", "b1", "b2"
        );
        for (name, g) in cases {
            let r = report(&g, P);
            eprintln!(
                "{:<6} {:>4} {:>7} {:>7} | {:>6} {:>6} {:>6} | {:>4} {:>4} {:>4}",
                name,
                r.c[0],
                r.c[1],
                r.c[2],
                r.rank_d[0],
                r.rank_d[1],
                r.rank_d[2],
                r.betti[0],
                r.betti[1],
                r.betti[2]
            );
        }
    }

    /// `H_1(C_n) = Z` for `n >= 4`: discrete homology fills triangles and
    /// squares, so the first cycle it leaves open has length 5.
    #[test]
    fn test_betti_cycles() {
        let c5 = betti_numbers(&c_n_graph(5), P);
        assert_eq!(c5[0], 1, "beta_0(C5), got {c5:?}");
        assert_eq!(c5[1], 1, "beta_1(C5), got {c5:?}");

        let c4 = betti_numbers(&c_n_graph(4), P);
        assert_eq!(c4[0], 1, "beta_0(C4), got {c4:?}");
        assert_eq!(c4[1], 0, "beta_1(C4) should vanish: 4-cycles bound, got {c4:?}");

        let c3 = betti_numbers(&c_n_graph(3), P);
        assert_eq!(c3[0], 1, "beta_0(C3), got {c3:?}");
        assert_eq!(c3[1], 0, "beta_1(C3) should vanish: 3-cycles bound, got {c3:?}");

        let c6 = betti_numbers(&c_n_graph(6), P);
        assert_eq!(c6[1], 1, "beta_1(C6), got {c6:?}");
    }
}
