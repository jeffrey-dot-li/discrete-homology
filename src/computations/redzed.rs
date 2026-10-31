//! Base RedZeD and active enumeration for non-persistent discrete homology.
//!
//! Reference: Ebel, Kapulkin, Kershaw, *Discrete homology computations by
//! reduction to zero differentials*, arXiv:2608.04262, Section 1 (RedZeD) and
//! Section 2 (RedZeD DH).
//!
//! The chain complex is `... -> C_p --d_p--> C_{p-1} -> ...`, with `C_p` free on
//! the non-degenerate `p`-cubes. [`ColumnReducer`] is the algebraic core: it
//! reduces the boundary matrix `d_p` one column at a time over `F_p`, keeping
//!
//! * `pivots`, the reduced columns indexed by their leading row (this is the
//!   `R_p^{-1}` side of the paper's dictionaries, and `pivots.len()` is
//!   `rank(d_p)`), and
//! * `births`, the generators whose column reduced to zero (generators with
//!   `r_p != 0`, i.e. the *active* ones, which are exactly `keys(R_p)`).
//!
//! Rows are not the degree-`(p-1)` generators themselves but the **births** of
//! degree `p-1`, i.e. the basis of `A_{p-1}`: the column of a degree-`p`
//! generator is `r_{p-1}(d g)`, so a face that is a death in degree `p-1`
//! contributes nothing. This is what makes every leading row a birth, so a
//! generator `g` in degree `p` is a **birth** if its column reduces to zero
//! and a **death** otherwise, in which case it kills the class born at
//! `births_{p-1}[row]`. Then `beta_p = |births_p| - rank(d_{p+1})`.
//!
//! Because `r_{p-1}` is injective on `ker d_{p-1}` and `im d_p` lies inside it,
//! `rank(r_{p-1} o d_p) = rank(d_p)`, so the ranks -- and hence the Betti
//! numbers -- agree with the direct computation in [`super::homology`]. That
//! equality is what the `test_redzed_matches_rank_betti` test pins down.
//!
//! Writing that out: `|births_p| = |C_p| - rank(d_p) = dim ker d_p`, so
//! `beta_p = dim ker d_p - rank(d_{p+1})`, the usual formula -- but computed
//! incrementally, which is what makes the second phase possible.
//!
//! ## Active enumeration
//!
//! A `(p+1)`-cube whose faces are all inactive has `r d = 0`, so it is a birth
//! one degree above the range of interest and cannot affect `H_<=p`. It is
//! therefore never enumerated: `run_phase2` only pairs *active* `p`-cubes with
//! the rest, and abandons a cube as soon as it ceases to be active. The result
//! is unchanged; only the number of top-dimensional cubes constructed drops.
//!
//! ## Caveat
//!
//! This is the unquotiented variant: cubes are not grouped into orbits under
//! `Aut(Q_n) = S_n |x (Z_2)^n` and semi-degenerate cubes are not discarded, so
//! `C_p` here is larger than the complex the paper reduces. It computes the same
//! homology, at proportionally higher cost.

use std::collections::{HashMap, HashSet};

use super::fp::{mod_inv, mod_norm, normalize, sub_scaled};
use super::homology::{gen_level, CubeMapN, Level};
use crate::computations::boundary::partial;
use crate::graph_maps::cube_maps::combined_cube_maps;
use crate::graph_maps::polynomial::Poly;
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;

/// Result of registering one column with a [`ColumnReducer`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// The column reduced to zero: this generator creates a class in its degree.
    Birth,
    /// The column reduced to a nonzero vector with leading row `row`: this
    /// generator kills the class born at the generator at that row.
    Death { row: u32 },
}

/// Incremental sparse reduction of one boundary matrix `d_p` over `F_p`.
///
/// Columns are degree-`p` generators (identified by their packed value) and rows
/// are degree-`(p-1)` generators (identified by their index in the level's
/// basis).
pub struct ColumnReducer {
    p: i64,
    /// Leading row -> (generator, reduced column normalised to lead 1).
    pivots: HashMap<u32, (u64, Vec<(u32, i64)>)>,
    /// Generators whose column reduced to zero, in processing order.
    births: Vec<u64>,
    /// Number of columns registered (births and deaths alike).
    pushed: usize,
}

impl ColumnReducer {
    pub fn new(p: i64) -> Self {
        Self {
            p,
            pivots: HashMap::new(),
            births: Vec::new(),
            pushed: 0,
        }
    }

    /// Reduce and register the column of generator `g`, whose boundary is
    /// `col` (rows as `(row, coefficient)`, arbitrary order).
    pub fn push(&mut self, g: u64, col: &[(u32, i64)]) -> Outcome {
        self.pushed += 1;
        let mut work = normalize(col, self.p);

        while let Some(&(lead, c)) = work.first() {
            match self.pivots.get(&lead) {
                Some((_, pivot)) => work = sub_scaled(&work, pivot, c, self.p),
                None => {
                    let inv = mod_inv(c, self.p);
                    let normalised: Vec<(u32, i64)> = work
                        .iter()
                        .map(|&(r, cc)| (r, mod_norm(cc * inv, self.p)))
                        .collect();
                    self.pivots.insert(lead, (g, normalised));
                    return Outcome::Death { row: lead };
                }
            }
        }

        self.births.push(g);
        Outcome::Birth
    }

    /// `rank(d_p)`, the number of death generators found so far.
    pub fn rank(&self) -> usize {
        self.pivots.len()
    }

    /// The active generators: those with a nonzero image, in processing order.
    pub fn births(&self) -> &[u64] {
        &self.births
    }

    pub fn columns_pushed(&self) -> usize {
        self.pushed
    }
}

/// Per-degree outcome of [`redzed_dh`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RedzedOutcome {
    /// `beta_0 ..= beta_n`.
    pub betti: Vec<i64>,
    /// `|births_p| = dim ker d_p` for `p = 0 ..= n`.
    pub births: Vec<usize>,
    /// `rank(d_p)` for `p = 1 ..= n+1`.
    pub rank_d: Vec<usize>,
    /// Number of `(n+1)`-cubes enumerated in the top phase. With active
    /// enumeration on, this drops by the number of top cubes whose faces are
    /// all inactive.
    pub top_enumerated: usize,
}

/// Row index for a set of generators: packed value -> row.
fn row_index(values: &[u64]) -> HashMap<u64, u32> {
    values
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i as u32))
        .collect()
}

/// Build the column of a degree-`deg` cube map as an element of `A_{deg-1}`.
///
/// `CAP` must equal `2 * deg`, matching the number of faces. `lower` supplies
/// the non-degenerate `(deg-1)`-cube membership set that `partial` needs to
/// decide which faces are generators at all; `rows` then keeps only the faces
/// that are *births* in degree `deg-1`, which is the map `r_{deg-1}`.
fn boundary_column<const CAP: usize, V: UGraph>(
    map: &CubeMapN<'_, '_, V>,
    lower: &Level,
    rows: &HashMap<u64, u32>,
) -> Vec<(u32, i64)> {
    let poly: Poly<u64, CAP> = partial(map.inner(), &lower.set);
    poly.terms()
        .filter_map(|(v, c)| rows.get(&v).map(|&r| (r, c as i64)))
        .collect()
}

/// Phase 1: register every non-degenerate degree-`deg` cube, in enumeration
/// order. `rows` is the row index of the births one degree below.
fn run_phase1<const CAP: usize, V: UGraph>(
    maps: &[CubeMapN<'_, '_, V>],
    lower: &Level,
    rows: &HashMap<u64, u32>,
    p: i64,
) -> ColumnReducer {
    let mut reducer = ColumnReducer::new(p);
    for m in maps.iter().filter(|m| !m.is_degenerate()) {
        let col = boundary_column::<CAP, V>(m, lower, rows);
        reducer.push(m.inner().value(), &col);
    }
    reducer
}

/// Phase 2: enumerate the top-degree cubes by pairing the degree below, and
/// register the resulting columns.
///
/// With `active` set, only pairs where the outer cube is currently active are
/// followed, and the outer loop moves on as soon as that cube dies. Returns the
/// reducer and the number of top cubes actually enumerated.
///
/// Every axis is tried, not just the appended one. The paper's enumeration
/// ranges the positive face over the whole hyperoctahedral orbit of `B`, which
/// lets it place the fixed active cube on *any* face; without the quotient, the
/// equivalent is to build the gluing with the new coordinate inserted at each
/// position in turn. Pairing along the appended axis only finds a top cube via
/// its two last-coordinate faces, so a cube whose only live face lies along an
/// earlier axis is lost -- and if it was a death, so is its pivot, which
/// inflates the Betti number by exactly that count.
fn run_phase2<const CAP: usize, V: UGraph>(
    maps: &[CubeMapN<'_, '_, V>],
    lower: &Level,
    births: &[u64],
    rows: &HashMap<u64, u32>,
    p: i64,
    active: bool,
) -> (ColumnReducer, usize) {
    // The cubes being paired are degree `CAP / 2 - 1`; their gluing adds one
    // coordinate, at any of `dim + 1` positions.
    let dim = (CAP as u32) / 2 - 1;

    let mut reducer = ColumnReducer::new(p);
    let mut processed: HashSet<u64> = HashSet::new();
    let mut enumerated = 0usize;
    let mut still_active: HashSet<u64> = births.iter().copied().collect();

    let index_of: HashMap<u64, usize> = maps
        .iter()
        .enumerate()
        .filter(|(_, m)| !m.is_degenerate())
        .map(|(i, m)| (m.inner().value(), i))
        .collect();

    let outers: Vec<usize> = if active {
        births.iter().map(|v| index_of[v]).collect()
    } else {
        (0..maps.len()).filter(|&i| !maps[i].is_degenerate()).collect()
    };

    'outer: for ai in outers {
        let a_val = maps[ai].inner().value();
        if active && !still_active.contains(&a_val) {
            continue;
        }

        // The partner is deliberately *not* filtered to non-degenerate cubes: a
        // non-degenerate (deg+1)-cube can have a degenerate face, and that face
        // still has to take part in the pairing for the cube to be found. A
        // degenerate face contributes nothing to the column anyway, since
        // `boundary_column` drops it.
        for bi in 0..maps.len() {
            let is_same = ai == bi;

            // The gluing condition is axis independent, so the appended-axis
            // attempt doubles as the check: if it fails, no axis glues, and
            // taking that fast path here keeps the cost of the non-glueable
            // majority of pairs to a single scan.
            if maps[ai]
                .try_combine_along(&maps[bi], dim, is_same)
                .is_err()
            {
                continue;
            }

            // Both orientations are needed: the active cube must be able to sit
            // on either the negative or the positive face.
            for axis in 0..=dim {
                for (neg, pos) in [(ai, bi), (bi, ai)] {
                    let Ok(tau) = maps[neg].try_combine_along(&maps[pos], axis, is_same) else {
                        continue;
                    };
                    if tau.is_degenerate() {
                        continue;
                    }
                    let tau_val = tau.inner().value();
                    if !processed.insert(tau_val) {
                        continue;
                    }
                    enumerated += 1;

                    let col = boundary_column::<CAP, V>(&tau, lower, rows);
                    if let Outcome::Death { row } = reducer.push(tau_val, &col) {
                        // Row indexes `births` (see `row_index` above), and the
                        // reduction loop guarantees the row was not already a
                        // pivot, so each birth dies at most once.
                        still_active.remove(&births[row as usize]);
                    }

                    if active && !still_active.contains(&a_val) {
                        // A has died, so every remaining cube reachable from
                        // this inner loop would need some other active face to
                        // justify enumerating it; those will be reached from
                        // their own outer iteration.
                        continue 'outer;
                    }
                }
            }
        }
    }

    (reducer, enumerated)
}

/// RedZeD DH: Betti numbers `beta_0 ..= beta_n` of `target` over `F_p`.
///
/// With `active` set, degree `n+1` is handled by active enumeration (Section 2);
/// otherwise every `(n+1)`-cube is enumerated. Both compute the same answer, and
/// the corresponding `top_enumerated` counts expose the difference.
///
/// `p` must be prime. Supported degrees are `n = 1, 2, 3`, because the `Poly`
/// capacity dispatch is hardcoded (see `boundary_column`). `n = 3` builds the
/// whole degree-3 level first and then pairs it for degree 4, which is what
/// makes it the first expensive case.
pub fn redzed_dh<V: UGraph>(target: &V, n: u32, p: i64, active: bool) -> RedzedOutcome {
    assert!(
        (1..=3).contains(&n),
        "only degrees 1 through 3 are wired up so far, got {n}"
    );

    let s0 = CubeGraph::new(0u32);
    let s1 = CubeGraph::new(1u32);
    let s2 = CubeGraph::new(2u32);

    let deg0: Vec<CubeMapN<'_, '_, V>> = gen_level(&s0, target);
    let deg1: Vec<CubeMapN<'_, '_, V>> = gen_level(&s1, target);
    let deg2: Vec<CubeMapN<'_, '_, V>> = gen_level(&s2, target);

    let l0 = Level::from_maps(&deg0);
    let l1 = Level::from_maps(&deg1);
    let l2 = Level::from_maps(&deg2);

    // Every degree-0 generator is a birth (d_0 = 0), so the rows of d_1 are C_0.
    let rows0 = row_index(&l0.basis);

    // d_1 : C_1 -> C_0
    let r1 = run_phase1::<2, V>(&deg1, &l0, &rows0, p);

    if n == 1 {
        // Top phase: enumerate 2-cubes by pairing 1-cubes, d_2 : C_2 -> C_1.
        let rows1 = row_index(r1.births());
        let (r2, top_enumerated) = run_phase2::<4, V>(&deg1, &l1, r1.births(), &rows1, p, active);
        return RedzedOutcome {
            betti: vec![
                l0.dim() as i64 - r1.rank() as i64,
                r1.births().len() as i64 - r2.rank() as i64,
            ],
            births: vec![l0.dim(), r1.births().len()],
            rank_d: vec![r1.rank(), r2.rank()],
            top_enumerated,
        };
    }

    // d_2 : C_2 -> C_1, with rows the births of degree 1.
    let rows1 = row_index(r1.births());
    let r2 = run_phase1::<4, V>(&deg2, &l1, &rows1, p);
    let rows2 = row_index(r2.births());

    if n == 2 {
        // Top phase: enumerate 3-cubes by pairing 2-cubes, d_3 : C_3 -> C_2.
        let (r3, top_enumerated) =
            run_phase2::<6, V>(&deg2, &l2, r2.births(), &rows2, p, active);
        return RedzedOutcome {
            betti: vec![
                l0.dim() as i64 - r1.rank() as i64,
                r1.births().len() as i64 - r2.rank() as i64,
                r2.births().len() as i64 - r3.rank() as i64,
            ],
            births: vec![l0.dim(), r1.births().len(), r2.births().len()],
            rank_d: vec![r1.rank(), r2.rank(), r3.rank()],
            top_enumerated,
        };
    }

    // n == 3: Phase 1 continues to degree 3, building the 3-cubes by pairing
    // 2-cubes, then the top phase pairs 3-cubes to reach degree 4.
    let deg3: Vec<CubeMapN<'_, '_, V>> = combined_cube_maps(&deg2).collect();
    let l3 = Level::from_maps(&deg3);
    let r3 = run_phase1::<6, V>(&deg3, &l2, &rows2, p);
    let rows3 = row_index(r3.births());
    let (r4, top_enumerated) = run_phase2::<8, V>(&deg3, &l3, r3.births(), &rows3, p, active);

    RedzedOutcome {
        betti: vec![
            l0.dim() as i64 - r1.rank() as i64,
            r1.births().len() as i64 - r2.rank() as i64,
            r2.births().len() as i64 - r3.rank() as i64,
            r3.births().len() as i64 - r4.rank() as i64,
        ],
        births: vec![
            l0.dim(),
            r1.births().len(),
            r2.births().len(),
            r3.births().len(),
        ],
        rank_d: vec![r1.rank(), r2.rank(), r3.rank(), r4.rank()],
        top_enumerated,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::computations::homology::{betti_numbers, report_to};
    use crate::graphs::extras::{c_n_graph, greene_sphere};

    const P: i64 = 1_000_003;

    #[test]
    fn test_reducer_basic_pairing() {
        // Two parallel edges between two vertices: rank 1, one birth, one death.
        let mut r = ColumnReducer::new(P);
        let e = [(0u32, 1i64), (1, 1)];
        assert_eq!(r.push(10, &e), Outcome::Death { row: 0 });
        assert_eq!(r.push(11, &e), Outcome::Birth);
        assert_eq!(r.rank(), 1);
        assert_eq!(r.births(), &[11]);
        assert_eq!(r.columns_pushed(), 2);
    }

    #[test]
    fn test_reducer_zero_and_empty_columns() {
        let mut r = ColumnReducer::new(P);
        assert_eq!(r.push(1, &[]), Outcome::Birth);
        // Cancels to zero under normalize.
        assert_eq!(r.push(2, &[(0u32, 1i64), (0, -1)]), Outcome::Birth);
        assert_eq!(r.rank(), 0);
        assert_eq!(r.births().len(), 2);
    }

    #[test]
    fn test_reducer_chain_of_deaths() {
        // Columns over rows {0,1,2}:
        //   c1 = [0,1]         -> pivot at row 0
        //   c2 = [1,2]         -> pivot at row 1
        //   c3 = c1 + c2 = [0,2,1] (coeffs 1,2,1) -> needs two subtractions,
        //        then reduces to zero.
        let mut r = ColumnReducer::new(P);
        assert_eq!(r.push(1, &[(0u32, 1i64), (1, 1)]), Outcome::Death { row: 0 });
        assert_eq!(r.push(2, &[(1u32, 1i64), (2, 1)]), Outcome::Death { row: 1 });
        assert_eq!(r.push(3, &[(0u32, 1i64), (1, 2), (2, 1)]), Outcome::Birth);
        assert_eq!(r.rank(), 2);
        assert_eq!(r.births(), &[3]);
    }

    /// RedZeD must reproduce the rank-based computation exactly.
    #[test]
    fn test_redzed_matches_rank_betti() {
        let cases = [
            ("C3", c_n_graph(3)),
            ("C4", c_n_graph(4)),
            ("C5", c_n_graph(5)),
            ("C6", c_n_graph(6)),
            ("Gsph", greene_sphere()),
        ];
        for (name, g) in cases {
            for active in [false, true] {
                let out = redzed_dh(&g, 2, P, active);
                let expected = betti_numbers(&g, P);
                assert_eq!(
                    out.betti, expected,
                    "{name} (active = {active}): RedZeD gave {:?}, rank gave {expected:?}",
                    out.betti
                );
                // beta_p = |births_p| - rank(d_{p+1})
                for p in 0..=2 {
                    assert_eq!(
                        out.betti[p],
                        out.births[p] as i64 - out.rank_d[p] as i64,
                        "{name}: beta_{p} identity"
                    );
                }
            }
        }
    }

    #[test]
    fn test_redzed_degree_one() {
        // beta_0 = 1 and beta_1 = 1 for C5, matching the degree-2 computation.
        let out = redzed_dh(&c_n_graph(5), 1, P, false);
        assert_eq!(out.betti, vec![1, 1]);

        // C3 and C4 have trivial H_1.
        assert_eq!(redzed_dh(&c_n_graph(3), 1, P, false).betti, vec![1, 0]);
        assert_eq!(redzed_dh(&c_n_graph(4), 1, P, false).betti, vec![1, 0]);
    }

    /// Pairing must yield every cube exactly once. `run_phase1` relies on this
    /// when it walks the combined degree-3 list: a generator pushed twice would
    /// reduce to zero the second time and be miscounted as a birth, inflating
    /// `|births_3|` and so `beta_2` and `beta_3`.
    #[test]
    fn test_pair_enumeration_has_no_duplicates() {
        for (name, g) in [
            ("C3", c_n_graph(3)),
            ("C4", c_n_graph(4)),
            ("C5", c_n_graph(5)),
            ("Gsph", greene_sphere()),
        ] {
            let s2 = CubeGraph::new(2u32);
            let deg2: Vec<CubeMapN<'_, '_, _>> = gen_level(&s2, &g);
            let deg3: Vec<CubeMapN<'_, '_, _>> = combined_cube_maps(&deg2).collect();

            let mut values: Vec<u64> = deg3.iter().map(|m| m.inner().value()).collect();
            let total = values.len();
            values.sort_unstable();
            values.dedup();

            assert_eq!(
                values.len(),
                total,
                "{name}: {total} degree-3 cubes but only {} distinct",
                values.len()
            );
        }
    }

    /// The paper's `dim H_2(Gsph) = 1`, via RedZeD rather than via ranks.
    #[test]
    fn test_redzed_greene_sphere_h2() {
        let out = redzed_dh(&greene_sphere(), 2, P, false);
        assert_eq!(out.betti, vec![1, 0, 1]);
    }

    /// Degree 3, cross-checked against `homology::report_to`, which reduces
    /// with the *full* non-degenerate basis as rows and no `r`-map. The two
    /// implementations share only `partial` and the F_p primitives, so
    /// agreement is real evidence.
    ///
    /// Ignored by default: degree 4 is built by pairing degree-3 cubes, so the
    /// cost grows as the square of that level, which a debug build cannot
    /// afford. Run it against a release build:
    ///
    /// `cargo test --release --lib -- --ignored --nocapture degree_three`
    ///
    /// The graph is chosen for a small degree-4 level, which is what actually
    /// sets the cost. C5 has 2475 degree-3 cubes and ~987,735 degree-4 cubes.
    /// Avoid the other obvious candidates: C4 has 256 degree-2 maps of which
    /// only 64 are non-degenerate, so its degree-3 level is an order of
    /// magnitude longer, and C3 is a complete graph, so *every* one of the
    /// 3^16 maps Q_4 -> C3 is a graph map.
    #[test]
    #[ignore = "degree 4 enumeration is far too slow in a debug build; run with --release"]
    fn test_redzed_degree_three_matches_rank() {
        for (name, g) in [("C5", c_n_graph(5))] {
            let expected = report_to(&g, 3, P).betti;
            let full = redzed_dh(&g, 3, P, false);
            let active = redzed_dh(&g, 3, P, true);

            for (mode, out) in [("full", &full), ("active", &active)] {
                assert_eq!(
                    out.betti, expected,
                    "{name} ({mode}): RedZeD gave {:?}, rank gave {expected:?}",
                    out.betti
                );
                for p in 0..=3 {
                    assert_eq!(
                        out.betti[p],
                        out.births[p] as i64 - out.rank_d[p] as i64,
                        "{name} ({mode}): beta_{p} identity"
                    );
                }
            }

            eprintln!(
                "{name}: beta = {:?}, |C3| = {}, top cubes enumerated {} -> {} ({:.1}% skipped)",
                full.betti,
                full.births[3] + full.rank_d[2],
                full.top_enumerated,
                active.top_enumerated,
                100.0 * (full.top_enumerated - active.top_enumerated) as f64
                    / full.top_enumerated.max(1) as f64
            );
        }
    }

    /// Active enumeration must give the same answer while constructing fewer
    /// top-dimensional cubes. `active_skipped` is reported rather than asserted
    /// as a fixed number, since the exact saving is graph dependent.
    #[test]
    fn test_active_enumeration_skips_top_cubes() {
        for (name, g) in [
            ("C4", c_n_graph(4)),
            ("C5", c_n_graph(5)),
            ("C6", c_n_graph(6)),
            ("Gsph", greene_sphere()),
        ] {
            let full = redzed_dh(&g, 2, P, false);
            let active = redzed_dh(&g, 2, P, true);
            assert_eq!(full.betti, active.betti, "{name}: Betti numbers differ");
            assert!(
                active.top_enumerated <= full.top_enumerated,
                "{name}: active enumeration enumerated MORE cubes: {} vs {}",
                active.top_enumerated,
                full.top_enumerated
            );
            eprintln!(
                "{name}: top cubes enumerated {} -> {} ({:.1}% skipped)",
                full.top_enumerated,
                active.top_enumerated,
                100.0 * (full.top_enumerated - active.top_enumerated) as f64
                    / full.top_enumerated.max(1) as f64
            );
        }
    }
}
