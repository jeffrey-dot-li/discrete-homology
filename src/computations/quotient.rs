//! The hyperoctahedral quotient of the discrete cubical chain complex.
//!
//! Greene-Welker-Wille (arXiv:2508.14639), as used in Ebel-Kapulkin-Kershaw
//! (arXiv:2608.04262, Theorem 1.3): inside `C_•(G)` sits the subcomplex `R_•(G)`
//! generated in degree `n` by the elements `A - sgn(s)*(A o s)` for
//! `A` a non-degenerate `n`-cube and `s` in `Aut(Q_n) = S_n |x (Z_2)^n`, and for
//! a prime `p` not dividing `(n+1)!`,
//!
//! ```text
//! H_n(G; F_p) ~= H_n(C_•(G)/R_•(G); F_p).
//! ```
//!
//! Two consequences drive the implementation.
//!
//! * The quotient identifies `A ~ sgn(s)*(A o s)`, so every element of an orbit
//!   is `+-` every other. A non-semi-degenerate orbit therefore contributes
//!   exactly one dimension, and the basis of `C_n/R_n` is indexed by orbits
//!   rather than by cubes -- a reduction of up to `n! * 2^n`.
//! * An `n`-cube is **semi-degenerate** when `A = A o s` for some `s` of
//!   negative sign. Then `A = -A`, so `A = 0` in odd characteristic and it is
//!   not a generator at all. Every degenerate cube is semi-degenerate (take `s`
//!   the reflection in a coordinate the cube ignores), so this *subsumes* the
//!   degeneracy filter the unquotiented code uses.
//!
//! Boundaries are signed at the level of classes: each face `A o d_i^+-` is
//! replaced by its representative and its coefficient multiplied by the sign of
//! that face within its class; semi-degenerate faces contribute nothing.
//!
//! [`quotient_report`] is deliberately written to mirror
//! [`crate::computations::homology::report_to`], which computes the same Betti
//! numbers without the quotient and with no `r`-map. The two are cross-checked
//! in tests, which is the point of keeping the unquotiented path around.

use std::collections::HashMap;

use crate::computations::homology::{gen_level, rank_mod_p, CubeMapN, Report};
use crate::graph_maps::cube_maps::{combined_cube_maps, d};
use crate::graph_maps::GraphMap;
use crate::graph_maps::polynomial::Poly;
use crate::graph_maps::stack_map::StackGraphMap;
use crate::graph_maps::symmetry::{AutQn, MAX_DIM};
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;

/// A cube map of degree `n`, with the vertex tuple packed base-`|V|`.
///
/// `u64` is enough while `|V|^(2^n) < 2^64`: for the Greene sphere that is
/// degree <= 4, and degree 5 needs `u128` (see `stack_map::check_fits`).
pub type QuotientCube<'u, 'v, V> = StackGraphMap<'u, 'v, CubeGraph<u32>, V, u64>;

/// The class of an `n`-cube in the quotient: the orbit representative it is set
/// equal to and the sign `epsilon` with `A = epsilon * rep`.
///
/// Returns sign `0` for a semi-degenerate cube, whose class vanishes in odd
/// characteristic.
pub fn rep_and_sign<'u, 'v, V: UGraph>(
    a: &QuotientCube<'u, 'v, V>,
    n: u32,
) -> (u64, i32) {
    assert!(n <= MAX_DIM, "n = {n} is out of range");
    let a_value = a.value();
    let domain_n = a.domain().n();
    let codomain_n = a.codomain().n();
    let digits: Vec<u32> = a.iter().collect();
    let base = u64::from(codomain_n);

    // The representative is the orbit element with the smallest packed value.
    // Any total order works; this one is a single integer comparison, which is
    // the whole reason the tuple is packed.
    let mut best: Option<(u64, i32)> = None;
    for k in 0..AutQn::count(n) {
        let s = AutQn::nth(n, k);
        let sign = s.sign(n);
        let perm = s.permutation(n);

        // Packed value of `A o s` without building the map: reindex the digits
        // and fold them back down in base `codomain_n`.
        let mut value = 0u64;
        let mut power = 1u64;
        for x in 0..domain_n {
            let source = AutQn::image_with(&perm, s.reflection(), x) as usize;
            value += u64::from(digits[source]) * power;
            power *= base;
        }

        if sign < 0 && value == a_value {
            // Fixed by a negative element, so `A = -A`.
            return (a_value, 0);
        }
        if best.is_none_or(|(b, _)| value < b) {
            best = Some((value, sign));
        }
    }

    // `rep = A o s` for the winning `s`, and `A o s = sgn(s) * A`, so
    // `A = sgn(s) * rep`.
    best.expect("the group is non-empty, so the orbit is non-empty")
}

/// One degree of the quotient complex: the orbit representatives that form the
/// basis, plus a lookup from any cube of that degree to its class.
pub struct QuotientLevel {
    /// Distinct representatives, ascending. Position = row index.
    pub basis: Vec<u64>,
    /// Representative -> row index.
    pub index: HashMap<u64, u32>,
    /// Every cube seen in this degree -> `Some((rep, sign))`, or `None` when it
    /// is semi-degenerate.
    entries: HashMap<u64, Option<(u64, i32)>>,
}

impl QuotientLevel {
    pub fn empty() -> Self {
        Self {
            basis: Vec::new(),
            index: HashMap::new(),
            entries: HashMap::new(),
        }
    }

    /// `dim C_n/R_n`, the number of non-semi-degenerate orbits.
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// The class of a cube of this degree.
    ///
    /// Panics if the cube was never enumerated, which would mean this degree
    /// was built incompletely -- a silent `None` here would quietly drop a
    /// boundary term.
    pub fn class_of(&self, value: u64) -> Option<(u64, i32)> {
        *self
            .entries
            .get(&value)
            .unwrap_or_else(|| panic!("cube {value} was not enumerated in this degree"))
    }
}

/// The boundary of a generator of `C_n/R_n`, as an element of `C_{n-1}/R_{n-1}`.
///
/// Same alternating sum as the unquotiented differential, except that each face
/// is replaced by its representative with its class sign folded into the
/// coefficient, and semi-degenerate faces are dropped.
pub fn quotient_boundary<'u, 'v, const CAP: usize, V: UGraph>(
    map: &QuotientCube<'u, 'v, V>,
    lower: &QuotientLevel,
    n: u32,
) -> Poly<u64, CAP> {
    debug_assert_eq!(CAP, 2 * n as usize, "capacity must be the number of faces");
    let mut sum = Poly::<u64, CAP>::zero();
    for i in 0..n {
        let sign = if i % 2 == 0 { 1 } else { -1 };
        for (positive, side) in [(false, -sign), (true, sign)] {
            let face = d(map, i, positive);
            if let Some((rep, face_sign)) = lower.class_of(face.value()) {
                sum = sum
                    .add(&Poly::new(rep, side * face_sign))
                    .expect("at most 2n distinct faces fit in CAP = 2n");
            }
        }
    }
    sum
}

/// Build the quotient level for degree `n` and `rank(d_n)` on it, in a single
/// pass over `maps`.
///
/// The columns are taken only from cubes that are their own representative, so
/// there is exactly one column per orbit rather than one per cube: this is where
/// the quotient pays off, and it is why the rank sees a matrix a factor of up
/// to `n! * 2^n` smaller.
fn quotient_level_and_rank<const CAP: usize, V: UGraph>(
    maps: &[CubeMapN<'_, '_, V>],
    lower: &QuotientLevel,
    n: u32,
    p: i64,
) -> (QuotientLevel, usize) {
    let mut entries: HashMap<u64, Option<(u64, i32)>> = HashMap::new();
    let mut basis: Vec<u64> = Vec::new();
    let mut columns: Vec<Vec<(u32, i64)>> = Vec::new();

    for m in maps {
        let inner = m.inner();
        let value = inner.value();
        if entries.contains_key(&value) {
            continue;
        }
        let (rep, sign) = rep_and_sign(inner, n);
        if sign == 0 {
            entries.insert(value, None);
            continue;
        }
        entries.insert(value, Some((rep, sign)));

        if value == rep {
            basis.push(value);
            columns.push(
                quotient_boundary::<CAP, V>(inner, lower, n)
                    .terms()
                    .map(|(v, c)| (lower.index[&v], c as i64))
                    .collect(),
            );
        }
    }

    basis.sort_unstable();
    let index: HashMap<u64, u32> = basis
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i as u32))
        .collect();
    let rank = rank_mod_p(&columns, p);

    (QuotientLevel { basis, index, entries }, rank)
}

/// Betti numbers of the quotient complex, mirroring
/// [`crate::computations::homology::report_to`].
///
/// `p` must be an odd prime with `p > (max_degree + 1)!`, so that Theorem 1.3
/// applies in the degrees computed.
pub fn quotient_report<V: UGraph>(target: &V, max_degree: u32, p: i64) -> Report {
    assert!(
        (1..=3).contains(&max_degree),
        "supported degrees are 1..=3, got {max_degree}"
    );
    let factorial: i64 = (1..=max_degree as i64 + 1).product();
    assert!(
        p > factorial,
        "Theorem 1.3 needs a prime p not dividing (n+1)! = {factorial}; take p > that"
    );

    let s0 = CubeGraph::new(0u32);
    let s1 = CubeGraph::new(1u32);
    let s2 = CubeGraph::new(2u32);

    // The enumeration machinery is shared with the unquotiented path; the
    // quotient only replaces how cubes are turned into generators and rows.
    let deg0: Vec<CubeMapN<'_, '_, V>> = gen_level(&s0, target);
    let deg1: Vec<CubeMapN<'_, '_, V>> = gen_level(&s1, target);
    let deg2: Vec<CubeMapN<'_, '_, V>> = gen_level(&s2, target);

    let (q0, _rank_d0) = quotient_level_and_rank::<0, V>(&deg0, &QuotientLevel::empty(), 0, p);
    let (q1, r1) = quotient_level_and_rank::<2, V>(&deg1, &q0, 1, p);
    let (q2, r2) = quotient_level_and_rank::<4, V>(&deg2, &q1, 2, p);

    let mut c = vec![q0.dim(), q1.dim(), q2.dim()];
    let mut rank_d = vec![r1, r2];

    if max_degree >= 2 {
        let deg3: Vec<CubeMapN<'_, '_, V>> = combined_cube_maps(&deg2).collect();
        let (q3, r3) = quotient_level_and_rank::<6, V>(&deg3, &q2, 3, p);
        rank_d.push(r3);

        if max_degree == 3 {
            let deg4: Vec<CubeMapN<'_, '_, V>> = combined_cube_maps(&deg3).collect();
            let (_q4, r4) = quotient_level_and_rank::<8, V>(&deg4, &q3, 4, p);
            c.push(q3.dim());
            rank_d.push(r4);
        }
    }

    c.truncate(max_degree as usize + 1);

    // beta_p = |C_p/R_p| - rank(d_p) - rank(d_{p+1}); the d_0 term is zero.
    let betti: Vec<i64> = (0..=max_degree as usize)
        .map(|k| {
            let rank_dk = if k == 0 { 0 } else { rank_d[k - 1] as i64 };
            c[k] as i64 - rank_dk - rank_d[k] as i64
        })
        .collect();

    Report { c, rank_d, betti }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::computations::homology::report_to;
    use crate::graph_maps::symmetry::compose;
    use crate::graphs::extras::{c_n_graph, greene_sphere};

    const P: i64 = 1_000_003;

    /// The class of a cube must not depend on which orbit element we start from,
    /// and the sign must move exactly as the group element does. This is the
    /// property that makes the quotient well defined, and it is what would break
    /// first if the sign convention were wrong.
    #[test]
    fn test_rep_and_sign_is_orbit_invariant() {
        let n = 2;
        let s2 = CubeGraph::new(n);
        let target = greene_sphere();
        let deg2: Vec<CubeMapN<'_, '_, _>> = gen_level(&s2, &target);

        let (mut checked, mut skipped) = (0usize, 0usize);
        for m in &deg2 {
            let inner = m.inner();
            let (rep, sign) = rep_and_sign(inner, n);
            if sign == 0 {
                skipped += 1;
                continue;
            }
            checked += 1;
            for k in 0..AutQn::count(n) {
                let s = AutQn::nth(n, k);
                let rotated = compose(inner, &s, n);
                let (rot_rep, rot_sign) = rep_and_sign(&rotated, n);
                assert_eq!(rot_rep, rep, "representative changed under {s:?}");
                assert_eq!(
                    rot_sign,
                    sign * s.sign(n),
                    "sign did not follow sgn(s) for {s:?}"
                );
            }
        }
        // Guard against the whole test going vacuous, which it would if the
        // semi-degeneracy test ever started rejecting everything.
        assert!(checked > 0, "no non-semi-degenerate cubes: this test proves nothing");
        assert!(skipped > 0, "expected some semi-degenerate cubes to be skipped");
    }

    /// Degenerate cubes must come out semi-degenerate, and there must be
    /// non-semi-degenerate ones too (otherwise the test proves nothing).
    #[test]
    fn test_degenerate_implies_semi_degenerate() {
        let n = 2;
        let target = greene_sphere();
        let s2 = CubeGraph::new(n);
        let deg2: Vec<CubeMapN<'_, '_, _>> = gen_level(&s2, &target);

        let (mut degenerate, mut semi) = (0usize, 0usize);
        for m in &deg2 {
            let (_, sign) = rep_and_sign(m.inner(), n);
            if sign == 0 {
                semi += 1;
            }
            if m.is_degenerate() {
                degenerate += 1;
                assert_eq!(sign, 0, "a degenerate cube must have vanishing class");
            }
        }
        assert!(degenerate > 0, "no degenerate cubes to check");
        assert!(
            semi >= degenerate,
            "semi-degenerate ({semi}) should subsume degenerate ({degenerate})"
        );
    }

    /// `d^2 = 0` on the *quotient* complex, which is the decisive check on the
    /// signed boundaries: it fails unless every face is replaced by its
    /// representative with the right sign.
    #[test]
    fn test_quotient_boundary_squared_zero() {
        let target = greene_sphere();
        let s0 = CubeGraph::new(0u32);
        let s1 = CubeGraph::new(1u32);
        let s2 = CubeGraph::new(2u32);

        let deg0: Vec<CubeMapN<'_, '_, _>> = gen_level(&s0, &target);
        let deg1: Vec<CubeMapN<'_, '_, _>> = gen_level(&s1, &target);
        let deg2: Vec<CubeMapN<'_, '_, _>> = gen_level(&s2, &target);
        let deg3: Vec<CubeMapN<'_, '_, _>> = combined_cube_maps(&deg2).collect();

        let (q0, _) = quotient_level_and_rank::<0, _>(&deg0, &QuotientLevel::empty(), 0, P);
        let (q1, _) = quotient_level_and_rank::<2, _>(&deg1, &q0, 1, P);
        let (q2, _) = quotient_level_and_rank::<4, _>(&deg2, &q1, 2, P);

        // A map realising each degree-2 representative, to apply d_2 to.
        let witness2: HashMap<u64, &CubeMapN<'_, '_, _>> = deg2
            .iter()
            .filter(|m| {
                let (rep, sign) = rep_and_sign(m.inner(), 2);
                sign != 0 && rep == m.inner().value()
            })
            .map(|m| (m.inner().value(), m))
            .collect();

        let mut checked = 0usize;
        for a in &deg3 {
            if rep_and_sign(a.inner(), 3).1 == 0 {
                continue;
            }
            let b3: Poly<u64, 6> = quotient_boundary(a.inner(), &q2, 3);
            let mut acc: HashMap<u64, i64> = HashMap::new();
            for (rep, c3) in b3.terms() {
                let m2 = witness2[&rep];
                let b2: Poly<u64, 4> = quotient_boundary(m2.inner(), &q1, 2);
                for (rep1, c2) in b2.terms() {
                    *acc.entry(rep1).or_insert(0) += i64::from(c3) * i64::from(c2);
                }
            }
            assert!(
                acc.values().all(|&c| c == 0),
                "d^2 != 0 on quotient generator {}: {acc:?}",
                a.inner().value()
            );
            checked += 1;
        }
        assert!(checked > 0, "no quotient generators of degree 3 were checked");
    }

    /// The whole point: the quotient must produce the same Betti numbers as the
    /// unquotiented complex, while working with far fewer generators.
    #[test]
    fn test_quotient_matches_unquotiented() {
        let cases = [
            ("C3", c_n_graph(3)),
            ("C4", c_n_graph(4)),
            ("C5", c_n_graph(5)),
            ("C6", c_n_graph(6)),
            ("Gsph", greene_sphere()),
        ];
        for (name, g) in cases {
            let quotient = quotient_report(&g, 2, P);
            let plain = report_to(&g, 2, P);
            assert_eq!(
                quotient.betti, plain.betti,
                "{name}: quotient gave {:?}, unquotiented gave {:?}",
                quotient.betti, plain.betti
            );
        }
    }

    /// The quotient can only shrink the complex, and it can shrink it by *more*
    /// than the group order. There is no bound of the form
    /// `|C_n| >= |C_n/R_n| * n! * 2^n`, for two separate reasons: a
    /// semi-degenerate cube vanishes *entirely* rather than joining an orbit, and
    /// an orbit can be smaller than the group when the cube has a stabiliser
    /// (necessarily of positive sign, or the cube would be semi-degenerate).
    ///
    /// It can also shrink a whole degree to nothing, and does. Every 2-cube of
    /// `C5` is semi-degenerate: the gluing constraints force `b = c` or `a = d`
    /// (a 5-cycle has no triangles, so adjacent `b, c` have common neighbourhood
    /// exactly `{b, c}`, and non-adjacent ones have a unique common neighbour),
    /// and each of those is realised by a negative element of `Aut(Q_2)`. So the
    /// 5-cycle has no degree-2 generators at all, yet still gets `H_1 = F_p` from
    /// the 1-cube level alone.
    ///
    /// That makes the `C5` case vacuous in degree 2, so the Greene sphere is
    /// pinned down too: it keeps a non-empty degree-2 level, where
    /// `beta_2 = 8 - rank(d_2) = 1` is a real computation. A regression in the
    /// semi-degeneracy test breaks one assertion or the other.
    #[test]
    fn test_quotient_shrinks_the_complex() {
        let c5_plain = report_to(&c_n_graph(5), 2, P);
        let c5_quotient = quotient_report(&c_n_graph(5), 2, P);
        let gsph_plain = report_to(&greene_sphere(), 2, P);
        let gsph_quotient = quotient_report(&greene_sphere(), 2, P);

        for (name, plain, quotient) in [
            ("C5", &c5_plain, &c5_quotient),
            ("Gsph", &gsph_plain, &gsph_quotient),
        ] {
            for k in 0..=2usize {
                assert!(
                    quotient.c[k] <= plain.c[k],
                    "{name}: degree {k} grew under the quotient: {} -> {}",
                    plain.c[k],
                    quotient.c[k]
                );
            }
            eprintln!(
                "{name}: |C_n| {:?} -> quotient {:?} (regardless of the group orders 2, 8, 48)",
                plain.c, quotient.c
            );
        }

        assert_eq!(
            c5_quotient.c[2], 0,
            "every 2-cube of C5 should be semi-degenerate"
        );
        assert!(
            gsph_quotient.c[2] > 0,
            "the Greene sphere should keep degree-2 generators, otherwise its beta_2 match is vacuous"
        );
    }

    /// Degree 3, and therefore degree 4: ignored by default because pairing the
    /// degree-3 level is expensive even though the quotient keeps the resulting
    /// matrix small.
    ///
    /// `cargo test --release --lib -- --ignored --nocapture quotient_degree_three`
    #[test]
    #[ignore = "degree 4 enumeration is far too slow in a debug build; run with --release"]
    fn test_quotient_degree_three_matches_unquotiented() {
        let c5 = c_n_graph(5);
        let quotient = quotient_report(&c5, 3, P);
        let plain = report_to(&c5, 3, P);
        eprintln!("C5: quotient |C_n| = {:?}, betti = {:?}", quotient.c, quotient.betti);
        eprintln!("C5: plain    |C_n| = {:?}, betti = {:?}", plain.c, plain.betti);
        assert_eq!(quotient.betti, plain.betti, "beta(C5) mismatch");
    }

    #[test]
    fn test_all_even_characteristic_is_rejected() {
        // The quotient argument needs A = -A to force A = 0, which fails in
        // characteristic 2.
        let result = std::panic::catch_unwind(|| quotient_report(&c_n_graph(3), 2, 2));
        assert!(result.is_err(), "p = 2 should be rejected");
    }
}
