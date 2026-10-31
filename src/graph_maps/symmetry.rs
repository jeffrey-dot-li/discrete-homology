//! The hyperoctahedral group acting on the discrete cube.
//!
//! `Aut(Q_n) = S_n |x (Z_2)^n`, of order `n! * 2^n`: a permutation of the
//! coordinates together with a choice of which coordinates to flip. Writing an
//! element `s = (sigma, r)`, its sign is
//!
//! ```text
//! sgn(s) = sgn(sigma) * prod_i (-1)^{r_i} = sgn(sigma) * (-1)^{|r|},
//! ```
//!
//! and it acts on the vertices of `Q_n` by
//!
//! ```text
//! (s . x)_{sigma(i)} = x_i xor r_{sigma(i)},
//! ```
//!
//! that is, the new bit in position `sigma(i)` is the old bit `i` flipped
//! according to `r`. On cube maps the action is `(A o s)(x) = A(s . x)`.
//!
//! The parametrisation is a bijection onto `Aut(Q_n)` and `sgn` is a
//! homomorphism for it: `|a xor b| = |a| + |b| (mod 2)`, so the reflection part
//! is multiplicative, and the permutation part is multiplicative because `sgn`
//! is a class function. `test_action_is_closed_and_sign_multiplicative` pins
//! both down.
//!
//! Only `n <= 12` is supported, since the permutation is stored as a factoradic
//! index and `12!` is the largest factorial below `2^32`.
//!
//! Reference: Ebel, Kapulkin, Kershaw, *Discrete homology computations by
//! reduction to zero differentials*, arXiv:2608.04262, Section 1 (the
//! hyperoctahedral quotient, after Greene-Welker-Wille, arXiv:2508.14639).

use arrayvec::ArrayVec;

use crate::graph_maps::GraphMap;
use crate::graphs::UGraph;

/// Largest dimension supported: `n!` must fit in a `u32` factoradic index.
pub const MAX_DIM: u32 = 12;

const fn factorial(n: u32) -> u64 {
    let mut acc = 1u64;
    let mut i = 2u32;
    while i <= n {
        acc *= i as u64;
        i += 1;
    }
    acc
}

/// An element of `Aut(Q_n) = S_n |x (Z_2)^n`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AutQn {
    /// Factoradic index of the coordinate permutation.
    perm: u32,
    /// Bit `i` set means coordinate `i` is flipped.
    refl: u32,
}

impl AutQn {
    pub fn id() -> Self {
        Self { perm: 0, refl: 0 }
    }

    /// The order of the group, `n! * 2^n`.
    pub fn count(n: u32) -> u64 {
        assert!(n <= MAX_DIM, "Aut(Q_{n}) is out of range, n <= {MAX_DIM}");
        factorial(n) * (1u64 << n)
    }

    /// The `k`-th element, for `k < Self::count(n)`.
    pub fn nth(n: u32, k: u64) -> Self {
        assert!(k < Self::count(n), "index {k} out of range for Aut(Q_{n})");
        Self {
            refl: (k % (1u64 << n)) as u32,
            perm: (k >> n) as u32,
        }
    }

    /// Every element, in index order.
    pub fn all(n: u32) -> impl Iterator<Item = Self> {
        (0..Self::count(n)).map(move |k| Self::nth(n, k))
    }

    /// The reflection bitmask `r`.
    pub fn reflection(&self) -> u32 {
        self.refl
    }

    /// The coordinate permutation `sigma`, with `sigma(i)` at index `i`.
    pub fn permutation(&self, n: u32) -> ArrayVec<u8, 12> {
        let mut available: ArrayVec<u8, 12> = (0..n as u8).collect();
        let mut out = ArrayVec::new();
        let mut code = self.perm as u64;
        for i in (0..n).rev() {
            let f = factorial(i);
            let digit = (code / f) as usize;
            code %= f;
            out.push(available.remove(digit));
        }
        out
    }

    /// The vertex `s . v`, given an already decoded permutation. Decoding once
    /// and reusing it matters: this runs once per vertex per group element.
    #[inline]
    pub fn image_with(perm: &[u8], refl: u32, v: u32) -> u32 {
        let mut out = 0u32;
        for (i, &target) in perm.iter().enumerate() {
            let bit = ((v >> i) & 1) ^ ((refl >> target) & 1);
            out |= bit << target;
        }
        out
    }

    /// The vertex `s . v`.
    #[inline]
    pub fn vertex_image(&self, n: u32, v: u32) -> u32 {
        Self::image_with(&self.permutation(n), self.refl, v)
    }

    /// `sgn(sigma) * (-1)^{|r|}`, as `+1` or `-1`.
    pub fn sign(&self, n: u32) -> i32 {
        let perm = self.permutation(n);
        let mut inversions = 0u32;
        for i in 0..perm.len() {
            for j in (i + 1)..perm.len() {
                if perm[i] > perm[j] {
                    inversions += 1;
                }
            }
        }
        if (inversions + self.refl.count_ones()).is_multiple_of(2) {
            1
        } else {
            -1
        }
    }
}

/// `A o s`, that is `(A o s)(x) = A(s . x)`.
///
/// The domain is unchanged and only the vertex tuple moves. The result is again
/// a cube map, because `s` is an automorphism of `Q_n`, so the unchecked
/// `change_domain` below is sound.
pub fn compose<U: UGraph, V: UGraph, M: GraphMap<U, V>>(a: &M, s: &AutQn, n: u32) -> M {
    let perm = s.permutation(n);
    let verts = (0..(1u32 << n)).map(|x| a.map(AutQn::image_with(&perm, s.reflection(), x)));
    unsafe { a.change_domain(a.domain().clone(), verts) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_order() {
        assert_eq!(AutQn::count(0), 1);
        assert_eq!(AutQn::count(1), 2);
        assert_eq!(AutQn::count(2), 8);
        assert_eq!(AutQn::count(3), 48);
        assert_eq!(AutQn::count(4), 384);
    }

    #[test]
    fn test_identity() {
        let id = AutQn::id();
        assert_eq!(id.sign(0), 1);
        assert_eq!(id.sign(3), 1);
        for v in 0..8u32 {
            assert_eq!(id.vertex_image(3, v), v);
        }
    }

    /// Every element must be a bijection of the vertices, and the elements must
    /// be pairwise distinct as maps -- otherwise the group would be
    /// under-counted and orbits would come out wrong.
    #[test]
    fn test_elements_are_distinct_bijections() {
        for n in 0..=3u32 {
            let verts = 1u32 << n;
            let mut seen: Vec<Vec<u32>> = Vec::new();
            for s in AutQn::all(n) {
                let images: Vec<u32> = (0..verts).map(|v| s.vertex_image(n, v)).collect();
                let mut sorted = images.clone();
                sorted.sort_unstable();
                sorted.dedup();
                assert_eq!(sorted.len(), verts as usize, "n = {n}: not injective");
                seen.push(images);
            }
            assert_eq!(seen.len() as u64, AutQn::count(n));
            seen.sort();
            let before = seen.len();
            seen.dedup();
            assert_eq!(seen.len(), before, "n = {n}: elements are not distinct");
        }
    }

    /// `sign` is a surjective homomorphism onto `{+1, -1}`, so it must split the
    /// group evenly. `n = 0` is excluded: the group is trivial there, so there is
    /// no negative element to find.
    #[test]
    fn test_sign_is_balanced() {
        for n in 1..=4u32 {
            let plus = AutQn::all(n).filter(|s| s.sign(n) == 1).count() as u64;
            assert_eq!(plus, AutQn::count(n) / 2, "n = {n}");
        }
    }

    #[test]
    fn test_single_reflections() {
        for n in 1..=3u32 {
            for axis in 0..n {
                let s = AutQn::nth(n, 1u64 << axis);
                assert_eq!(s.sign(n), -1, "n = {n}, axis = {axis}");
                assert_eq!(
                    s.permutation(n).as_slice(),
                    &(0..n as u8).collect::<Vec<_>>()[..],
                    "a reflection must not permute"
                );
                for v in 0..(1u32 << n) {
                    assert_eq!(s.vertex_image(n, v), v ^ (1 << axis));
                }
            }
        }
    }

    /// The action must be closed, and the element realising a composite action
    /// must have the product of the signs -- which is exactly the multiplicativ-
    /// ity the quotient theory needs.
    #[test]
    fn test_action_is_closed_and_sign_multiplicative() {
        for n in 0..=3u32 {
            let all: Vec<AutQn> = AutQn::all(n).collect();
            let verts = 1u32 << n;
            for s in &all {
                for t in &all {
                    let composed: Vec<u32> = (0..verts)
                        .map(|v| s.vertex_image(n, t.vertex_image(n, v)))
                        .collect();
                    let u = all
                        .iter()
                        .find(|u| (0..verts).all(|v| u.vertex_image(n, v) == composed[v as usize]))
                        .expect("the action must be closed");
                    assert_eq!(
                        u.sign(n),
                        s.sign(n) * t.sign(n),
                        "n = {n}: sign is not multiplicative for {s:?}, {t:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_compose_identity_is_identity() {
        use crate::graph_maps::stack_map::generate_maps_naive_stack;
        use crate::graphs::cube::{CubeGraph, Newable};
        use crate::graphs::extras::c_n_graph;

        let source = CubeGraph::new(1u32);
        let target = c_n_graph(5);
        let (maps, _) = generate_maps_naive_stack(&source, &target);
        for m in &maps {
            assert_eq!(compose(m, &AutQn::id(), 1).value(), m.value());
        }
    }

    /// The orbit of any cube must have size dividing the group order, which
    /// fails if the action is not a genuine action.
    #[test]
    fn test_orbit_size_divides_group_order() {
        use crate::graph_maps::stack_map::generate_maps_naive_stack;
        use crate::graphs::cube::{CubeGraph, Newable};
        use crate::graphs::extras::c_n_graph;

        let n = 2;
        let source = CubeGraph::new(n);
        let target = c_n_graph(5);
        let (maps, _) = generate_maps_naive_stack(&source, &target);
        let order = AutQn::count(n);
        for m in &maps {
            let mut orbit: Vec<u64> = AutQn::all(n).map(|s| compose(m, &s, n).value()).collect();
            orbit.sort_unstable();
            orbit.dedup();
            assert_eq!(
                order % orbit.len() as u64,
                0,
                "orbit of size {} does not divide {order}",
                orbit.len()
            );
        }
    }
}
