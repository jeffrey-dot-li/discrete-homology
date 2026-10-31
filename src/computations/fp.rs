//! Arithmetic over the prime field `F_p`, plus the sparse-column primitives
//! shared by the rank computation ([`super::homology`]) and the RedZeD
//! reduction ([`super::redzed`]).
//!
//! A sparse column is a `Vec<(row, coefficient)>`; [`normalize`] puts one into
//! the canonical form the rest of the module assumes: coefficients reduced mod
//! `p`, duplicate rows merged, zero entries dropped, rows ascending.

use std::collections::HashMap;

/// Reduce `x` into `0..p`.
#[inline]
pub fn mod_norm(x: i64, p: i64) -> i64 {
    let m = x % p;
    if m < 0 {
        m + p
    } else {
        m
    }
}

/// Modular inverse of `a` mod the prime `p`, via the extended Euclidean
/// algorithm. `a` must be nonzero mod `p`.
pub fn mod_inv(a: i64, p: i64) -> i64 {
    let a = mod_norm(a, p);
    let (mut t, mut new_t) = (0i64, 1i64);
    let (mut r, mut new_r) = (p, a);
    while new_r != 0 {
        let q = r / new_r;
        let t_next = t - q * new_t;
        t = new_t;
        new_t = t_next;
        let r_next = r - q * new_r;
        r = new_r;
        new_r = r_next;
    }
    debug_assert_eq!(r, 1, "mod_inv requires a prime modulus and a != 0");
    mod_norm(t, p)
}

/// Canonicalise a sparse column: reduce mod `p`, merge duplicate rows, drop
/// zeros, sort ascending by row.
pub fn normalize(col: &[(u32, i64)], p: i64) -> Vec<(u32, i64)> {
    let mut acc: HashMap<u32, i64> = HashMap::new();
    for &(row, c) in col {
        let e = acc.entry(row).or_insert(0);
        *e = mod_norm(*e + mod_norm(c, p), p);
    }
    let mut out: Vec<(u32, i64)> = acc.into_iter().filter(|&(_, c)| c != 0).collect();
    out.sort_unstable_by_key(|&(row, _)| row);
    out
}

/// `a - scale * b`, termwise mod `p`. Both inputs must be in the canonical
/// form produced by [`normalize`]; the output is too. Zero terms are dropped.
pub fn sub_scaled(a: &[(u32, i64)], b: &[(u32, i64)], scale: i64, p: i64) -> Vec<(u32, i64)> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        let (ra, ca) = a[i];
        let rb = b[j].0;
        if ra < rb {
            out.push((ra, ca));
            i += 1;
        } else if ra > rb {
            let cb = mod_norm(-(b[j].1 * scale), p);
            if cb != 0 {
                out.push((rb, cb));
            }
            j += 1;
        } else {
            let v = mod_norm(ca - b[j].1 * scale, p);
            if v != 0 {
                out.push((ra, v));
            }
            i += 1;
            j += 1;
        }
    }
    while i < a.len() {
        out.push(a[i]);
        i += 1;
    }
    while j < b.len() {
        let v = mod_norm(-(b[j].1 * scale), p);
        if v != 0 {
            out.push((b[j].0, v));
        }
        j += 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: i64 = 1_000_003;

    #[test]
    fn test_mod_norm_and_inv() {
        assert_eq!(mod_norm(-1, P), P - 1);
        assert_eq!(mod_norm(P + 3, P), 3);
        assert_eq!(mod_norm(0, P), 0);
        for a in [1i64, 2, 3, 17, P - 1] {
            let inv = mod_inv(a, P);
            assert_eq!(mod_norm(a * inv, P), 1, "inverse of {a}");
        }
    }

    #[test]
    fn test_normalize_merges_and_sorts() {
        // Same row twice, summing to zero -> dropped.
        assert_eq!(normalize(&[(3, 5), (3, -5), (1, 2)], P), vec![(1, 2)]);
        // Same row twice, nonzero sum.
        assert_eq!(normalize(&[(2, 4), (1, 1), (2, 3)], P), vec![(1, 1), (2, 7)]);
        // Coefficients reduced, negatives wrapped.
        assert_eq!(normalize(&[(0, P + 2), (1, -1)], P), vec![(0, 2), (1, P - 1)]);
        // Three copies of a row cancelling out completely.
        assert_eq!(normalize(&[(7, 1), (7, 1), (7, -2)], P), vec![]);
        assert_eq!(normalize(&[], P), vec![]);
    }

    #[test]
    fn test_sub_scaled() {
        let a = normalize(&[(0, 1), (1, 2), (3, 4)], P);
        let b = normalize(&[(0, 1), (3, 1)], P);
        // a - 1*b == 2*(1) + 3*(3)
        assert_eq!(sub_scaled(&a, &b, 1, P), vec![(1, 2), (3, 3)]);
        // Cancelling the leading row is the point of the operation.
        assert_eq!(sub_scaled(&a, &b, 1, P).first().unwrap().0, 1);
        // Subtracting twice the pivot leaves row 0 at -1 mod p, so it stays.
        assert_eq!(
            sub_scaled(&a, &b, 2, P),
            vec![(0, mod_norm(-1, P)), (1, 2), (3, mod_norm(4 - 2, P))]
        );
        // Empty result when the two agree.
        assert_eq!(sub_scaled(&b, &b, 1, P), vec![]);
    }
}
