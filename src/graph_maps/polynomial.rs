use arrayvec::ArrayVec;
use num_traits::{PrimInt, Unsigned, Zero};
use std::cmp::Ordering;

// R is Z (i32). Basis elements are u64.
// Capacity is strictly limited to 16.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Poly<T: PrimInt + Unsigned + Zero + Default, const CAP: usize> {
    // Invariants:
    // 1. Both vectors must always have the same length.
    // 2. 'indices' is sorted in ascending order.
    indices: ArrayVec<T, CAP>,
    coeffs: ArrayVec<i32, CAP>,
}
impl<T: PrimInt + Unsigned + Default + Zero, const CAP: usize> Poly<T, CAP> {
    /// Create a new element with a single term
    pub fn new(index: T, coeff: i32) -> Self {
        let mut poly = Self::default();
        if coeff != 0 {
            poly.indices.push(index);
            poly.coeffs.push(coeff);
        }
        poly
    }
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn mul_scalar(mut self, scalar: i32) -> Self {
        if scalar == 0 {
            return Self::zero();
        }
        for coeff in self.coeffs.iter_mut() {
            *coeff *= scalar;
        }
        self
    }

    /// Adds two polynomials.
    /// Returns None if the result exceeds the capacity of 16.
    pub fn add(&self, other: &Self) -> Option<Self> {
        let mut result = Self::default();
        let mut i = 0;
        let mut j = 0;

        // Helper to keep pushes in sync and check bounds
        // Using a closure here captures 'result' to reduce boilerplate
        let mut push_term = |idx: T, val: i32| -> bool {
            if result.indices.is_full() {
                return false;
            }
            result.indices.push(idx);
            result.coeffs.push(val);
            true
        };

        // Standard linear merge algorithm
        while (i < self.indices.len()) && (j < other.indices.len()) {
            let idx_lhs = self.indices[i];
            let idx_rhs = other.indices[j];

            match idx_lhs.cmp(&idx_rhs) {
                Ordering::Less => {
                    if !push_term(idx_lhs, self.coeffs[i]) {
                        return None;
                    }
                    i += 1;
                }
                Ordering::Greater => {
                    if !push_term(idx_rhs, other.coeffs[j]) {
                        return None;
                    }
                    j += 1;
                }
                Ordering::Equal => {
                    let sum = self.coeffs[i] + other.coeffs[j];
                    if sum != 0 {
                        if !push_term(idx_lhs, sum) {
                            return None;
                        }
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        // Append remaining elements from Self
        while i < self.indices.len() {
            if !push_term(self.indices[i], self.coeffs[i]) {
                return None;
            }
            i += 1;
        }

        // Append remaining elements from Other
        while j < other.indices.len() {
            if !push_term(other.indices[j], other.coeffs[j]) {
                return None;
            }
            j += 1;
        }

        Some(result)
    }

    /// Returns the number of non-zero terms
    pub fn len(&self) -> usize {
        self.indices.len() // Coeffs.len() is guaranteed to be the same
    }
}
