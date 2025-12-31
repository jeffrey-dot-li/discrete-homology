use num_traits::{PrimInt, Unsigned};
/// Increment a number represented as digits in a given base.
/// Returns false if overflow occurs (i.e., all digits wrap around to 0).
pub fn increment_mod_base(digits: &mut [u32], base: u32) -> bool {
    for d in digits.iter_mut() {
        *d += 1;
        if *d >= base {
            *d = 0;
        } else {
            return true;
        }
    }
    false
}

#[derive(Clone, Debug, Copy, PartialEq, Eq)]
pub struct PermutationGenerator {
    seed: u64,
    current: Option<u64>,
    num_permutations: u64,
    codomain: u32,
    domain: u32,
}

#[derive(Clone, Debug, Copy, PartialEq, Eq)]
pub struct PermutationIterator<T: PrimInt + Unsigned> {
    // Misaligned rip
    val: T,
    original_val: T,
    index: u32,
    codomain: u32,
    domain: u32,
}

impl<T: PrimInt + Unsigned> Iterator for PermutationIterator<T> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.domain {
            return None;
        }
        let codomain = T::from(self.codomain).unwrap();
        let remainder = (self.val % codomain)
            .to_u32()
            .expect("codomain fits in u32");
        self.val = self.val / codomain;
        self.index += 1;
        Some(remainder)
    }
}

impl<T: PrimInt + Unsigned> PermutationIterator<T> {
    pub fn new(val: T, codomain: u32, domain: u32) -> Self {
        Self {
            original_val: val,
            val,
            index: 0,
            codomain,
            domain,
        }
    }
}

impl Iterator for PermutationGenerator {
    type Item = PermutationIterator<u64>;
    fn next(&mut self) -> Option<Self::Item> {
        // This is so ugly but w/e
        if self.current.map(|c| c == self.seed).unwrap_or(false) {
            return None;
        }
        let current = self.current.unwrap_or(self.seed);
        let iter = PermutationIterator::new(current, self.codomain, self.domain);
        self.current = Some((current + 1) % self.num_permutations);
        Some(iter)
    }
}

impl PermutationGenerator {
    pub fn new(domain: u32, codomain: u32, val: u64) -> Self {
        // Assert codomain^domain fits in u64
        // Equivalent to (log2 codomain) * domain < 64
        let permutations = (codomain as u64).checked_pow(domain).unwrap();

        Self {
            codomain,
            domain,
            seed: val,
            current: None,
            num_permutations: permutations,
        }
    }
}
