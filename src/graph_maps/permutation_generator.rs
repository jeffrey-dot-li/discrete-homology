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
pub struct PermutationIterator {
    // Misaligned rip
    val: u64,
    original_val: u64,
    index: u32,
    codomain: u32,
    domain: u32,
}

impl Iterator for PermutationIterator {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.domain {
            return None;
        }
        let remainder = (self.val % self.codomain as u64) as u32;
        self.val /= self.codomain as u64;
        self.index += 1;
        Some(remainder)
    }
}

impl PermutationIterator {
    pub fn new(val: u64, codomain: u32, domain: u32) -> Self {
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
    type Item = PermutationIterator;
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
