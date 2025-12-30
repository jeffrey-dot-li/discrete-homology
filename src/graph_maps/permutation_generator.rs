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

pub struct PermutationGenerator {
    codomain: u32,
    pub current: Vec<u32>,
    start: Vec<u32>,
}

impl PermutationGenerator {
    pub fn next(&mut self) -> Option<()> {
        increment_mod_base(&mut self.current, self.codomain);
        if self.current == self.start && !self.current.is_empty() {
            return None;
        }
        Some(())
    }

    pub fn new(domain: u32, codomain: u32, val: u64) -> Self {
        // Assert codomain^domain fits in u64
        // Equivalent to (log2 codomain) * domain < 64
        assert!(
            (codomain as u64).checked_pow(domain).is_some(),
            "{codomain}^{domain} must fit in u64"
        );

        let mut val = val;
        let mut start = Vec::with_capacity(domain as usize);
        for _ in 0..domain {
            let remainder = (val % codomain as u64) as u32;
            val /= codomain as u64;
            start.push(remainder);
        }
        Self {
            codomain,
            current: start.clone(),
            start,
        }
    }
}

// #[test]
// fn all_numbers_are_even() {
//     use arbtest::arbtest;

//     arbtest(|u| {
//         let number: u32 = u.arbitrary()?;
//         assert!(number % 2 == 0);
//         Ok(())
//     });
// }
