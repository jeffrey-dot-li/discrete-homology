use super::uint::UINT;
use std::fmt;

pub trait Hashable {
    fn hash(&self) -> u64;
}

/// A simple open-addressing HashSet for `u128` using linear probing.
/// - No deletions (so we don't need tombstones).
/// - Power-of-two capacity for fast indexing.
/// - Separate occupancy bitmap because all `u128` values are valid keys.
///
/// This is optimized for *membership / insert* workloads.
#[derive(Clone)]
pub struct UIntOpenHashSet<T: UINT + Hashable> {
    keys: Vec<T>,
    // TODO: Get rid of occupied, use 0 for "empty", then search i + 1.
    occupied: Vec<u8>, // 0 = empty, 1 = full
    len: usize,
    mask: usize, // capacity - 1
    max_load: f64,
}

pub trait OpenHashSet<K> {
    fn get(&self, key: &K) -> bool;
}

impl<T: UINT + Hashable> OpenHashSet<T> for UIntOpenHashSet<T> {
    fn get(&self, key: &T) -> bool {
        self.contains(*key)
    }
}
impl<T: UINT + Hashable> fmt::Debug for UIntOpenHashSet<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("UIntOpenHashSet")
            .field("len", &self.len)
            .field("capacity", &self.capacity())
            .field("load_factor", &self.load_factor())
            .finish()
    }
}

impl<T: UINT + Hashable> UIntOpenHashSet<T> {
    /// Create with a minimum capacity (will be rounded up to power-of-two, at least 16).
    pub fn with_capacity(min_capacity: usize) -> Self {
        let cap: usize = next_pow2(min_capacity.max(16));
        let keys: Vec<T> = vec![T::zero(); cap];
        let occupied: Vec<u8> = vec![0u8; cap];
        Self {
            keys,
            occupied,
            len: 0,
            mask: cap - 1,
            max_load: 0.60,
        }
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.keys.len()
    }

    #[inline]
    pub fn load_factor(&self) -> f64 {
        self.len as f64 / self.capacity() as f64
    }

    /// Returns true if inserted (was not already present).
    pub fn insert(&mut self, key: T) -> bool {
        if (self.len + 1) as f64 > (self.capacity() as f64 * self.max_load) {
            self.rehash(self.capacity() * 2);
        }

        let mut idx: usize = self.index_for(key);
        loop {
            if self.occupied[idx] == 0 {
                self.occupied[idx] = 1;
                self.keys[idx] = key;
                self.len += 1;
                return true;
            }

            if self.keys[idx] == key {
                return false;
            }

            idx = (idx + 1) & self.mask;
        }
    }

    /// Returns true if present.
    pub fn contains(&self, key: T) -> bool {
        let mut idx: usize = self.index_for(key);
        loop {
            if self.occupied[idx] == 0 {
                return false;
            }
            if self.keys[idx] == key {
                return true;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    /// Insert many keys.
    pub fn extend<I: IntoIterator<Item = T>>(&mut self, iter: I) {
        for k in iter {
            let _ = self.insert(k);
        }
    }

    /// Rebuild into a new table with `new_capacity` (rounded up to power-of-two).
    fn rehash(&mut self, new_capacity: usize) {
        let cap: usize = next_pow2(new_capacity.max(16));
        let mut new_set = UIntOpenHashSet {
            keys: vec![T::zero(); cap],
            occupied: vec![0u8; cap],
            len: 0,
            mask: cap - 1,
            max_load: self.max_load,
        };

        for i in 0..self.capacity() {
            if self.occupied[i] != 0 {
                let _ = new_set.insert(self.keys[i]);
            }
        }

        *self = new_set;
    }

    #[inline]
    fn index_for(&self, key: T) -> usize {
        // Hash u128 -> u64, then take low bits via mask.
        let h: u64 = Hashable::hash(&key);
        (h as usize) & self.mask
    }
}

impl Default for UIntOpenHashSet<u128> {
    fn default() -> Self {
        Self::with_capacity(16)
    }
}

impl<T: UINT> Hashable for T {
    #[inline]
    fn hash(&self) -> u64 {
        hash_u128_to_u64(self.to_u128().unwrap())
    }
}

/// Hash a u128 into a u64 using a strong 64-bit mixer (SplitMix64),
/// with folding of high/low halves.
///
/// This is fast and has good low-bit behavior for power-of-two tables.
#[inline]
fn hash_u128_to_u64(x: u128) -> u64 {
    let lo: u64 = x as u64;
    let hi: u64 = (x >> 64) as u64;

    // Two mixes with distinct perturbations, then combine.
    let a: u64 = splitmix64(lo ^ 0x9E37_79B9_7F4A_7C15u64);
    let b: u64 = splitmix64(hi ^ 0xBF58_476D_1CE4_E5B9u64);

    // Combine (xor is fine; could also use wrapping_add/rotate).
    a ^ b
}

#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15u64);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9u64);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EBu64);
    x ^ (x >> 31)
}

#[inline]
fn next_pow2(n: usize) -> usize {
    // Next power of two for usize, with n>=1.
    // For n already power of two, returns n.
    n.next_power_of_two()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_insert_contains() {
        let mut s: UIntOpenHashSet<u128> = UIntOpenHashSet::default();
        assert!(!s.contains(123u128));
        assert!(s.insert(123u128));
        assert!(s.contains(123u128));
        assert!(!s.insert(123u128));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn many_keys() {
        let mut s: UIntOpenHashSet<u128> = UIntOpenHashSet::with_capacity(1_000);
        let base: u128 = (1u128 << 100) + 7;
        for i in 0..50_000u128 {
            let k: u128 = base.wrapping_add(i.wrapping_mul(1_000_003u128));
            assert!(s.insert(k));
        }
        for i in 0..50_000u128 {
            let k: u128 = base.wrapping_add(i.wrapping_mul(1_000_003u128));
            assert!(s.contains(k));
        }
        assert!(!s.contains(base.wrapping_add(999_999_999u128)));
    }

    #[test]
    fn rehash_preserves() {
        let mut s: UIntOpenHashSet<u128> = UIntOpenHashSet::with_capacity(16);
        let keys: Vec<u128> = (0..10_000u128).map(|i| (i << 64) ^ (i * 17)).collect();
        for &k in &keys {
            let _ = s.insert(k);
        }
        for &k in &keys {
            assert!(s.contains(k));
        }
    }
}
