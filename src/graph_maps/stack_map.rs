use crate::graph_maps::permutation_generator::PermutationIterator;
use crate::graph_maps::{GraphMapError, VertGraphMap};
use crate::graphs::UGraph;
use crate::prelude::*;
use num_traits::{PrimInt, Unsigned};
use std::borrow::Cow;

#[derive(Debug, Clone)]
pub struct StackGraphMap<'u, 'v, U: UGraph, V: UGraph, T>
where
    T: PrimInt + Unsigned,
{
    domain: Cow<'u, U>,
    codomain: Cow<'v, V>,
    vert_maps: T,
}

pub fn check_fits<T: PrimInt + Unsigned>(domain_n: u32, codomain_n: u32) -> bool {
    // Check if codomain_n^domain_n fits in T
    // We need: domain_n * log2(codomain_n) <= bits_in_T

    let bits_in_t = (std::mem::size_of::<T>() * 8) as u32;

    // Number of bits needed to represent values [0, codomain_n)
    let bits_per_digit = if codomain_n == 0 {
        return true; // Empty codomain, any domain fits
    } else {
        32 - (codomain_n - 1).leading_zeros()
    };

    // Check if domain_n * bits_per_digit <= bits_in_t
    // Use saturating_mul to avoid overflow in the check itself
    bits_per_digit.saturating_mul(domain_n) <= bits_in_t
}

impl<'u, 'v, U: UGraph, V: UGraph, T: PrimInt + Unsigned> StackGraphMap<'u, 'v, U, V, T> {
    pub fn iter(&self) -> PermutationIterator<T> {
        self.into_iter()
    }
    /// # Safety
    ///
    /// Yeah this doesn't actually check if vert_maps are legit
    pub unsafe fn new_unchecked(
        domain: impl Into<Cow<'u, U>>,
        codomain: impl Into<Cow<'v, V>>,
        vert_maps: T,
    ) -> Self {
        Self {
            domain: domain.into(),
            codomain: codomain.into(),
            vert_maps,
        }
    }

    pub fn try_from(
        domain: impl Into<Cow<'u, U>>,
        codomain: impl Into<Cow<'v, V>>,
        mapped_verts: impl IntoIterator<Item = u32>,
        // Workspace is a length n mutable array
        workspace: &mut [u32],
    ) -> Result<Self, GraphMapError> {
        let domain = domain.into();
        let codomain = codomain.into();
        let n = domain.n() as usize;
        debug_assert!(check_fits::<T>(domain.n(), codomain.n()));
        debug_assert!(workspace.len() >= n);

        use GraphMapError as E;
        let mut max = 0;

        for (i, mapped_i) in mapped_verts.into_iter().enumerate() {
            assert!(
                mapped_i < codomain.n(),
                "invalid map: codomain node {mapped_i} out of range {:?}",
                codomain.n()
            );
            let neighbors = domain.neighbors(i as u32);
            for neighbor in neighbors {
                let neighbor = neighbor as usize;
                if neighbor < i {
                    let mapped_neighbor = workspace[neighbor];
                    if !codomain.is_edge(mapped_neighbor, mapped_i) {
                        return Err(E::BadEdge(
                            i as u32,
                            neighbor as u32,
                            mapped_i,
                            mapped_neighbor,
                        ));
                    }
                }
            }
            // We have confirmed that edges j - i where j <= i are good
            workspace[i] = mapped_i;
            max = i;
        }
        assert!(max == (n - 1), "max should be n - 1, got {max}");
        Ok(Self {
            vert_maps: Self::slice_convert(workspace, codomain.n()),
            codomain,
            domain,
        })
    }

    pub fn slice_convert(slice: &[u32], codomain_n: u32) -> T {
        debug_assert!(check_fits::<T>(slice.len() as u32, codomain_n));
        let codomain_n = T::from(codomain_n).unwrap();
        slice
            .iter()
            .copied()
            .enumerate()
            .fold::<T, _>(T::zero(), |s, (i, v)| {
                s + (T::from(v).unwrap() * codomain_n.pow(i as u32))
            })
    }
}

impl<'u, 'v, 'a, U: UGraph, V: UGraph, T: PrimInt + Unsigned> IntoIterator
    for &'a StackGraphMap<'u, 'v, U, V, T>
{
    type Item = u32;
    type IntoIter = PermutationIterator<T>;
    fn into_iter(self) -> Self::IntoIter {
        PermutationIterator::new(self.vert_maps, self.codomain.n(), self.domain.n())
    }
}

impl<'u, 'v, U: UGraph, V: UGraph, T: PrimInt + Unsigned> From<StackGraphMap<'u, 'v, U, V, T>>
    for VertGraphMap<'u, 'v, U, V>
{
    fn from(value: StackGraphMap<'u, 'v, U, V, T>) -> Self {
        let vert_maps = (&value).into_iter().collect::<Vec<_>>();
        unsafe { VertGraphMap::new_unchecked(value.domain, value.codomain, Cow::Owned(vert_maps)) }
    }
}

impl<'u, 'v, U: UGraph, V: UGraph, T: PrimInt + Unsigned> From<VertGraphMap<'u, 'v, U, V>>
    for StackGraphMap<'u, 'v, U, V, T>
{
    fn from(value: VertGraphMap<'u, 'v, U, V>) -> Self {
        let val = Self::slice_convert(value.vert_maps.as_slice(), value.codomain.n());
        unsafe { StackGraphMap::new_unchecked(value.domain, value.codomain, val) }
    }
}

impl<'u, 'v, 'a, U: UGraph, V: UGraph, T: PrimInt + Unsigned> From<&'a VertGraphMap<'u, 'v, U, V>>
    for StackGraphMap<'u, 'v, U, V, T>
{
    fn from(value: &VertGraphMap<'u, 'v, U, V>) -> Self {
        let val = Self::slice_convert(value.vert_maps.as_slice(), value.codomain.n());
        unsafe { StackGraphMap::new_unchecked(value.domain.clone(), value.codomain.clone(), val) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_maps::cube_maps::get_valid_graph_map;
    use crate::graphs::cube::{CubeGraph, Newable};
    use crate::graphs::extras;

    #[test]
    fn test_vert_to_stack_to_vert_u32() {
        use arbtest::arbtest;
        arbtest(|u| {
            let dim = 2;
            let source = CubeGraph::new(dim);
            let target = extras::greene_sphere();

            // Check that the conversion fits in u32
            if !check_fits::<u32>(source.n(), target.n()) {
                return Ok(()); // Skip if doesn't fit
            }

            let original_map = get_valid_graph_map(&source, &target, u.arbitrary()?);

            // Convert VertGraphMap -> StackGraphMap -> VertGraphMap
            let stack_map: StackGraphMap<'_, '_, _, _, u32> = (&original_map).into();
            let roundtrip_map: VertGraphMap<'_, '_, _, _> = stack_map.into();

            assert_eq!(
                original_map.vert_maps, roundtrip_map.vert_maps,
                "Roundtrip conversion failed: original {:?} != roundtrip {:?}",
                original_map.vert_maps, roundtrip_map.vert_maps
            );

            Ok(())
        });
    }

    #[test]
    fn test_vert_to_stack_to_vert_u64() {
        use arbtest::arbtest;
        arbtest(|u| {
            let dim = 2;
            let source = CubeGraph::new(dim);
            let target = extras::greene_sphere();

            // u64 should always fit for these small graphs
            assert!(check_fits::<u64>(source.n(), target.n()));

            let original_map = get_valid_graph_map(&source, &target, u.arbitrary()?);

            // Convert VertGraphMap -> StackGraphMap -> VertGraphMap
            let stack_map: StackGraphMap<'_, '_, _, _, u64> = (&original_map).into();
            let roundtrip_map: VertGraphMap<'_, '_, _, _> = stack_map.into();

            assert_eq!(
                original_map.vert_maps, roundtrip_map.vert_maps,
                "Roundtrip conversion failed: original {:?} != roundtrip {:?}",
                original_map.vert_maps, roundtrip_map.vert_maps
            );

            Ok(())
        });
    }

    #[test]
    fn test_check_fits() {
        // Test check_fits with various parameters

        // u32 has 32 bits
        // codomain_n=2, domain_n=32 needs 32 bits (1 bit per vertex) - should fit
        assert!(check_fits::<u32>(32, 2));

        // codomain_n=2, domain_n=33 needs 33 bits - should NOT fit
        assert!(!check_fits::<u32>(33, 2));

        // codomain_n=4 needs 2 bits per vertex, domain_n=16 needs 32 bits - should fit
        assert!(check_fits::<u32>(16, 4));

        // codomain_n=4, domain_n=17 needs 34 bits - should NOT fit
        assert!(!check_fits::<u32>(17, 4));

        // Need this for greene_sphere 5-cube vertices
        assert!(check_fits::<u128>(10, 32));

        // u64 has 64 bits
        // codomain_n=2, domain_n=64 needs 64 bits - should fit
        assert!(check_fits::<u64>(64, 2));

        // codomain_n=2, domain_n=65 needs 65 bits - should NOT fit
        assert!(!check_fits::<u64>(65, 2));

        // Edge cases
        assert!(check_fits::<u32>(0, 0));
        assert!(check_fits::<u32>(0, 100));
        assert!(check_fits::<u32>(100, 0));

        // u8 tests
        assert!(check_fits::<u8>(8, 2)); // 8 bits for binary
        assert!(!check_fits::<u8>(9, 2)); // 9 bits doesn't fit
        assert!(check_fits::<u8>(4, 4)); // 4 vertices, 4 possible values = 2 bits each = 8 bits
        assert!(!check_fits::<u8>(5, 4)); // 5 vertices would need 10 bits

        // u16 tests
        assert!(check_fits::<u16>(16, 2));
        assert!(!check_fits::<u16>(17, 2));
    }
}
