use crate::graph_maps::cube_maps::{check_is_degenerate_naive, d};
use crate::prelude::*;
use crate::{
    graph_maps::{
        cube_maps::CubeMap,
        hashset::{Hashable, OpenHashSet, UIntOpenHashSet},
        polynomial::Poly,
        stack_map::StackGraphMap,
        GraphMap,
    },
    graphs::cube::CubeGraph,
};

pub fn partial<'u, 'v, V: UGraph, T: UINT, const CAP: usize, H: OpenHashSet<T>>(
    map: &StackGraphMap<'u, 'v, CubeGraph<u32>, V, T>,
    non_degen_maps: &H,
) -> Poly<T, CAP> {
    let mut sum = Poly::<T, CAP>::zero();
    debug_assert!(CAP == (2 * map.domain().dim().size()) as usize);
    for i in 0..map.domain().dim().size() {
        let sign = if (i % 2) == 0 { 1 } else { -1 };
        let neg_side = d(map, i, false);
        let pos_side = d(map, i, true);
        if non_degen_maps.get(&neg_side.value()) {
            sum = sum.add(&Poly::new(neg_side.value(), -1 * sign)).unwrap();
        }
        if non_degen_maps.get(&pos_side.value()) {
            sum = sum.add(&Poly::new(pos_side.value(), 1 * sign)).unwrap();
        }
    }
    sum
}

pub fn build_hashmap<'u, 'v, 'w, V: UGraph + 'v, T: UINT + 'w, I>(iter: I) -> UIntOpenHashSet<T>
where
    I: IntoIterator<Item = &'w StackGraphMap<'u, 'v, CubeGraph<u32>, V, T>> + 'w,
    'u: 'w,
    'v: 'w,
{
    let mut set = UIntOpenHashSet::with_capacity(1024);
    set.extend(iter.into_iter().map(|m| m.value()));
    set
}

pub fn boundary_map<'u, 'v, V, T: UINT + Hashable, H: OpenHashSet<T>, const CAP: usize>(
    map: &CubeMap<u32, V, StackGraphMap<'u, 'v, CubeGraph<u32>, V, T>>,
    non_degen_set: &'u H,
) -> Option<Poly<T, CAP>>
where
    V: UGraph + 'v,
    'u: 'v,
{
    if map.is_degenerate() {
        return None;
    }
    let boundary = partial(map.inner(), non_degen_set);
    if boundary.is_zero() {
        None
    } else {
        Some(boundary)
    }
}

pub fn boundary_maps<'u, 'v, V, T: UINT + Hashable, H: OpenHashSet<T>, const CAP: usize, Item, I>(
    combined_cube_maps: I,
    non_degen_set: &'u H,
) -> impl Iterator<Item = Poly<T, CAP>> + 'u + 'v
where
    V: UGraph + 'v,
    I: IntoIterator<Item = CubeMap<u32, V, StackGraphMap<'u, 'v, CubeGraph<u32>, V, T>>>,
    I::IntoIter: 'u + 'v,
    'u: 'v,
{
    combined_cube_maps
        .into_iter()
        .filter_map(move |map| boundary_map(&map, non_degen_set))
}

// should this be on graph maps or on CubeMap?
impl<V: UGraph, M: GraphMap<CubeGraph<u32>, V>> CubeMap<u32, V, M> {
    pub fn partial_naive(&self) -> Vec<(i32, M)> {
        let mut maps = Vec::new();
        for i in 0..self.inner().domain().dim().size() {
            // False is negative, True is positive.
            let sign = if (i % 2) == 0 { -1 } else { 1 };
            let neg_side = d(self.inner(), i, false);
            let pos_side = d(self.inner(), i, true);
            if !check_is_degenerate_naive(&neg_side) {
                maps.push((sign, neg_side));
            }
            if !check_is_degenerate_naive(&pos_side) {
                maps.push((-sign, pos_side));
            }
        }
        maps
    }
}

#[cfg(test)]
mod tests {

    use crate::graph_maps::{
        cube_maps::{combined_cube_maps, get_valid_graph_map},
        stack_map::generate_maps_naive_stack,
    };

    use super::*;

    #[test]
    fn test_partial_cube_map() {
        const N: usize = 2;

        use cube::CubeGraph;
        let source = CubeGraph::new(N as u32);
        let target = extras::greene_sphere();
        let (cube_maps, _) = generate_maps_naive_stack(&source, &target);
        let cube_maps = cube_maps.into_iter().map(CubeMap::from).collect::<Vec<_>>();

        let non_degen_set = build_hashmap(
            cube_maps
                .iter()
                .filter(|&map| !map.is_degenerate())
                .map(|m| m.inner()),
        );
        //
        let combined_iter = combined_cube_maps(&cube_maps);

        for combined_map in combined_iter {
            let naive_partial = combined_map.partial_naive();
            let mut poly: Poly<u64, { 2 * (N + 1) }> =
                partial(combined_map.inner(), &non_degen_set);
            let original_poly = poly.clone();
            for (sign, element) in naive_partial.clone() {
                let poly_component = Poly::new(element.value(), sign);
                poly = poly.add(&poly_component.mul_scalar(-1)).unwrap();
            }
            assert!(
                poly == Poly::zero(),
                "Partial from naive {:?} does not match partial from hashset {:?} {:?}",
                original_poly,
                naive_partial,
                poly,
            );
        }
    }

    #[test]
    fn test_d_i_cube() {
        use arbtest::arbtest;
        arbtest(|u| {
            let dim = 2;
            let source = CubeGraph::new(dim);
            let target = extras::greene_sphere();

            let map = get_valid_graph_map(&source, &target, u.arbitrary()?);
            let dn_map_pos = CubeMap::from(d(&map, dim - 1, true));
            let dn_map_neg = CubeMap::from(d(&map, dim - 1, false));
            let recombined_map = dn_map_neg.try_combine(&dn_map_pos, false);

            if recombined_map.is_err() {
                panic!(
                    "Failed to recombine maps: {:?} {:?} {:?}",
                    recombined_map.err().unwrap(),
                    dn_map_neg.inner().vert_maps(),
                    dn_map_pos.inner().vert_maps(),
                );
            }
            let recombined_map = recombined_map.unwrap().0;

            assert!(
                recombined_map.inner().vert_maps() == map.vert_maps(),
                "Recombined map does not match original map {:?} vs {:?}",
                recombined_map.inner().vert_maps(),
                map.vert_maps()
            );
            Ok(())
        });
    }
}
