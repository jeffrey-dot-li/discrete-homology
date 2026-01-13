use std::borrow::Cow;

use crate::graph_maps::hashset::{OpenHashSet, UIntOpenHashSet};
use crate::graph_maps::polynomial::Poly;
use crate::graph_maps::stack_map::{StackGraphMap, UINT};
use crate::graph_maps::{GraphMap, GraphMapError, VertGraphMap};
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;
use crate::prelude::*;

use std::fmt::Debug;

#[derive(Debug)]
pub struct CubeMap<D: Dim, V: UGraph, M>
where
    M: GraphMap<CubeGraph<D>, V>,
{
    map: M,
    // TODO: Maybe use stack map here
    degenerate_indices: Vec<bool>, // -> should just be a u32 bitmask
    _marker: std::marker::PhantomData<(D, V)>,
}

// Inclusion / Forgetful functor
impl<V: UGraph, M: GraphMap<CubeGraph<u32>, V>> From<M> for CubeMap<u32, V, M> {
    fn from(value: M) -> Self {
        let degenerate_indices = (0..value.domain().dim().size() as usize)
            .map(|i| d(&value, i as u32, false) == d(&value, i as u32, true))
            .collect();

        Self {
            map: value,
            degenerate_indices,
            _marker: std::marker::PhantomData,
        }
    }
}

fn put_bit(x: u32, pos: u32, value: u32) -> u32 {
    debug_assert!(pos < 32);
    debug_assert!(value == 0 || value == 1);

    let lower_mask = (1u32 << pos) - 1;
    let lower = x & lower_mask;

    let upper = x & !lower_mask;
    let upper_shifted = upper << 1;

    lower | (value << pos) | upper_shifted
}

fn d<V: UGraph, M: GraphMap<CubeGraph<u32>, V>>(map: &M, i: u32, sign: bool) -> M {
    let dim = map.domain().dim().size();
    assert!(dim != 0u32);
    debug_assert!(i < dim);
    let new_dim = dim.checked_sub(1).unwrap();
    let num_verts = 2_u32.checked_pow(new_dim).unwrap();
    let vert_maps = (0..num_verts).map(|v| map.map(put_bit(v, i, if sign { 1 } else { 0 })));
    unsafe { map.change_domain(CubeGraph::new(new_dim), vert_maps) }
}

fn check_is_degenerate_naive<V: UGraph, M: GraphMap<CubeGraph<u32>, V>>(map: &M) -> bool {
    (0..map.domain().dim().size()).any(|i| d(map, i, false) == d(map, i, true))
}

// should this be on graph maps or on CubeMap?
impl<V: UGraph, M: GraphMap<CubeGraph<u32>, V>> CubeMap<u32, V, M> {
    pub fn partial_naive(&self) -> Vec<(i32, M)> {
        let mut maps = Vec::new();
        for i in 0..self.map.domain().dim().size() {
            // False is negative, True is positive.
            let sign = if (i % 2) == 0 { -1 } else { 1 };
            let neg_side = d(&self.map, i, false);
            let pos_side = d(&self.map, i, true);
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

impl<D: Dim, V: UGraph, M> CubeMap<D, V, M>
where
    M: GraphMap<CubeGraph<D>, V>,
{
    // Consumes self and returns the inner map
    pub fn into_inner(self) -> M {
        self.map
    }
    pub fn is_degenerate(&self) -> bool {
        self.degenerate_indices.iter().any(|b| *b)
    }
}

impl<V: UGraph, M: GraphMap<CubeGraph<u32>, V>> CubeMap<u32, V, M> {
    pub fn dim(&self) -> u32 {
        self.map.domain().dim()
    }

    pub fn try_combine<'w>(
        &self,
        other: &CubeMap<u32, V, M>,
        is_same: bool,
    ) -> Result<(CubeMap<u32, V, M>, CubeMap<u32, V, M>), GraphMapError> {
        use GraphMapError as E;
        assert!(self.map.codomain() == other.map.codomain());
        let n = self.map.domain().n();
        let dim = self.map.domain().dim();
        for (i, (x, y)) in self
            .map
            .mapped_vertices()
            .zip(other.map.mapped_vertices())
            .enumerate()
        {
            if !self.map.codomain().is_edge(x, y) {
                return Err(E::BadEdge(i as u32, i as u32 + n, x, y));
            }
        }

        let combined_verts = self
            .map
            .mapped_vertices()
            .chain(other.map.mapped_vertices());
        let combined_verts_2 = other
            .map
            .mapped_vertices()
            .chain(self.map.mapped_vertices());

        let map = unsafe {
            self.map
                .change_domain(CubeGraph::new(dim + 1), combined_verts)
        };
        let map2 = unsafe {
            other
                .map
                .change_domain(CubeGraph::new(dim + 1), combined_verts_2)
        };
        let degenerate_indices_iter = self
            .degenerate_indices
            .iter()
            .zip(other.degenerate_indices.iter())
            .map(|(a, b)| *a && *b);
        let mut degenerate_indices = Vec::with_capacity((n + 1) as usize);
        degenerate_indices.extend(degenerate_indices_iter);
        degenerate_indices.push(is_same);

        Ok((
            Self {
                map,
                degenerate_indices: degenerate_indices.clone(),
                _marker: std::marker::PhantomData,
            },
            Self {
                map: map2,
                degenerate_indices,
                _marker: std::marker::PhantomData,
            },
        ))
    }
}

pub fn combined_cube_maps<V: UGraph, M: GraphMap<CubeGraph<u32>, V>>(
    maps: &[CubeMap<u32, V, M>],
) -> impl Iterator<Item = CubeMap<u32, V, M>> + '_ {
    (0..maps.len()).flat_map(move |i| {
        // Self-combine: maps[i] with itself
        let self_combined = std::iter::once(maps[i].try_combine(&maps[i], true).unwrap().0);

        // Cross-combine: maps[i] with maps[j] for all j > i
        let cross_combines = (i + 1..maps.len())
            .filter_map(move |j| maps[i].try_combine(&maps[j], false).ok())
            .flat_map(|(combined_0, combined_1)| [combined_0, combined_1]);

        self_combined.chain(cross_combines)
    })
}

use crate::graph_maps::permutation_generator::PermutationGenerator;
pub fn get_valid_graph_map<'u, 'v, U: UGraph, V: UGraph>(
    source: &'u U,
    target: &'v V,
    seed: u64,
) -> VertGraphMap<'u, 'v, U, V> {
    // TODO: Write proper generator for valid graph maps
    assert!(target.n() > 0);
    let mut generator = PermutationGenerator::new(source.n(), target.n(), seed);
    for _ in 0..1_000_000 {
        let next_iter = generator.next().unwrap();
        let candidate_map: Result<VertGraphMap<'u, 'v, U, V>, GraphMapError> =
            VertGraphMap::try_from(
                Cow::Borrowed(source),
                Cow::Borrowed(target),
                next_iter,
                &mut vec![0; source.n() as usize],
            );
        if let Ok(map) = candidate_map {
            return map;
        }
        // else {
        //     return unsafe {
        //         VertGraphMap::new_unchecked(
        //             Cow::Borrowed(source),
        //             Cow::Borrowed(target),
        //             Cow::Owned(generator.current.iter().copied().collect::<Vec<_>>()),
        //         )
        //     };
        // }
    }

    panic!(
        "Exceeded maximum iterations searching for valid graph map from {} to {}",
        source.n(),
        target.n()
    );
}

#[cfg(test)]
mod tests {
    use crate::{
        graph_maps::{generate_maps_naive, stack_map::generate_maps_naive_stack},
        graphs::extras::greene_sphere,
    };

    use super::*;

    #[test]
    fn test_put_bit() {
        // put_bit inserts a bit at position `pos`, shifting higher bits left

        // Insert 0 at position 0 of 0b101 (5) -> 0b1010 (10)
        assert_eq!(
            put_bit(0b101, 0, 0),
            0b1010,
            "Inserting 0 at position 0 should shift all bits left"
        );

        // Insert 1 at position 0 of 0b101 (5) -> 0b1011 (11)
        assert_eq!(
            put_bit(0b101, 0, 1),
            0b1011,
            "Inserting 1 at position 0 should add bit at position 0"
        );

        // Insert 0 at position 1 of 0b101 (5) -> 0b1001 (9)
        // Original: _ 1 0 1
        // Result:   1 0 0 1
        assert_eq!(
            put_bit(0b101, 1, 0),
            0b1001,
            "Inserting 0 at position 1 should preserve lower bit, insert 0, shift upper bits"
        );

        // Insert 1 at position 1 of 0b101 (5) -> 0b1011 (11)
        // Original: _ 1 0 1
        // Result:   1 0 1 1
        assert_eq!(
            put_bit(0b101, 1, 1),
            0b1011,
            "Inserting 1 at position 1 should preserve lower bit, insert 1, shift upper bits"
        );

        // Insert 1 at position 2 of 0b11 (3) -> 0b111 (7)
        // Original: _ _ 1 1
        // Result:   _ 1 1 1
        assert_eq!(
            put_bit(0b11, 2, 1),
            0b111,
            "Inserting 1 at position 2 should preserve lower 2 bits and add bit at position 2"
        );

        // Insert 0 at position 0 of 0 -> 0
        assert_eq!(put_bit(0, 0, 0), 0, "Inserting 0 into 0 should give 0");

        // Insert 1 at position 0 of 0 -> 1
        assert_eq!(put_bit(0, 0, 1), 1, "Inserting 1 into 0 should give 1");

        // Insert at higher positions
        // Insert 1 at position 3 of 0b111 (7) -> 0b1111 (15)
        assert_eq!(
            put_bit(0b111, 3, 1),
            0b1111,
            "Inserting 1 at position 3 of 0b111 should give 0b1111"
        );
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
                    dn_map_neg.map.vert_maps,
                    dn_map_pos.map.vert_maps,
                );
            }
            let recombined_map = recombined_map.unwrap().0;

            assert!(
                recombined_map.map.vert_maps == map.vert_maps,
                "Recombined map does not match original map {:?} vs {:?}",
                recombined_map.map.vert_maps,
                map.vert_maps
            );
            Ok(())
        });
    }

    fn test_num_graph_maps<D: UGraph, C: UGraph>(domain: &D, codomain: &C, expected_num: usize) {
        let (cube_maps, _) = generate_maps_naive(domain, codomain);

        assert!(
            cube_maps.len() == expected_num,
            "num maps was {}, expected {}",
            cube_maps.len(),
            expected_num
        );
    }

    #[test]
    fn test_simple_graph_maps() {
        use cube::CubeGraph;
        let cube2 = CubeGraph::new(2);
        let cube3 = CubeGraph::new(3);
        let c5 = extras::c_n_graph(5);
        test_num_graph_maps(&cube2, &c5, 95);
        test_num_graph_maps(&cube2, &cube3, 320);
        test_num_graph_maps(&cube2, &greene_sphere(), 442);
    }

    #[test]
    fn test_2cube_gsphere_combined() {
        let n = 2;

        use cube::CubeGraph;
        let source = CubeGraph::new(n);
        let cube_prev = CubeGraph::new(n - 1);
        let target = extras::greene_sphere();
        // let target = extras::c_n_graph(5);
        let (cube_prev_maps, _) = generate_maps_naive(&cube_prev, &target);
        let cube_prev_maps = cube_prev_maps
            .into_iter()
            .map(CubeMap::from)
            .collect::<Vec<_>>();
        // Collect the iterator to count results
        let cube_n_combined_maps: Vec<_> = combined_cube_maps(&cube_prev_maps).collect();

        let cube_n_naive_maps = generate_maps_naive(&source, &target).0;

        assert!(
            cube_n_combined_maps.len() == cube_n_naive_maps.len(),
            "num maps combined was {}, but naive was {}",
            cube_n_combined_maps.len(),
            cube_n_naive_maps.len()
        );
    }

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
                .map(|m| &m.map),
        );
        let combined_iter = combined_cube_maps(&cube_maps);

        for combined_map in combined_iter {
            let naive_partial = combined_map.partial_naive();
            let mut poly: Poly<u64, { 2 * (N + 1) }> = partial(&combined_map.map, &non_degen_set);
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
}
