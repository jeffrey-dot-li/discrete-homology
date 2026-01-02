use std::borrow::Cow;

use crate::graph_maps::{GraphMap, GraphMapError, VertGraphMap};
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;
use crate::prelude::*;

#[derive(Debug)]
pub struct CubeMap<D: Dim, V: UGraph, M>
where
    M: GraphMap<CubeGraph<D>, V>,
{
    map: M,
    degenerate_indices: Vec<bool>, // -> should just be a u32 bitmask
    _marker: std::marker::PhantomData<(D, V)>,
}

// Inclusion / Forgetful functor
impl<V: UGraph, M: GraphMap<CubeGraph<u32>, V>> From<M> for CubeMap<u32, V, M> {
    fn from(value: M) -> Self {
        let degenerate_indices = (0..value.domain().dim().size() as usize)
            .map(|i| d(&value, i as u32, false).vert_maps == d(&value, i as u32, true).vert_maps)
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
// TODO: Implement d_i as stack graph map
fn d<V: UGraph>(
    map: &impl GraphMap<CubeGraph<u32>, V>,
    i: u32,
    sign: bool,
) -> VertGraphMap<'_, '_, CubeGraph<u32>, V> {
    let dim = map.domain().dim().size();
    assert!(dim != 0u32);
    debug_assert!(i < dim);
    let new_dim = dim.checked_sub(1).unwrap();
    let num_verts = 2_u32.checked_pow(new_dim).unwrap();
    let vert_maps = (0..num_verts).map(|v| map.map(put_bit(v, i, if sign { 1 } else { 0 })));
    unsafe {
        VertGraphMap::new_unchecked(
            Cow::Owned(CubeGraph::new(new_dim)),
            Cow::Borrowed(map.codomain()),
            Cow::Owned(vert_maps.collect()),
        )
    }
}

// impl<D: Dim, V: UGraph, M: GraphMap<CubeGraph<D>, V>> From<CubeMap<D, V, M>> for M {
//     fn from(value: CubeMap<D, V, M>) -> Self {
//         value.map
//     }
// }

impl<D: Dim, V: UGraph, M> CubeMap<D, V, M>
where
    M: GraphMap<CubeGraph<D>, V>,
{
    // Consumes self and returns the inner map
    pub fn into_inner(self) -> M {
        self.map
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
                .change_domain(CubeGraph::new(n + 1), combined_verts)
        };
        let map2 = unsafe {
            other
                .map
                .change_domain(CubeGraph::new(n + 1), combined_verts_2)
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
) -> Vec<CubeMap<u32, V, M>> {
    let mut combined_maps = Vec::new();
    let len = maps.len();
    for i in 0..len {
        let combined = maps[i].try_combine(&maps[i], true).unwrap().0;
        // Generate combined cube maps
        // Then filter for only non-degenerate
        // then each map needs to check what the operator is on it.
        // Then check if that is non zero,
        // then if it is non zero add it to the matrix
        combined_maps.push(combined);

        for j in i + 1..len {
            if let Ok(combined) = maps[i].try_combine(&maps[j], false) {
                combined_maps.push(combined.0);
                combined_maps.push(combined.1);
            }
        }
    }
    combined_maps
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
    use crate::graph_maps::generate_maps_naive;

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
        // black_box prevents the compiler from optimizing away inputs/outputs
        let cube_n_combined_maps = combined_cube_maps(&cube_prev_maps);

        let cube_n_naive_maps = generate_maps_naive(&source, &target).0;

        assert!(
            cube_n_combined_maps.len() == cube_n_naive_maps.len(),
            "num maps combined was {}, but naive was {}",
            cube_n_combined_maps.len(),
            cube_n_naive_maps.len()
        );
    }
}
