use std::borrow::Cow;

use crate::graph_maps::cube_maps::CubeMap;
use crate::graph_maps::{GraphMap, GraphMapError, VertGraphMap};
use crate::graphs::cube::{CubeGraph, Newable};
use crate::graphs::UGraph;
use crate::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CubeIso<D: Dim> {
    graph: CubeGraph<D>,
    // Symmetric group S_dim has size dim!
    // This works for dim! < 2^32 i.e. dim <= 12.
    // If use u64, then it works for dim <= 20.
    permutation: u32,

    // Reflection group has size 2^dim
    // The first (dim) bits of this integer indicate whether we are reflecting
    // in that dimension
    reflection: u32,
}

#[derive(Clone, Debug, Copy, PartialEq, Eq)]
pub struct PermutationIter {
    // Current permutation remainder
    p: u64,
    // Current factorial divisor ((k-1)!)
    current_factorial: u64,
    // Bitmask: 0 = available, 1 = used.
    // We use u32 because n <= 20.
    used_mask: u32,
    // Total size of the permutation
    n: usize,
    // How many numbers are left to emit
    remaining: usize,
}

impl CubeIso<u32> {
    pub fn id(n: u32) -> Self {
        Self {
            graph: CubeGraph::new(n),
            permutation: 0,
            reflection: 0,
        }
    }

    pub fn d(&self, i: u32, sign: bool) -> CubeMap<'_, '_, u32, CubeGraph<u32>> {
        assert!(self.graph.dim() != 0);
        // We need to get the index of d_i+- for each i.

        CubeMap::from(
            // unsafe { VertGraphMap::new_unchecked(domain, codomain, vert_maps)}
            {
                VertGraphMap::try_from(
                    Cow::Owned(CubeGraph::new(
                        self.graph.dim().checked_sub(self.graph.dim()).unwrap(),
                    )),
                    Cow::Owned(self.graph),
                    [0u32].iter().copied(),
                    &mut [0u32],
                )
                .unwrap()
            },
        )
    }
}

impl From<CubeIso<u32>> for VertGraphMap<'_, '_, CubeGraph<u32>, CubeGraph<u32>> {
    fn from(value: CubeIso<u32>) -> Self {
        Self {
            domain: Cow::Owned(value.graph),
            codomain: Cow::Owned(value.graph),
            vert_maps: vec![],
        }
    }
}
