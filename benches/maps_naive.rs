use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use discrete_homology::graph_maps::cube_maps::{combined_cube_maps, CubeMap};
use discrete_homology::graph_maps::generate_maps_naive;
use discrete_homology::graph_maps::stack_map::generate_maps_naive_stack;
use discrete_homology::prelude::*;

const CUBE3_CUBE3_NUM_MAPS: usize = 15488;
const CUBE3_TO_GSPHERE_NUM_MAPS: usize = 22762;

fn bench_my_cpu_bound(c: &mut Criterion) {
    let mut group = c.benchmark_group("high");

    group
        .sample_size(10) // default is 100
        .measurement_time(std::time::Duration::from_secs(5))
        .warm_up_time(std::time::Duration::from_secs(1));

    group.bench_function(BenchmarkId::new("3cube_3cube_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [2.0000 s 2.0242 s 2.0491 s]
            // 648.20 ms on desktop
            let n = 3;
            let expected_checked = 2u64.pow(n).pow(2u32.pow(n));
            use cube::CubeGraph;
            let cube = CubeGraph::new(n);
            let (maps, num_checked) = generate_maps_naive(&cube, &cube);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == CUBE3_CUBE3_NUM_MAPS);
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });
    group.bench_function(BenchmarkId::new("3cube_gsphere_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [4.2128 s 4.2246 s 4.2355] on desktop
            let n = 3;
            use cube::CubeGraph;
            use extras::greene_sphere;
            let cube = CubeGraph::new(n);
            let gsphere_graph = greene_sphere();

            let expected_checked: u64 = (gsphere_graph.n() as u64).checked_pow(cube.n()).unwrap();

            let (maps, num_checked) = generate_maps_naive(&cube, &gsphere_graph);
            // print!("num maps: {} num_checked: {}\n", maps.len(), num_checked);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(
                maps.len() == CUBE3_TO_GSPHERE_NUM_MAPS,
                "num maps was {}",
                maps.len()
            );
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });
    group.finish();

    let mut group_non_constrained = c.benchmark_group("low");

    group_non_constrained.bench_function(BenchmarkId::new("3cube_c5_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [60.776 ms 61.818 ms 63.138 ms]
            let n = 3;
            use cube::CubeGraph;
            use extras::c_n_graph;
            let cube = CubeGraph::new(n);
            let c5_graph = c_n_graph(5);

            let expected_checked: u64 = (c5_graph.n() as u64).checked_pow(cube.n()).unwrap();

            let (maps, num_checked) = generate_maps_naive(&cube, &c5_graph);
            // print!("num maps: {} num_checked: {}\n", maps.len(), num_checked);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 2475, "num maps was {}", maps.len());
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });

    group_non_constrained.bench_function(BenchmarkId::new("2cube_c5_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [81.899 µs 81.899 µs 81.89 µs]
            let n = 2;
            use cube::CubeGraph;
            use extras::c_n_graph;
            let cube = CubeGraph::new(n);
            let c5_graph = c_n_graph(5);

            let expected_checked: u64 = (c5_graph.n() as u64).checked_pow(cube.n()).unwrap();

            let (maps, num_checked) = generate_maps_naive(&cube, &c5_graph);
            // print!("num maps: {} num_checked: {}\n", maps.len(), num_checked);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 95, "num maps was {}", maps.len());
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });

    group_non_constrained.bench_function(BenchmarkId::new("2cube_2cube_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [45.306 µs 45.938 µs 46.950 µs]
            let n = 2;
            let expected_checked = 2u64.pow(n).pow(2u32.pow(n));
            use cube::CubeGraph;
            let cube = CubeGraph::new(n);
            let (maps, num_checked) = generate_maps_naive(&cube, &cube);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 84);
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });

    group_non_constrained.bench_function(BenchmarkId::new("2cube_gsphere_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [524.13 µs 527.23 µs 532.70 µs]
            let n = 2;
            use cube::CubeGraph;
            use extras::greene_sphere;
            let cube = CubeGraph::new(n);
            let gsphere_graph = greene_sphere();

            let expected_checked: u64 = (gsphere_graph.n() as u64).checked_pow(cube.n()).unwrap();

            let (maps, num_checked) = generate_maps_naive(&cube, &gsphere_graph);
            // print!("num maps: {} num_checked: {}\n", maps.len(), num_checked);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 442, "num maps was {}", maps.len());
            // 442 / 10^4 = ~4%
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });

    group_non_constrained.finish();

    bench_non_naive(c);
}

fn bench_non_naive(c: &mut Criterion) {
    let mut group_high = c.benchmark_group("high_non_naive");
    group_high
        .sample_size(10) // default is 100
        .measurement_time(std::time::Duration::from_secs(5))
        .warm_up_time(std::time::Duration::from_secs(1));
    group_high.bench_function(BenchmarkId::new("3cube_from_2cube_gsphere", "1e6"), |b| {
        let n = 3;

        use cube::CubeGraph;
        let cube2 = CubeGraph::new(n - 1);
        let gsphere = extras::greene_sphere();
        let (cube2_maps, _) = generate_maps_naive(&cube2, &gsphere);
        let cube2_maps = cube2_maps
            .into_iter()
            .map(CubeMap::from)
            .collect::<Vec<_>>();
        // black_box prevents the compiler from optimizing away inputs/outputs
        b.iter(|| {
            // time:  [4.7626 ms 4.7679 ms 4.7732 ms]

            // O(2475^2) checks = O(~6_000_000) vs 10^8 naive = ~6%
            let cube3_maps = combined_cube_maps(&cube2_maps);

            assert!(
                cube3_maps.len() == CUBE3_TO_GSPHERE_NUM_MAPS,
                "num maps was {}",
                cube3_maps.len()
            );
            std::hint::black_box(cube3_maps);
        })
    });

    group_high.bench_function(
        BenchmarkId::new("3cube_from_2cube_gsphere_stackmap", "1e6"),
        |b| {
            let n = 3;

            use cube::CubeGraph;
            let cube2 = CubeGraph::new(n - 1);
            let gsphere = extras::greene_sphere();
            let (cube2_maps, _) = generate_maps_naive_stack(&cube2, &gsphere);
            let cube2_maps = cube2_maps
                .into_iter()
                .map(CubeMap::from)
                .collect::<Vec<_>>();

            b.iter(|| {
                // time:  [4.7626 ms 4.7679 ms 4.7732 ms]
                // black_box prevents the compiler from optimizing away inputs/outputs
                let cube3_maps = combined_cube_maps(&cube2_maps);
                // O(2475^2) checks = O(~6_000_000) vs 10^8 naive = ~6%
                assert!(
                    cube3_maps.len() == CUBE3_TO_GSPHERE_NUM_MAPS,
                    "num maps was {}",
                    cube3_maps.len()
                );
                std::hint::black_box(cube3_maps);
            })
        },
    );

    group_high.finish();
}

criterion_group!(benches, bench_my_cpu_bound);
criterion_main!(benches);
