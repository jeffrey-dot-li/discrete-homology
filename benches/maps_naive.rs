use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use discrete_homology::graph_maps::generate_maps_naive;
use discrete_homology::prelude::*;

fn bench_my_cpu_bound(c: &mut Criterion) {
    let mut group = c.benchmark_group("high");

    group
        .sample_size(10) // default is 100
        .measurement_time(std::time::Duration::from_secs(5))
        .warm_up_time(std::time::Duration::from_secs(1));

    group.bench_function(BenchmarkId::new("3cube_3cube_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [2.0000 s 2.0242 s 2.0491 s]
            let n = 3;
            let expected_checked = 2u64.pow(n).pow(2u32.pow(n));
            // Assert it fits in usize
            assert!(expected_checked <= u64::MAX);
            use cube::CubeGraph;
            let cube = CubeGraph::new(n);
            let (maps, num_checked) = generate_maps_naive(&cube, &cube);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 15488);
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });
    group.finish();

    let mut group_non_constrained = c.benchmark_group("low");

    group_non_constrained.bench_function(BenchmarkId::new("2cube_2cube_naive", "1e6"), |b| {
        b.iter(|| {
            // time:   [45.306 µs 45.938 µs 46.950 µs]
            let n = 2;
            let expected_checked = 2u64.pow(n).pow(2u32.pow(n));
            // Assert it fits in usize
            assert!(expected_checked <= u64::MAX);
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
            // Assert it fits in usize
            assert!(expected_checked <= u64::MAX);

            let (maps, num_checked) = generate_maps_naive(&cube, &gsphere_graph);
            // print!("num maps: {} num_checked: {}\n", maps.len(), num_checked);
            assert!(num_checked == expected_checked); // there are many maps from cube to itself
            assert!(maps.len() == 442, "num maps was {}", maps.len());
            // black_box prevents the compiler from optimizing away inputs/outputs
            std::hint::black_box((maps, num_checked));
        })
    });
    group_non_constrained.finish();
}

criterion_group!(benches, bench_my_cpu_bound);
criterion_main!(benches);
