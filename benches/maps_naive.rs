use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn my_cpu_bound(x: &[u64]) -> u64 {
    // your code here
    x.iter()
        .fold(0u64, |acc, &v| acc.wrapping_add(v.rotate_left(13)))
}

fn bench_my_cpu_bound(c: &mut Criterion) {
    let data: Vec<u64> = (0..1_000_000).collect();

    c.bench_function("my_cpu_bound 1e6", |b| {
        b.iter(|| {
            // black_box prevents the compiler from optimizing away inputs/outputs
            let out = my_cpu_bound(black_box(&data));
            black_box(out);
        })
    });
}

criterion_group!(benches, bench_my_cpu_bound);
criterion_main!(benches);
