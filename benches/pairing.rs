use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::Field;
use halo2curves::bn256::*;
use halo2curves::pairing::Engine;
use halo2curves::pluto_eris::*;
use rand_core::OsRng;
use std::ops::MulAssign;

pub fn pairing_pluto(c: &mut Criterion) {
    let mut g1 = halo2curves::pluto_eris::G1::generator();
    let mut g2 = halo2curves::pluto_eris::G2::generator();
    let mut rng = OsRng;
    let a = halo2curves::pluto_eris::Fq::random(&mut rng);
    let b = halo2curves::pluto_eris::Fq::random(&mut rng);
    g1.mul_assign(a);
    g2.mul_assign(b);

    c.bench_function("pairing on Pluto", move |b| {
        b.iter(|| {
            Pluto::pairing(
                black_box(&halo2curves::pluto_eris::G1Affine::from(g1)),
                black_box(&halo2curves::pluto_eris::G2Affine::from(g2)),
            )
        })
    });
}

pub fn pairing_bn256(c: &mut Criterion) {
    let mut g1 = halo2curves::bn256::G1::generator();
    let mut g2 = halo2curves::bn256::G2::generator();
    let mut rng = OsRng;
    let a = halo2curves::bn256::Fr::random(&mut rng);
    let b = halo2curves::bn256::Fr::random(&mut rng);
    g1.mul_assign(a);
    g2.mul_assign(b);
    c.bench_function("pairing on BN256", move |b| {
        b.iter(|| {
            Bn256::pairing(
                black_box(&halo2curves::bn256::G1Affine::from(g1)),
                black_box(&halo2curves::bn256::G2Affine::from(g2)),
            )
        })
    });
}

criterion_group!(benches, pairing_pluto, pairing_bn256);
criterion_main!(benches);
