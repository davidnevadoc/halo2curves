use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::Field;
use group::prime::PrimeCurveAffine;
//use halo2curves::pluto_eris::G1;
use pasta_curves::arithmetic::CurveExt;
use rand_core::OsRng;

fn criterion_benchmark<G: CurveExt>(c: &mut Criterion, name: &'static str) {
    // G1Projective
    {
        
        let p1 = G::random(OsRng);
        let p2 = G::random(OsRng);
        let p1_affine = G::AffineExt::from(p1);
        let s = G::ScalarExt::random(OsRng);

        const N: usize = 1000;
        let v = vec![G::generator(); N];
        let mut q = vec![G::AffineExt::identity(); N];

        c.bench_function(&format!("G1 Projective of {} check on curve", name), move |b| {
            b.iter(|| black_box(p1).is_on_curve())
        });
        c.bench_function(&format!("G1 Projective of {} check equality", name), move |b| {
            b.iter(|| black_box(p1) == black_box(p1))
        });
        c.bench_function(&format!("G1 Projective of {} to affine", name), move |b| {
            b.iter(|| G::AffineExt::from(black_box(p1)))
        });
        c.bench_function(&format!("G1 Projective of {} doubling", name), move |b| {
            b.iter(|| black_box(p1).double())
        });
        c.bench_function(&format!("G1 Projective of {} addition", name), move |b| {
            b.iter(|| black_box(p1).add(&p2))
        });
        c.bench_function(&format!("G1 Projective of {} mixed addition", name), move |b| {
            b.iter(|| black_box(p2).add(&p1_affine))
        });
        c.bench_function(&format!("G1 Projective of {} scalar multiplication", name), move |b| {
            b.iter(|| black_box(p1) * black_box(s))
        });
        c.bench_function(&format!("G1 Projective of {} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }
}

fn bench_over_pluto(c: &mut Criterion) {
    criterion_benchmark::<halo2curves::pluto_eris::G1>(c, "Pluto")
}

fn bench_over_bn256(c: &mut Criterion) {
    criterion_benchmark::<halo2curves::bn256::G1>(c, "BN256")
}

fn ben_over_bls12_381(c: &mut Criterion) {
    // Pairings
    // {
    //     let g = G1Affine::generator();
    //     let h = G2Affine::generator();
    //     c.bench_function("full pairing", move |b| {
    //         b.iter(|| pairing(black_box(&g), black_box(&h)))
    //     });
    //     c.bench_function("G2 preparation for pairing", move |b| {
    //         b.iter(|| G2Prepared::from(h))
    //     });
    //     let prep = G2Prepared::from(h);
    //     c.bench_function("miller loop for pairing", move |b| {
    //         b.iter(|| multi_miller_loop(&[(&g, &prep)]))
    //     });
    //     let prep = G2Prepared::from(h);
    //     let r = multi_miller_loop(&[(&g, &prep)]);
    //     c.bench_function("final exponentiation for pairing", move |b| {
    //         b.iter(|| r.final_exponentiation())
    //     });
    // }
    // G1Affine
    // {
    //     let name = "G1Affine";
    //     let a = G1Affine::generator();
    //     let s = Scalar::from_raw([1, 2, 3, 4]);
    //     let compressed = [0u8; 48];
    //     let uncompressed = [0u8; 96];
    //     c.bench_function(&format!("{} check on curve", name), move |b| {
    //         b.iter(|| black_box(a).is_on_curve())
    //     });
    //     c.bench_function(&format!("{} check equality", name), move |b| {
    //         b.iter(|| black_box(a) == black_box(a))
    //     });
    //     c.bench_function(&format!("{} scalar multiplication", name), move |b| {
    //         b.iter(|| black_box(a) * black_box(s))
    //     });
    //     c.bench_function(&format!("{} subgroup check", name), move |b| {
    //         b.iter(|| black_box(a).is_torsion_free())
    //     });
    //     c.bench_function(
    //         &format!("{} deserialize compressed point", name),
    //         move |b| b.iter(|| G1Affine::from_compressed(black_box(&compressed))),
    //     );
    //     c.bench_function(
    //         &format!("{} deserialize uncompressed point", name),
    //         move |b| b.iter(|| G1Affine::from_uncompressed(black_box(&uncompressed))),
    //     );
    // }

    // G1Projective
    {
        let name = "BLS12-381";
        let a = halo2curves::bls12_381::G1Projective::generator();
        let a_affine = halo2curves::bls12_381::G1Affine::generator();
        let s = halo2curves::bls12_381::Scalar::from_raw([1, 2, 3, 4]);

        const N: usize = 10000;
        let v = vec![halo2curves::bls12_381::G1Projective::generator(); N];
        let mut q = vec![halo2curves::bls12_381::G1Affine::identity(); N];

        c.bench_function(&format!("G1 Projective of {} check on curve", name), move |b| {
            b.iter(|| black_box(a).is_on_curve())
        });
        c.bench_function(&format!("G1 Projective of {} check equality", name), move |b| {
            b.iter(|| black_box(a) == black_box(a))
        });
        c.bench_function(&format!("G1 Projective of {} to affine", name), move |b| {
            b.iter(|| halo2curves::bls12_381::G1Affine::from(black_box(a)))
        });
        c.bench_function(&format!("G1 Projective of {} doubling", name), move |b| {
            b.iter(|| black_box(a).double())
        });
        c.bench_function(&format!("G1 Projective of {} addition", name), move |b| {
            b.iter(|| black_box(a).add(&a))
        });
        c.bench_function(&format!("G1 Projective of {} mixed addition", name), move |b| {
            b.iter(|| black_box(a).add_mixed(&a_affine))
        });
        c.bench_function(&format!("G1 Projective of {} scalar multiplication", name), move |b| {
            b.iter(|| black_box(a) * black_box(s))
        });
        c.bench_function(&format!("G1 Projective of {} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                halo2curves::bls12_381::G1Projective::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }

    // G2Affine
    // {
    //     let name = "G2Affine";
    //     let a = G2Affine::generator();
    //     let s = Scalar::from_raw([1, 2, 3, 4]);
    //     let compressed = [0u8; 96];
    //     let uncompressed = [0u8; 192];
    //     c.bench_function(&format!("{} check on curve", name), move |b| {
    //         b.iter(|| black_box(a).is_on_curve())
    //     });
    //     c.bench_function(&format!("{} check equality", name), move |b| {
    //         b.iter(|| black_box(a) == black_box(a))
    //     });
    //     c.bench_function(&format!("{} scalar multiplication", name), move |b| {
    //         b.iter(|| black_box(a) * black_box(s))
    //     });
    //     c.bench_function(&format!("{} subgroup check", name), move |b| {
    //         b.iter(|| black_box(a).is_torsion_free())
    //     });
    //     c.bench_function(
    //         &format!("{} deserialize compressed point", name),
    //         move |b| b.iter(|| G2Affine::from_compressed(black_box(&compressed))),
    //     );
    //     c.bench_function(
    //         &format!("{} deserialize uncompressed point", name),
    //         move |b| b.iter(|| G2Affine::from_uncompressed(black_box(&uncompressed))),
    //     );
    // }

    // G2Projective
    // {
    //     let name = "G2Projective";
    //     let a = G2Projective::generator();
    //     let a_affine = G2Affine::generator();
    //     let s = Scalar::from_raw([1, 2, 3, 4]);

    //     const N: usize = 10000;
    //     let v = vec![G2Projective::generator(); N];
    //     let mut q = vec![G2Affine::identity(); N];

    //     c.bench_function(&format!("{} check on curve", name), move |b| {
    //         b.iter(|| black_box(a).is_on_curve())
    //     });
    //     c.bench_function(&format!("{} check equality", name), move |b| {
    //         b.iter(|| black_box(a) == black_box(a))
    //     });
    //     c.bench_function(&format!("{} to affine", name), move |b| {
    //         b.iter(|| G2Affine::from(black_box(a)))
    //     });
    //     c.bench_function(&format!("{} doubling", name), move |b| {
    //         b.iter(|| black_box(a).double())
    //     });
    //     c.bench_function(&format!("{} addition", name), move |b| {
    //         b.iter(|| black_box(a).add(&a))
    //     });
    //     c.bench_function(&format!("{} mixed addition", name), move |b| {
    //         b.iter(|| black_box(a).add_mixed(&a_affine))
    //     });
    //     c.bench_function(&format!("{} scalar multiplication", name), move |b| {
    //         b.iter(|| black_box(a) * black_box(s))
//         });
//         c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
//             b.iter(|| {
//                 G2Projective::batch_normalize(black_box(&v), black_box(&mut q));
//                 black_box(&q)[0]
//             })
//         });
//     }
// }
}


criterion_group!(benches, bench_over_pluto, bench_over_bn256, ben_over_bls12_381);
criterion_main!(benches);
