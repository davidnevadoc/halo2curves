use super::fq12::Fq12;
use super::fq2::Fq2;
use super::{Fr, G1Affine, G2Affine, BLS_X, G1, G2};
use crate::ff_ext::quadratic::QuadSparseMul;
use crate::ff_ext::ExtField;
use core::borrow::Borrow;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use ff::PrimeField;
use group::prime::PrimeCurveAffine;
use group::Group;
use pairing::{Engine, MillerLoopResult, MultiMillerLoop, PairingCurveAffine};
use rand::RngCore;
use std::ops::MulAssign;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

crate::impl_gt!(Gt, Fq12, Fr);
crate::impl_miller_loop_components!(Bls12381, G1, G1Affine, G2, G2Affine, Fq12, Gt, Fr);

impl MillerLoopResult for Fq12 {
    type Gt = Gt;

    fn final_exponentiation(&self) -> Gt {
        #[must_use]
        fn exp_by_x(f: Fq12) -> Fq12 {
            let mut acc = Fq12::one();
            for (i, b) in BLS_X.into_iter().enumerate() {
                (i != 0).then(|| acc.cyclotomic_square());
                (b == 1).then(|| acc *= f);
            }
            acc.conjugate();
            acc
        }

        let mut t0 = *self;
        t0.frobenius_map(6);

        Gt(self
            .invert()
            .map(|mut t1| {
                let mut t2 = t0 * t1;
                t1 = t2;
                t2.frobenius_map(2);
                t2 *= t1;
                t1 = t2;
                t1.cyclotomic_square();
                t1.conjugate();
                let mut t3 = exp_by_x(t2);
                let mut t4 = t3;
                t4.cyclotomic_square();
                let mut t5 = t1 * t3;
                t1 = exp_by_x(t5);
                t0 = exp_by_x(t1);
                let mut t6 = exp_by_x(t0) * t4;
                t4 = exp_by_x(t6);
                t5.conjugate();
                t4 *= t5 * t2;
                t1 *= t2;
                t1.frobenius_map(3);
                t2.conjugate();
                t6 *= t2;
                t6.frobenius_map(1);
                t3 *= t0;
                t3.frobenius_map(2);
                t3 * t4 * t1 * t6
            })
            .unwrap())
    }
}

pub fn multi_miller_loop(terms: &[(&G1Affine, &G2Affine)]) -> Fq12 {
    let terms = terms
        .iter()
        .filter_map(|&(p, q)| {
            if bool::from(p.is_identity()) || bool::from(q.is_identity()) {
                None
            } else {
                Some((p, q))
            }
        })
        .collect::<Vec<_>>();

    let mut f = Fq12::one();
    let mut r = terms.iter().map(|(_, q)| q.to_curve()).collect::<Vec<_>>();

    for (i, x) in BLS_X.iter().map(|&b| b == 1).skip(1).enumerate() {
        if i != 0 {
            f.square_assign();
        }

        terms.iter().zip(r.iter_mut()).for_each(|((p, _), r)| {
            double(&mut f, r, p);
        });

        if x {
            for ((p, q), r) in terms.iter().zip(r.iter_mut()) {
                add(&mut f, r, q, p);
            }
        }
    }

    f.conjugate();
    f
}

fn ell(f: &mut Fq12, coeffs: &(Fq2, Fq2, Fq2), p: &G1Affine) {
    let mut c0 = coeffs.0;
    let mut c1 = coeffs.1;
    c0.c0.mul_assign(&p.y);
    c0.c1.mul_assign(&p.y);
    c1.c0.mul_assign(&p.x);
    c1.c1.mul_assign(&p.x);
    Fq12::mul_by_014(f, &coeffs.2, &c1, &c0);
}

#[cfg(test)]
mod test {
    use super::super::{Bls12381, Fr, G1, G2};
    use super::{multi_miller_loop, Fq12, G1Affine, G2Affine, Gt};
    use ff::Field;
    use group::{prime::PrimeCurveAffine, Curve, Group};
    use pairing::{Engine as _, MillerLoopResult, PairingCurveAffine};
    use rand_core::OsRng;
    crate::test_pairing!(Bls12381, G1, G1Affine, G2, G2Affine, Fq12, Gt, Fr);
}
