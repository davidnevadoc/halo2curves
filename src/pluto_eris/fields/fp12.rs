use super::fp::Fp;
use super::fp2::Fp2;
use super::fp6::Fp6;
use crate::ff::Field;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// -GAMMA is a quadratic non-residue in Fp6. Fp12 = Fp6[X]/(X^2 + GAMMA)
/// We introduce the variable w such that w^2 = -GAMMA
// GAMMA = - v
#[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl ConditionallySelectable for Fp12 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp12 {
            c0: Fp6::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp6::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl ConstantTimeEq for Fp12 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl Neg for Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        -&self
    }
}

impl<'a> Neg for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        self.neg()
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn sub(self, rhs: &'b Fp12) -> Fp12 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn add(self, rhs: &'b Fp12) -> Fp12 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, rhs: &'b Fp12) -> Fp12 {
        self.mul(rhs)
    }
}

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    impl_sum_prod,
};
impl_binops_additive!(Fp12, Fp12);
impl_binops_multiplicative!(Fp12, Fp12);
impl_sum_prod!(Fp12);

impl Fp12 {
    #[inline]
    pub const fn zero() -> Self {
        Fp12 {
            c0: Fp6::ZERO,
            c1: Fp6::ZERO,
        }
    }

    #[inline]
    pub const fn one() -> Self {
        Fp12 {
            c0: Fp6::ONE,
            c1: Fp6::ZERO,
        }
    }

    pub fn mul_assign(&mut self, other: &Self) {
        let t0 = self.c0 * other.c0;
        let mut t1 = self.c1 * other.c1;
        let t2 = other.c0 + other.c1;

        self.c1 += &self.c0;
        self.c1 *= &t2;
        self.c1 -= &t0;
        self.c1 -= &t1;

        t1.mul_by_nonresidue();
        self.c0 = t0 + t1;
    }

    pub fn square_assign(&mut self) {
        let mut ab = self.c0 * self.c1;

        let c0c1 = self.c0 + self.c1;

        let mut c0 = self.c1;
        c0.mul_by_nonresidue();
        c0 += &self.c0;
        c0 *= &c0c1;
        c0 -= &ab;
        self.c1 = ab;
        self.c1 += &ab;
        ab.mul_by_nonresidue();
        c0 -= &ab;
        self.c0 = c0;
    }

    pub fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }

    pub fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
    }

    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut t = *other;
        t.mul_assign(self);
        t
    }

    pub fn square(&self) -> Self {
        let mut t = *self;
        t.square_assign();
        t
    }

    #[inline(always)]
    pub fn neg(&self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
        }
    }

    #[inline(always)]
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    pub fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);

        self.c1.c0.mul_assign(&FROBENIUS_COEFF_FP12_C1[power % 12]);
        self.c1.c1.mul_assign(&FROBENIUS_COEFF_FP12_C1[power % 12]);
        self.c1.c2.mul_assign(&FROBENIUS_COEFF_FP12_C1[power % 12]);
    }

    pub fn mul_by_014(&mut self, c0: &Fp2, c1: &Fp2, c4: &Fp2) {
        let mut aa = self.c0;
        aa.mul_by_01(c0, c1);
        let mut bb = self.c1;
        bb.mul_by_1(c4);
        let o = c1 + c4;
        self.c1 += &self.c0;
        self.c1.mul_by_01(c0, &o);
        self.c1 -= &aa;
        self.c1 -= &bb;
        self.c0 = bb;
        self.c0.mul_by_nonresidue();
        self.c0 += &aa;
    }

    pub fn mul_by_034(&mut self, c0: &Fp2, c3: &Fp2, c4: &Fp2) {
        let t0 = Fp6 {
            c0: self.c0.c0 * c0,
            c1: self.c0.c1 * c0,
            c2: self.c0.c2 * c0,
        };
        let mut t1 = self.c1;
        t1.mul_by_01(c3, c4);
        let o = c0 + c3;
        let mut t2 = self.c0 + self.c1;
        t2.mul_by_01(&o, c4);
        t2 -= t0;
        self.c1 = t2 - t1;
        t1.mul_by_nonresidue();
        self.c0 = t0 + t1;
    }

    pub fn invert(&self) -> CtOption<Self> {
        let mut c0s = self.c0;
        c0s.square_assign();
        let mut c1s = self.c1;
        c1s.square_assign();
        c1s.mul_by_nonresidue();
        c0s -= &c1s;

        c0s.invert().map(|t| {
            let mut tmp = Fp12 { c0: t, c1: t };
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1 = tmp.c1.neg();

            tmp
        })
    }

    pub fn cyclotomic_square(&mut self) {
        fn fp4_square(c0: &mut Fp2, c1: &mut Fp2, a0: &Fp2, a1: &Fp2) {
            let t0 = a0.square();
            let t1 = a1.square();
            let mut t2 = t1;
            t2.mul_by_nonresidue();
            *c0 = t2 + t0;
            t2 = a0 + a1;
            t2.square_assign();
            t2 -= t0;
            *c1 = t2 - t1;
        }

        let mut t3 = Fp2::zero();
        let mut t4 = Fp2::zero();
        let mut t5 = Fp2::zero();
        let mut t6 = Fp2::zero();

        fp4_square(&mut t3, &mut t4, &self.c0.c0, &self.c1.c1);
        let mut t2 = t3 - self.c0.c0;
        t2.double_assign();
        self.c0.c0 = t2 + t3;

        t2 = t4 + self.c1.c1;
        t2.double_assign();
        self.c1.c1 = t2 + t4;

        fp4_square(&mut t3, &mut t4, &self.c1.c0, &self.c0.c2);
        fp4_square(&mut t5, &mut t6, &self.c0.c1, &self.c1.c2);

        t2 = t3 - self.c0.c1;
        t2.double_assign();
        self.c0.c1 = t2 + t3;
        t2 = t4 + self.c1.c2;
        t2.double_assign();
        self.c1.c2 = t2 + t4;
        t3 = t6;
        t3.mul_by_nonresidue();
        t2 = t3 + self.c1.c0;
        t2.double_assign();
        self.c1.c0 = t2 + t3;
        t2 = t5 - self.c0.c2;
        t2.double_assign();
        self.c0.c2 = t2 + t5;
    }
}

impl Field for Fp12 {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        Fp12 {
            c0: Fp6::random(&mut rng),
            c1: Fp6::random(&mut rng),
        }
    }

    fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero()
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn double(&self) -> Self {
        self.double()
    }

    fn sqrt(&self) -> CtOption<Self> {
        unimplemented!()
    }

    fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
        unimplemented!()
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }
}

// non_residue^((modulus^i-1)/6) for i=0,...,11
pub const FROBENIUS_COEFF_FP12_C1: [Fp2; 12] = [
    // Fp2(v)**(((p^0) - 1) / 6)
    Fp2::ONE,
    // Fp2(v)**(((p^1) - 1) / 6)
    Fp2 {
        // 0x1b753f8b1c790596306643e8d56a36e8a2ab4f2aaaed2b24e4982068bdf5e6fa81754882f0b9fea9ae5247e985509a8faefbc19a0ebb9c80
        c0: Fp::from_raw([
            0xaefbc19a0ebb9c80,
            0xae5247e985509a8f,
            0x81754882f0b9fea9,
            0xe4982068bdf5e6fa,
            0xa2ab4f2aaaed2b24,
            0x306643e8d56a36e8,
            0x1b753f8b1c790596,
        ]),

        // 0x39a53badefef2817c8f793c5b3d3d674964841b3cf9e8921fde920acb190ffde389f817427cbd34d576f049a35a7126f4bee638825d52d
        c1: Fp::from_raw([
            0x6f4bee638825d52d,
            0x4d576f049a35a712,
            0xde389f817427cbd3,
            0x21fde920acb190ff,
            0x74964841b3cf9e89,
            0x17c8f793c5b3d3d6,
            0x0039a53badefef28,
        ]),
    },
    // Fp2(v)**(((p^2) - 1) / 6)
    Fp2 {
        // 0x24000000000024000130e0000d7f28e4a803ca76be3924a5f43f8cddf9a5c4781b50d5e1ff708dc8d9fa5d8a200bc4398ffff80f80000002
        c0: Fp::from_raw([
            0x8ffff80f80000002,
            0xd9fa5d8a200bc439,
            0x1b50d5e1ff708dc8,
            0xf43f8cddf9a5c478,
            0xa803ca76be3924a5,
            0x0130e0000d7f28e4,
            0x2400000000002400,
        ]),
        c1: Fp::ZERO,
    },
    // Fp2(v)**(((p^3) - 1) / 6)
    Fp2 {
        // 0x199280177f37106928d7340ed9afd69177eabe26c082e02d94386c823bbd0793471c80a54a510f935dd43eebcbb8f8958f5a203e3a9b283c
        c0: Fp::from_raw([
            0x8f5a203e3a9b283c,
            0x5dd43eebcbb8f895,
            0x471c80a54a510f93,
            0x94386c823bbd0793,
            0x77eabe26c082e02d,
            0x28d7340ed9afd691,
            0x199280177f371069,
        ]),

        // 0x1d57dc6eb9af4c9bc8a17db46d26ee52161a3b40a3e6218677467860304c33e21ecdc1bf397519bf52c465ca50b05be07889e7c734730407
        c1: Fp::from_raw([
            0x7889e7c734730407,
            0x52c465ca50b05be0,
            0x1ecdc1bf397519bf,
            0x77467860304c33e2,
            0x161a3b40a3e62186,
            0xc8a17db46d26ee52,
            0x1d57dc6eb9af4c9b,
        ]),
    },
    // Fp2(v)**(((p^4) - 1) / 6)
    Fp2 {
        // 0x480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe
        c0: Fp::from_raw([
            0x100004c37ffffffe,
            0xc8ad8b38dffaf50c,
            0xc956d01c903d720d,
            0x50000d7ee0e4a803,
            0x00000000360001c9,
            0x0000000000004800,
            0x0000000000000000,
        ]),
        c1: Fp::ZERO,
    },
    // Fp2(v)**(((p^5) - 1) / 6)
    Fp2 {
        // 0x12f8405d64503200a92448086be4d44f3571879c7d02418c0faea7cebb61ea6a00bd82d4e450f17039294ab0af03df6601aa17cdb6a93b46
        c0: Fp::from_raw([
            0x01aa17cdb6a93b46,
            0x39294ab0af03df66,
            0x00bd82d4e450f170,
            0x0faea7cebb61ea6a,
            0x3571879c7d02418c,
            0xa92448086be4d44f,
            0x12f8405d64503200,
        ]),

        // 0x66e7e559860e83c20c66ab7daa4aebc1d5346f49c83665faafb38dbfd8ca799e7a144bde2111a44028c13f41520b652b82a26a8436726cd
        c1: Fp::from_raw([
            0xb82a26a8436726cd,
            0x028c13f41520b652,
            0xe7a144bde2111a44,
            0xaafb38dbfd8ca799,
            0x1d5346f49c83665f,
            0x20c66ab7daa4aebc,
            0x066e7e559860e83c,
        ]),
    },
    // Fp2(v)**(((p^6) - 1) / 6)
    Fp2::ONE,
    // Fp2(v)**(((p^7) - 1) / 6)
    Fp2 {
        // 0x1b753f8b1c790596306643e8d56a36e8a2ab4f2aaaed2b24e4982068bdf5e6fa81754882f0b9fea9ae5247e985509a8faefbc19a0ebb9c80
        c0: Fp::from_raw([
            0xaefbc19a0ebb9c80,
            0xae5247e985509a8f,
            0x81754882f0b9fea9,
            0xe4982068bdf5e6fa,
            0xa2ab4f2aaaed2b24,
            0x306643e8d56a36e8,
            0x1b753f8b1c790596,
        ]),

        // 0x39a53badefef2817c8f793c5b3d3d674964841b3cf9e8921fde920acb190ffde389f817427cbd34d576f049a35a7126f4bee638825d52d
        c1: Fp::from_raw([
            0x6f4bee638825d52d,
            0x4d576f049a35a712,
            0xde389f817427cbd3,
            0x21fde920acb190ff,
            0x74964841b3cf9e89,
            0x17c8f793c5b3d3d6,
            0x0039a53badefef28,
        ]),
    },
    // Fp2(v)**(((p^8) - 1) / 6)
    Fp2 {
        // 0x24000000000024000130e0000d7f28e4a803ca76be3924a5f43f8cddf9a5c4781b50d5e1ff708dc8d9fa5d8a200bc4398ffff80f80000002
        c0: Fp::from_raw([
            0x8ffff80f80000002,
            0xd9fa5d8a200bc439,
            0x1b50d5e1ff708dc8,
            0xf43f8cddf9a5c478,
            0xa803ca76be3924a5,
            0x0130e0000d7f28e4,
            0x2400000000002400,
        ]),
        c1: Fp::ZERO,
    },
    // Fp2(v)**(((p^9) - 1) / 6)
    Fp2 {
        //0x199280177f37106928d7340ed9afd69177eabe26c082e02d94386c823bbd0793471c80a54a510f935dd43eebcbb8f8958f5a203e3a9b283c
        c0: Fp::from_raw([
            0x8f5a203e3a9b283c,
            0x5dd43eebcbb8f895,
            0x471c80a54a510f93,
            0x94386c823bbd0793,
            0x77eabe26c082e02d,
            0x28d7340ed9afd691,
            0x199280177f371069,
        ]),

        // 0x1d57dc6eb9af4c9bc8a17db46d26ee52161a3b40a3e6218677467860304c33e21ecdc1bf397519bf52c465ca50b05be07889e7c734730407
        c1: Fp::from_raw([
            0x7889e7c734730407,
            0x52c465ca50b05be0,
            0x1ecdc1bf397519bf,
            0x77467860304c33e2,
            0x161a3b40a3e62186,
            0xc8a17db46d26ee52,
            0x1d57dc6eb9af4c9b,
        ]),
    },
    // Fp2(v)**(((p^10) - 1) / 6)
    Fp2 {
        // 0x480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe
        c0: Fp::from_raw([
            0x100004c37ffffffe,
            0xc8ad8b38dffaf50c,
            0xc956d01c903d720d,
            0x50000d7ee0e4a803,
            0x00000000360001c9,
            0x0000000000004800,
            0x0000000000000000,
        ]),
        c1: Fp::ZERO,
    },
    // Fp2(v)**(((p^11) - 1) / 6)
    Fp2 {
        // 0x12f8405d64503200a92448086be4d44f3571879c7d02418c0faea7cebb61ea6a00bd82d4e450f17039294ab0af03df6601aa17cdb6a93b46
        c0: Fp::from_raw([
            0x01aa17cdb6a93b46,
            0x39294ab0af03df66,
            0x00bd82d4e450f170,
            0x0faea7cebb61ea6a,
            0x3571879c7d02418c,
            0xa92448086be4d44f,
            0x12f8405d64503200,
        ]),

        // 0x66e7e559860e83c20c66ab7daa4aebc1d5346f49c83665faafb38dbfd8ca799e7a144bde2111a44028c13f41520b652b82a26a8436726cd
        c1: Fp::from_raw([
            0xb82a26a8436726cd,
            0x028c13f41520b652,
            0xe7a144bde2111a44,
            0xaafb38dbfd8ca799,
            0x1d5346f49c83665f,
            0x20c66ab7daa4aebc,
            0x066e7e559860e83c,
        ]),
    },
];

#[cfg(test)]
use rand::SeedableRng;
#[cfg(test)]
use rand_xorshift::XorShiftRng;

#[test]
fn test_fp12_mul_by_014() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..1000 {
        let c0 = Fp2::random(&mut rng);
        let c1 = Fp2::random(&mut rng);
        let c5 = Fp2::random(&mut rng);
        let mut a = Fp12::random(&mut rng);
        let mut b = a;

        a.mul_by_014(&c0, &c1, &c5);
        b.mul_assign(&Fp12 {
            c0: Fp6 {
                c0,
                c1,
                c2: Fp2::zero(),
            },
            c1: Fp6 {
                c0: Fp2::zero(),
                c1: c5,
                c2: Fp2::zero(),
            },
        });

        assert_eq!(a, b);
    }
}

#[test]
fn test_fp12_mul_by_034() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..1000 {
        let c0 = Fp2::random(&mut rng);
        let c3 = Fp2::random(&mut rng);
        let c4 = Fp2::random(&mut rng);
        let mut a = Fp12::random(&mut rng);
        let mut b = a;

        a.mul_by_034(&c0, &c3, &c4);
        b.mul_assign(&Fp12 {
            c0: Fp6 {
                c0,
                c1: Fp2::zero(),
                c2: Fp2::zero(),
            },
            c1: Fp6 {
                c0: c3,
                c1: c4,
                c2: Fp2::zero(),
            },
        });

        assert_eq!(a, b);
    }
}

#[test]
fn test_squaring() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..1000 {
        let mut a = Fp12::random(&mut rng);
        let mut b = a;
        b.mul_assign(&a);
        a.square_assign();
        assert_eq!(a, b);
    }
}

#[test]
fn test_frobenius() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..50 {
        for i in 0..13 {
            let mut a = Fp12::random(&mut rng);
            let mut b = a;

            for _ in 0..i {
                a = a.pow_vartime(&[
                    0x9ffffcd300000001,
                    0xa2a7e8c30006b945,
                    0xe4a7a5fe8fadffd6,
                    0x443f9a5cda8a6c7b,
                    0xa803ca76f439266f,
                    0x0130e0000d7f70e4,
                    0x2400000000002400,
                ]);
            }
            b.frobenius_map(i);

            assert_eq!(a, b);
        }
    }
}

#[test]
fn test_field() {
    crate::tests::field::random_field_tests::<Fp12>("fp12".to_string());
}
