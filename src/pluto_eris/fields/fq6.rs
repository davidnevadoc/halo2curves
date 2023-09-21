use super::fq::Fq;
use super::fq2::Fq2;
use crate::ff::Field;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// -BETA is a cubic non-residue in Fq2. Fq6 = Fq2[X]/(X^3 + BETA)
/// We introduce the variable v such that v^3 = -ALPHA
// BETA = - (u+2)

// V_CUBE = u + 2
// const V_CUBE: Fq2 = Fq2 {
//     c0: Fq::from_raw([0x02, 0, 0, 0, 0, 0, 0]),
//     c1: Fq::from_raw([0x01, 0, 0, 0, 0, 0, 0]),
// };

#[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
pub struct Fq6 {
    pub c0: Fq2,
    pub c1: Fq2,
    pub c2: Fq2,
}

impl ConditionallySelectable for Fq6 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fq6 {
            c0: Fq2::conditional_select(&a.c0, &b.c0, choice),
            c1: Fq2::conditional_select(&a.c1, &b.c1, choice),
            c2: Fq2::conditional_select(&a.c2, &b.c2, choice),
        }
    }
}

impl ConstantTimeEq for Fq6 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1) & self.c2.ct_eq(&other.c2)
    }
}

impl Neg for Fq6 {
    type Output = Fq6;

    #[inline]
    fn neg(self) -> Fq6 {
        -&self
    }
}

impl<'a> Neg for &'a Fq6 {
    type Output = Fq6;

    #[inline]
    fn neg(self) -> Fq6 {
        self.neg()
    }
}

impl<'a, 'b> Sub<&'b Fq6> for &'a Fq6 {
    type Output = Fq6;

    #[inline]
    fn sub(self, rhs: &'b Fq6) -> Fq6 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fq6> for &'a Fq6 {
    type Output = Fq6;

    #[inline]
    fn add(self, rhs: &'b Fq6) -> Fq6 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fq6> for &'a Fq6 {
    type Output = Fq6;

    #[inline]
    fn mul(self, rhs: &'b Fq6) -> Fq6 {
        self.mul(rhs)
    }
}

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    impl_sum_prod,
};
impl_binops_additive!(Fq6, Fq6);
impl_binops_multiplicative!(Fq6, Fq6);
impl_sum_prod!(Fq6);

impl Fq6 {
    #[inline]
    pub const fn zero() -> Self {
        Fq6 {
            c0: Fq2::ZERO,
            c1: Fq2::ZERO,
            c2: Fq2::ZERO,
        }
    }

    #[inline]
    pub const fn one() -> Self {
        Fq6 {
            c0: Fq2::ONE,
            c1: Fq2::ZERO,
            c2: Fq2::ZERO,
        }
    }

    pub fn mul_assign(&mut self, other: &Self) {
        let mut a_a = self.c0;
        let mut b_b = self.c1;
        let mut c_c = self.c2;
        a_a *= &other.c0;
        b_b *= &other.c1;
        c_c *= &other.c2;

        let mut t1 = other.c1;
        t1 += &other.c2;
        {
            let mut tmp = self.c1;
            tmp += &self.c2;

            t1 *= &tmp;
            t1 -= &b_b;
            t1 -= &c_c;
            t1.mul_by_nonresidue();
            t1 += &a_a;
        }

        let mut t3 = other.c0;
        t3 += &other.c2;
        {
            let mut tmp = self.c0;
            tmp += &self.c2;

            t3 *= &tmp;
            t3 -= &a_a;
            t3 += &b_b;
            t3 -= &c_c;
        }

        let mut t2 = other.c0;
        t2 += &other.c1;
        {
            let mut tmp = self.c0;
            tmp += &self.c1;

            t2 *= &tmp;
            t2 -= &a_a;
            t2 -= &b_b;
            c_c.mul_by_nonresidue();
            t2 += &c_c;
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = t3;
    }

    pub fn square_assign(&mut self) {
        // s0 = a^2
        let mut s0 = self.c0;
        s0.square_assign();
        // s1 = 2ab
        let mut ab = self.c0;
        ab *= &self.c1;
        let mut s1 = ab;
        s1.double_assign();
        // s2 = (a - b + c)^2
        let mut s2 = self.c0;
        s2 -= &self.c1;
        s2 += &self.c2;
        s2.square_assign();
        // bc
        let mut bc = self.c1;
        bc *= &self.c2;
        // s3 = 2bc
        let mut s3 = bc;
        s3.double_assign();
        // s4 = c^2
        let mut s4 = self.c2;
        s4.square_assign();

        // new c0 = 2bc.mul_by_xi + a^2
        self.c0 = s3;
        self.c0.mul_by_nonresidue();
        // self.c0.mul_by_xi();
        self.c0 += &s0;

        // new c1 = (c^2).mul_by_xi + 2ab
        self.c1 = s4;
        self.c1.mul_by_nonresidue();
        // self.c1.mul_by_xi();
        self.c1 += &s1;

        // new c2 = 2ab + (a - b + c)^2 + 2bc - a^2 - c^2 = b^2 + 2ac
        self.c2 = s1;
        self.c2 += &s2;
        self.c2 += &s3;
        self.c2 -= &s0;
        self.c2 -= &s4;
    }

    pub fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
            c2: self.c2.double(),
        }
    }

    pub fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
        self.c2 = self.c2.double();
    }

    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
            c2: self.c2 + other.c2,
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
            c2: self.c2 - other.c2,
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

    pub fn neg(&self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }

    pub fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c2.frobenius_map(power);

        self.c1.mul_assign(&FROBENIUS_COEFF_FQ6_C1[power % 6]);
        self.c2.mul_assign(&FROBENIUS_COEFF_FQ6_C2[power % 6]);
    }

    /// Multiply by cubic nonresidue v.
    pub fn mul_by_nonresidue(&mut self) {
        use std::mem::swap;
        swap(&mut self.c0, &mut self.c1);
        swap(&mut self.c0, &mut self.c2);
        // c0, c1, c2 -> c2, c0, c1
        self.c0.mul_by_nonresidue();
    }

    // /// Multiply by cubic nonresidue v.
    // pub fn mul_by_v(&mut self) {
    //     use std::mem::swap;
    //     swap(&mut self.c0, &mut self.c1);
    //     swap(&mut self.c0, &mut self.c2);

    //     self.c0.mul_by_xi();
    // }

    pub fn mul_by_1(&mut self, c1: &Fq2) {
        let mut b_b = self.c1;
        b_b *= c1;

        let mut t1 = *c1;
        {
            let mut tmp = self.c1;
            tmp += &self.c2;

            t1 *= &tmp;
            t1 -= &b_b;
            t1.mul_by_nonresidue();
        }

        let mut t2 = *c1;
        {
            let mut tmp = self.c0;
            tmp += &self.c1;

            t2 *= &tmp;
            t2 -= &b_b;
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = b_b;
    }

    pub fn mul_by_01(&mut self, c0: &Fq2, c1: &Fq2) {
        let mut a_a = self.c0;
        let mut b_b = self.c1;
        a_a *= c0;
        b_b *= c1;

        let mut t1 = *c1;
        {
            let mut tmp = self.c1;
            tmp += &self.c2;

            t1 *= &tmp;
            t1 -= &b_b;
            t1.mul_by_nonresidue();
            t1 += &a_a;
        }

        let mut t3 = *c0;
        {
            let mut tmp = self.c0;
            tmp += &self.c2;

            t3 *= &tmp;
            t3 -= &a_a;
            t3 += &b_b;
        }

        let mut t2 = *c0;
        t2 += c1;
        {
            let mut tmp = self.c0;
            tmp += &self.c1;

            t2 *= &tmp;
            t2 -= &a_a;
            t2 -= &b_b;
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = t3;
    }

    fn invert(&self) -> CtOption<Self> {
        let mut c0 = self.c2;
        c0.mul_by_nonresidue();
        c0 *= &self.c1;
        c0 = -c0;
        {
            let mut c0s = self.c0;
            c0s.square_assign();
            c0 += &c0s;
        }
        let mut c1 = self.c2;
        c1.square_assign();
        c1.mul_by_nonresidue();
        {
            let mut c01 = self.c0;
            c01 *= &self.c1;
            c1 -= &c01;
        }
        let mut c2 = self.c1;
        c2.square_assign();
        {
            let mut c02 = self.c0;
            c02 *= &self.c2;
            c2 -= &c02;
        }

        let mut tmp1 = self.c2;
        tmp1 *= &c1;
        let mut tmp2 = self.c1;
        tmp2 *= &c2;
        tmp1 += &tmp2;
        tmp1.mul_by_nonresidue();
        tmp2 = self.c0;
        tmp2 *= &c0;
        tmp1 += &tmp2;

        tmp1.invert().map(|t| {
            let mut tmp = Fq6 {
                c0: t,
                c1: t,
                c2: t,
            };
            tmp.c0 *= &c0;
            tmp.c1 *= &c1;
            tmp.c2 *= &c2;

            tmp
        })
    }
}

impl Field for Fq6 {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        Fq6 {
            c0: Fq2::random(&mut rng),
            c1: Fq2::random(&mut rng),
            c2: Fq2::random(&mut rng),
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

pub const FROBENIUS_COEFF_FQ6_C1: [Fq2; 6] = [
    // Fq2(v^3)**(((q^0) - 1) / 3)
    Fq2::ONE,
    // Fq2(v^3)**(((q^1) - 1) / 3)
    Fq2 {
        // 0xa12a7a5bb16c30e70ed772cd3a465fb9a95b4ac5841e111d3ee22591d2c9f387870fd58bb208ac369b0926203e0d8b2672a2f804be5c3d2
        c0: Fq::from_raw([
            0x672a2f804be5c3d2,
            0x69b0926203e0d8b2,
            0x7870fd58bb208ac3,
            0xd3ee22591d2c9f38,
            0x9a95b4ac5841e111,
            0x70ed772cd3a465fb,
            0x0a12a7a5bb16c30e,
        ]),

        // 0x169c69ad87060cd2f94b547d64e48b8eb2b3f55438c0bc3e38a1914bdb01e208918d3b6fbd6061efa04dc91e9dc401c5b0aa20a5bf27d84b
        c1: Fq::from_raw([
            0xb0aa20a5bf27d84b,
            0xa04dc91e9dc401c5,
            0x918d3b6fbd6061ef,
            0x38a1914bdb01e208,
            0xb2b3f55438c0bc3e,
            0xf94b547d64e48b8e,
            0x169c69ad87060cd2,
        ]),
    },
    // Fq2(v^3)**(((q^2) - 1) / 3)
    Fq2 {
        // 0x480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe
        c0: Fq::from_raw([
            0x100004c37ffffffe,
            0xc8ad8b38dffaf50c,
            0xc956d01c903d720d,
            0x50000d7ee0e4a803,
            0x00000000360001c9,
            0x0000000000004800,
            0x0000000000000000,
        ]),
        c1: Fq::ZERO,
    },
    // Fq2(v^3)**(((q^3) - 1) / 3)
    Fq2 {
        // 0x1fdb6f538b54ca12ddd30422cf76537d2a39e3fc90d2f7c1b94fd59ed356516c03ee97d0838d20874b647fb3feaa9e1269546ccd30584139
        c0: Fq::from_raw([
            0x69546ccd30584139,
            0x4b647fb3feaa9e12,
            0x03ee97d0838d2087,
            0xb94fd59ed356516c,
            0x2a39e3fc90d2f7c1,
            0xddd30422cf76537d,
            0x1fdb6f538b54ca12,
        ]),

        // 0x1c1612fc50dc8de03e08ba7d3bff3bbc8c0a66e2ad6fb668d641c9c0dec7251d2c3d56b69469165b567458d579b33ac78024a5443680656e
        c1: Fq::from_raw([
            0x8024a5443680656e,
            0x567458d579b33ac7,
            0x2c3d56b69469165b,
            0xd641c9c0dec7251d,
            0x8c0a66e2ad6fb668,
            0x3e08ba7d3bff3bbc,
            0x1c1612fc50dc8de0,
        ]),
    },
    // Fq2(v^3)**(((q^4) - 1) / 3)
    Fq2 {
        // 0x24000000000024000130e0000d7f28e4a803ca76be3924a5f43f8cddf9a5c4781b50d5e1ff708dc8d9fa5d8a200bc4398ffff80f80000002
        c0: Fq::from_raw([
            0x8ffff80f80000002,
            0xd9fa5d8a200bc439,
            0x1b50d5e1ff708dc8,
            0xf43f8cddf9a5c478,
            0xa803ca76be3924a5,
            0x0130e0000d7f28e4,
            0x2400000000002400,
        ]),

        c1: Fq::ZERO,
    },
    // Fq2(v^3)**(((q^5) - 1) / 3)
    Fq2 {
        // 0x1e11e906b994badeb3a144b077e428508b37fc44ff5d740afb413cc1c491e8534cefb6d3e0ae5462903abf6ffd81fbc66f815d5883c1faf7
        c0: Fq::from_raw([
            0x6f815d5883c1faf7,
            0x903abf6ffd81fbc6,
            0x4cefb6d3e0ae5462,
            0xfb413cc1c491e853,
            0x8b37fc44ff5d740a,
            0xb3a144b077e42850,
            0x1e11e906b994bade,
        ]),

        // 0x154d8356281dad4ccb0db1057a1b1a7e114938b70241da37799bd9acfb4bd1d20b84b9d6cd9287624e8daf91e89635fe0f3133bc0a57c249
        c1: Fq::from_raw([
            0x0f3133bc0a57c249,
            0x4e8daf91e89635fe,
            0x0b84b9d6cd928762,
            0x799bd9acfb4bd1d2,
            0x114938b70241da37,
            0xcb0db1057a1b1a7e,
            0x154d8356281dad4c,
        ]),
    },
];

pub const FROBENIUS_COEFF_FQ6_C2: [Fq2; 6] = [
    // Fq2(v^3)**(((2q^0) - 2) / 3)
    Fq2::ONE,
    // Fq2(v^3)**(((2q^1) - 2) / 3)
    Fq2 {
        // 0x1d51c23fcb1f9dae6458b6193a3766305bedb614a7feb8ec89f2de54da4204fe6832f879b2d457a239174d60a6247988fc8442b768993461
        c0: Fq::from_raw([
            0xfc8442b768993461,
            0x39174d60a6247988,
            0x6832f879b2d457a2,
            0x89f2de54da4204fe,
            0x5bedb614a7feb8ec,
            0x6458b6193a376630,
            0x1d51c23fcb1f9dae,
        ]),

        // 0x1500a70695886e6cab571195c7d6000acb2b830a12825c2b4c494d816a5f107ce3393a7f8808224892191b8bc5a401e57acbf5e1ff7aa60d
        c1: Fq::from_raw([
            0x7acbf5e1ff7aa60d,
            0x92191b8bc5a401e5,
            0xe3393a7f88082248,
            0x4c494d816a5f107c,
            0xcb2b830a12825c2b,
            0xab571195c7d6000a,
            0x1500a70695886e6c,
        ]),
    },
    // Fq2(v^3)**(((2q^2) - 2) / 3)
    Fq2 {
        // 0x24000000000024000130e0000d7f28e4a803ca76be3924a5f43f8cddf9a5c4781b50d5e1ff708dc8d9fa5d8a200bc4398ffff80f80000002
        c0: Fq::from_raw([
            0x8ffff80f80000002,
            0xd9fa5d8a200bc439,
            0x1b50d5e1ff708dc8,
            0xf43f8cddf9a5c478,
            0xa803ca76be3924a5,
            0x0130e0000d7f28e4,
            0x2400000000002400,
        ]),
        c1: Fq::ZERO,
    },
    // Fq2(v^3)**(((2q^3) - 2) / 3)
    Fq2 {
        // 0x120be5e275b7740b992ac11ed87340811340a18ffdae8c5a75a401e1859b4c46081d3ffa689cca53b3084b8347191c1df788b62181838924
        c0: Fq::from_raw([
            0xf788b62181838924,
            0xb3084b8347191c1d,
            0x081d3ffa689cca53,
            0x75a401e1859b4c46,
            0x1340a18ffdae8c5a,
            0x992ac11ed8734081,
            0x120be5e275b7740b,
        ]),
        //0xf64268ae9d89108841b2e2044b239b513e060979a5412ff8a3b20a3c8faec3a4b83661cd06f3f91e511eb36374fe08bf6cc480f091b565
        c1: Fq::from_raw([
            0xbf6cc480f091b565,
            0x1e511eb36374fe08,
            0xa4b83661cd06f3f9,
            0xf8a3b20a3c8faec3,
            0x513e060979a5412f,
            0x8841b2e2044b239b,
            0x00f64268ae9d8910,
        ]),
    },
    // Fq2(v^3)**(((2q^4) - 2) / 3)
    Fq2 {
        // 0x480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe
        c0: Fq::from_raw([
            0x100004c37ffffffe,
            0xc8ad8b38dffaf50c,
            0xc956d01c903d720d,
            0x50000d7ee0e4a803,
            0x00000000360001c9,
            0x0000000000004800,
            0x0000000000000000,
        ]),
        c1: Fq::ZERO,
    },
    // Fq2(v^3)**(((2q^5) - 2) / 3)
    Fq2 {
        //0x18a257ddbf29364604de48c808543b17e0d93d4942c5079788e85483553787b358ff138903eaddb7593038a212cfdce44bf300cd15e3427d
        c0: Fq::from_raw([
            0x4bf300cd15e3427d,
            0x593038a212cfdce4,
            0x58ff138903eaddb7,
            0x88e85483553787b3,
            0xe0d93d4942c50797,
            0x04de48c808543b17,
            0x18a257ddbf293646,
        ]),

        c1: Fq::from_raw(
            // 0xe091690bbda2c82cd981b88415e4d3e8b9a416368118913ff529ad1339bad3b5cb6351d3a9ee994f23dae83d6edb95765c742700ff3a48f
            [
                0x65c742700ff3a48f,
                0xf23dae83d6edb957,
                0x5cb6351d3a9ee994,
                0xff529ad1339bad3b,
                0x8b9a416368118913,
                0xcd981b88415e4d3e,
                0x0e091690bbda2c82,
            ],
        ),
    },
];

#[cfg(test)]
use rand::SeedableRng;
#[cfg(test)]
use rand_xorshift::XorShiftRng;

#[test]
fn test_fq6_mul_nonresidue() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let nqr = Fq6 {
        c0: Fq2::zero(),
        c1: Fq2::one(),
        c2: Fq2::zero(),
    };

    for _ in 0..1000 {
        let mut a = Fq6::random(&mut rng);
        let mut b = a;
        a.mul_by_nonresidue();
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq6_mul_by_1() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..1000 {
        let c1 = Fq2::random(&mut rng);
        let mut a = Fq6::random(&mut rng);
        let mut b = a;

        a.mul_by_1(&c1);
        b.mul_assign(&Fq6 {
            c0: Fq2::zero(),
            c1,
            c2: Fq2::zero(),
        });

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq6_mul_by_01() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..1000 {
        let c0 = Fq2::random(&mut rng);
        let c1 = Fq2::random(&mut rng);
        let mut a = Fq6::random(&mut rng);
        let mut b = a;

        a.mul_by_01(&c0, &c1);
        b.mul_assign(&Fq6 {
            c0,
            c1,
            c2: Fq2::zero(),
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
        let mut a = Fq6::random(&mut rng);
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
        for i in 0..8 {
            let mut a = Fq6::random(&mut rng);
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
    crate::tests::field::random_field_tests::<Fq6>("fq6".to_string());
}
