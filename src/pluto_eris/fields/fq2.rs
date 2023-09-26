use super::fq::{Fq, MODULUS_STR};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::legendre::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use std::ops::MulAssign;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

/// -ALPHA is a quadratic non-residue in Fq. Fq2 = Fq[X]/(X^2 + ALPHA)
/// We introduce the variable u such that u^2 = -ALPHA
// const ALPHA: u64 = 5;

// U_SQUARE = -5
// 0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5cda8a6c7be4a7a5fe8fadffd6a2a7e8c30006b9459ffffcd2fffffffc
const U_SQUARE: Fq = Fq::from_raw([
    0x9ffffcd2fffffffc,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

const NEG_ONE: Fq2 = Fq2 {
    c0: super::fq::NEG_ONE,
    c1: Fq::ZERO,
};

/// An element of Fq2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    #[inline(always)]
    fn cmp(&self, other: &Fq2) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fq2) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ConditionallySelectable for Fq2 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fq2 {
            c0: Fq::conditional_select(&a.c0, &b.c0, choice),
            c1: Fq::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl ConstantTimeEq for Fq2 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl Default for Fq2 {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl From<Fq2> for [u8; 112] {
    fn from(value: Fq2) -> [u8; 112] {
        value.to_bytes()
    }
}

impl<'a> From<&'a Fq2> for [u8; 112] {
    fn from(value: &'a Fq2) -> [u8; 112] {
        value.to_bytes()
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        -&self
    }
}

impl<'a> Neg for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        self.neg()
    }
}

impl<'a, 'b> Sub<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn sub(self, rhs: &'b Fq2) -> Fq2 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn add(self, rhs: &'b Fq2) -> Fq2 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn mul(self, rhs: &'b Fq2) -> Fq2 {
        self.mul(rhs)
    }
}

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    impl_sum_prod,
};
impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
impl_sum_prod!(Fq2);

const SIZE: usize = 112;
const COORD_SIZE: usize = 56;

impl Fq2 {
    #[inline]
    pub const fn zero() -> Fq2 {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    #[inline]
    pub const fn one() -> Fq2 {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    pub const fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0, c1 }
    }

    pub const fn size() -> usize {
        SIZE
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fq`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; SIZE]) -> CtOption<Fq2> {
        let c0 = Fq::from_bytes(bytes[0..COORD_SIZE].try_into().unwrap());
        let c1 = Fq::from_bytes(bytes[COORD_SIZE..SIZE].try_into().unwrap());
        CtOption::new(
            Fq2 {
                c0: c0.unwrap(),
                c1: c1.unwrap(),
            },
            c0.is_some() & c1.is_some(),
        )
    }

    /// Converts an element of `Fq` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; SIZE] {
        let mut res = [0u8; SIZE];
        let c0_bytes = self.c0.to_bytes();
        let c1_bytes = self.c1.to_bytes();
        res[0..COORD_SIZE].copy_from_slice(&c0_bytes[..]);
        res[COORD_SIZE..SIZE].copy_from_slice(&c1_bytes[..]);
        res
    }

    // TODO: This is a naive method using 4 multiplications
    pub fn mul_assign(&mut self, other: &Self) {
        // r0 = s0 * s0 + U_SQUARE * s1 * o1
        // r1 = s0 * o1 - s1 * o0

        let t0 = self.c0 * other.c0;
        let t1 = self.c0 * other.c1;
        let t2 = self.c1 * other.c0;
        let t3 = self.c1 * other.c1;

        self.c0 = t0 + U_SQUARE * t3;
        self.c1 = t1 + t2
    }

    // TODO: This is a naive method using 3 multiplications
    pub fn square_assign(&mut self) {
        // r0 = s0^2 + U_SQUARE * s1^2
        // r1 = 2* s0s1

        let ab = self.c0 * self.c1;
        let a2 = self.c0 * self.c0;
        let b2 = self.c1 * self.c1;

        self.c1 = ab.double();
        self.c0 = a2 + U_SQUARE * b2;
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
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.sub(&other.c0),
            c1: self.c1.sub(&other.c1),
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
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }

    // conjucate by negating c1
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    pub fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            self.conjugate()
        }
    }

    /// Multiply this element by cubic nonresidue: u+2
    pub fn mul_by_nonresidue(&mut self) {
        // (x + y * u)(u + 2) = (2*x + u_square * y) + u (x + 2y)

        let t0 = self.c0;
        let t1 = self.c1;

        self.c0 = t0.double() + U_SQUARE * t1;
        self.c1 = t1.double() + t0;
    }

    pub fn invert(&self) -> CtOption<Self> {
        let mut t1 = self.c1;
        t1 = t1.square();
        t1 *= U_SQUARE;
        let mut t0 = self.c0;
        t0 = t0.square();
        //t0 = c0^2 - U_SQUARE c1^2
        t0 -= &t1;
        t0.invert().map(|t| {
            let mut tmp = Fq2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0 *= &t;
            tmp.c1 *= &t;
            tmp.c1 = -tmp.c1;

            tmp
        })
    }
}

impl Legendre for Fq2 {
    type BasePrimeField = Fq;
    fn legendre_exp() -> &'static [u64] {
        Self::BasePrimeField::legendre_exp()
    }

    /// Norm of Fq2 as extension field in u over Fq
    fn norm(&self) -> Self::BasePrimeField {
        // norm = self * self.cojungate()
        let t0 = self.c0.square();
        let t1 = self.c1.square() * U_SQUARE;
        t1 - t0
    }
}

impl Field for Fq2 {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        Fq2 {
            c0: Fq::random(&mut rng),
            c1: Fq::random(&mut rng),
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
        // Algorithm 10, https://eprint.iacr.org/2012/685.pdf

        // Constants and aux variables

        // Aux elements. Described in PRECOMPUTATION of Algorithm 10.
        // As element of Fq2: E = 0 +  U *
        // 0x13e275a1fa6a13af7a82a3d83bc9e63a667c70cf991a36e603b21f15823a404a021848271d63f0875d232408689b4c6c67153f9701e19938
        const E: Fq2 = Fq2 {
            c0: Fq::ZERO,
            c1: Fq::from_raw([
                0x67153f9701e19938,
                0x5d232408689b4c6c,
                0x021848271d63f087,
                0x03b21f15823a404a,
                0x667c70cf991a36e6,
                0x7a82a3d83bc9e63a,
                0x13e275a1fa6a13af,
            ]),
        };

        // As element of Fq2: f = 5 + 0 * U
        // 0x5
        const F: Fq2 = Fq2 {
            c0: Fq::from_raw([0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            c1: Fq::ZERO,
        };

        // Algorithm (not constant time)
        let b = self.pow_vartime(&[
            // (q-1)/4 =
            // 0x900000000000900004c3800035fdc392a00f29dbd0e499bd10fe69736a29b1ef929e97fa3eb7ff5a8a9fa30c001ae5167ffff34c0000000
            0x67ffff34c0000000,
            0xa8a9fa30c001ae51,
            0xf929e97fa3eb7ff5,
            0xd10fe69736a29b1e,
            0x2a00f29dbd0e499b,
            0x004c3800035fdc39,
            0x0900000000000900,
        ]);

        let b_2 = b.square();
        let mut b_2_q = b_2.clone();
        b_2_q.frobenius_map(1);

        let a0 = b_2_q * b_2;
        if a0 == NEG_ONE {
            CtOption::new(a0, Choice::from(0))
        } else {
            let mut x = b.clone();
            x.frobenius_map(1);
            if x * b == Fq2::ONE {
                let x0 = (b_2 * self).c0.sqrt().unwrap();
                x.c0.mul_assign(x0);
                x.c1.mul_assign(x0);
                CtOption::new(x, Choice::from(1))
            } else {
                let x0 = (self * b_2 * F).sqrt().unwrap();
                x *= x0 * E;
                CtOption::new(x, Choice::from(1))
            }
        }
    }

    // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

    // if self.is_zero().into() {
    //     CtOption::new(Self::ZERO, Choice::from(1))
    // } else {
    //     // a1 = self^((q - 3) / 4)
    //     // 0xc19139cb84c680a6e14116da060561765e05aa45a1c72a34f082305b61f3f51
    //     let u: [u64; 4] = [
    //         0x4f082305b61f3f51,
    //         0x65e05aa45a1c72a3,
    //         0x6e14116da0605617,
    //         0x0c19139cb84c680a,
    //     ];
    //     let mut a1 = self.pow(&u);
    //     let mut alpha = a1;

    //     alpha.square_assign();
    //     alpha.mul_assign(self);
    //     let mut a0 = alpha;
    //     a0.frobenius_map(1);
    //     a0.mul_assign(&alpha);

    //     let neg1 = Fq2 {
    //         // c0: NEGATIVE_ONE,
    //         c0: -Fq::one(),
    //         c1: Fq::zero(),
    //     };

    //     if a0 == neg1 {
    //         CtOption::new(a0, Choice::from(0))
    //     } else {
    //         a1.mul_assign(self);

    //         if alpha == neg1 {
    //             a1.mul_assign(&Fq2 {
    //                 c0: Fq::zero(),
    //                 c1: Fq::one(),
    //             });
    //         } else {
    //             alpha += &Fq2::ONE;
    //             // alpha = alpha^((q - 1) / 2)
    //             // 0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3
    //             let u: [u64; 4] = [
    //                 0x9e10460b6c3e7ea3,
    //                 0xcbc0b548b438e546,
    //                 0xdc2822db40c0ac2e,
    //                 0x183227397098d014,
    //             ];
    //             alpha = alpha.pow(&u);
    //             a1.mul_assign(&alpha);
    //         }
    //         CtOption::new(a1, Choice::from(1))
    //     }
    // }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }
}

impl From<bool> for Fq2 {
    fn from(bit: bool) -> Fq2 {
        if bit {
            Fq2::ONE
        } else {
            Fq2::ZERO
        }
    }
}

impl From<u64> for Fq2 {
    fn from(val: u64) -> Self {
        Fq2 {
            c0: Fq::from(val),
            c1: Fq::zero(),
        }
    }
}

// This trait is only implemented to satisfy the requirement of CurveExt
impl PrimeField for Fq2 {
    type Repr = Fq2Bytes;

    const MODULUS: &'static str = MODULUS_STR;
    const MULTIPLICATIVE_GENERATOR: Self = Fq2 {
        c0: Fq::MULTIPLICATIVE_GENERATOR,
        c1: Fq::ZERO,
    };
    const NUM_BITS: u32 = 446;
    const CAPACITY: u32 = 445;
    const S: u32 = 0;
    // TODO: Check that we can just 0 this and forget.
    const ROOT_OF_UNITY: Self = Fq2::zero();
    const ROOT_OF_UNITY_INV: Self = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };
    const DELTA: Self = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };
    const TWO_INV: Self = Fq2 {
        c0: Fq::TWO_INV,
        c1: Fq::zero(),
    };

    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        let c0 = Fq::from_bytes(&repr.0[..COORD_SIZE].try_into().unwrap());
        let c1 = Fq::from_bytes(&repr.0[COORD_SIZE..].try_into().unwrap());
        // Disallow overflow representation
        CtOption::new(Fq2::new(c0.unwrap(), c1.unwrap()), Choice::from(1))
    }

    fn to_repr(&self) -> Self::Repr {
        Fq2Bytes(self.to_bytes())
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr().as_ref()[0] & 1)
    }
}

impl FromUniformBytes<64> for Fq2 {
    fn from_uniform_bytes(bytes: &[u8; 64]) -> Self {
        Self::new(Fq::from_uniform_bytes(bytes), Fq::zero())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Fq2Bytes([u8; SIZE]);

impl Default for Fq2Bytes {
    fn default() -> Self {
        Self([0u8; SIZE])
    }
}

impl AsMut<[u8]> for Fq2Bytes {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl AsRef<[u8]> for Fq2Bytes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl crate::serde::SerdeObject for Fq2 {
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        debug_assert_eq!(bytes.len(), 112);
        let [c0, c1] = [0, 56].map(|i| Fq::from_raw_bytes_unchecked(&bytes[i..i + 56]));
        Self { c0, c1 }
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != SIZE {
            return None;
        }
        let [c0, c1] = [0, COORD_SIZE].map(|i| Fq::from_raw_bytes(&bytes[i..i + COORD_SIZE]));
        c0.zip(c1).map(|(c0, c1)| Self { c0, c1 })
    }
    fn to_raw_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(SIZE);
        for limb in self.c0.0.iter().chain(self.c1.0.iter()) {
            res.extend_from_slice(&limb.to_le_bytes());
        }
        res
    }
    fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
        let [c0, c1] = [(); 2].map(|_| Fq::read_raw_unchecked(reader));
        Self { c0, c1 }
    }
    fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let c0 = Fq::read_raw(reader)?;
        let c1 = Fq::read_raw(reader)?;
        Ok(Self { c0, c1 })
    }
    fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.c0.write_raw(writer)?;
        self.c1.write_raw(writer)
    }
}

impl WithSmallOrderMulGroup<3> for Fq2 {
    const ZETA: Self = Fq2 {
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
        c1: Fq::zero(),
    };
}

#[cfg(test)]
use rand::SeedableRng;
#[cfg(test)]
use rand_xorshift::XorShiftRng;

#[test]
fn test_ser() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a0 = Fq2::random(&mut rng);
    let a_bytes = a0.to_bytes();
    let a1 = Fq2::from_bytes(&a_bytes).unwrap();
    assert_eq!(a0, a1);
}

#[test]
fn test_fq2_ordering() {
    let mut a = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };

    let mut b = a;

    assert!(a.cmp(&b) == Ordering::Equal);
    b.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Equal);
    b.c1 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c1 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Greater);
    b.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Equal);
}

#[test]
fn test_fq2_basics() {
    assert_eq!(
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        },
        Fq2::ZERO
    );
    assert_eq!(
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        },
        Fq2::ONE
    );
    assert_eq!(Fq2::ZERO.is_zero().unwrap_u8(), 1);
    assert_eq!(Fq2::ONE.is_zero().unwrap_u8(), 0);
    assert_eq!(
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::one(),
        }
        .is_zero()
        .unwrap_u8(),
        0
    );
}

#[test]
fn test_fq2_squaring() {
    // u + 1
    let mut a = Fq2 {
        c0: Fq::one(),
        c1: Fq::one(),
    };
    // (u + 1) ^2 = 1 + u^2 + 2u = -4 + 2u
    a.square_assign();
    let minus_4 = -Fq::from(4u64);
    assert_eq!(
        a,
        Fq2 {
            c0: minus_4,
            c1: Fq::one() + Fq::one(),
        }
    );

    // u
    let mut a = Fq2 {
        c0: Fq::zero(),
        c1: Fq::one(),
    };
    // u^2
    a.square_assign();
    assert_eq!(
        a,
        Fq2 {
            c0: U_SQUARE,
            c1: Fq::zero(),
        }
    );
}

#[test]
fn test_fq2_mul_nonresidue() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    let nqr = Fq2 {
        c0: Fq::one() + Fq::one(),
        c1: Fq::one(),
    };

    for _ in 0..1000 {
        let mut a = Fq2::random(&mut rng);
        let mut b = a;
        a.mul_by_nonresidue();
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
pub fn test_sqrt() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    const N_ITER: usize = 1000;
    for _ in 0..N_ITER {
        let a = Fq2::random(&mut rng);
        if a.legendre() == -Fq::ONE {
            assert!(bool::from(a.sqrt().is_none()));
        }
    }

    for _ in 0..N_ITER {
        let a = Fq2::random(&mut rng);
        let mut b = a;
        b.square_assign();
        assert_eq!(b.legendre(), Fq::ONE);

        let b = b.sqrt().unwrap();
        let mut negb = b;
        negb = negb.neg();

        assert!(a == b || a == negb);
    }

    let mut c = Fq2::ONE;
    for _ in 0..N_ITER {
        let mut b = c;
        b.square_assign();
        assert_eq!(b.legendre(), Fq::ONE);

        b = b.sqrt().unwrap();

        if b != c {
            b = b.neg();
        }

        assert_eq!(b, c);

        c += &Fq2::ONE;
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
            let mut a = Fq2::random(&mut rng);
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
    crate::tests::field::random_field_tests::<Fq2>("fq2".to_string());
}

#[test]
fn test_serialization() {
    crate::tests::field::random_serialization_test::<Fq2>("fq2".to_string());
    #[cfg(feature = "derive_serde")]
    crate::tests::field::random_serde_test::<Fq2>("fq2".to_string());
}
