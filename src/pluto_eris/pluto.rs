use crate::arithmetic::EndoParameters;
use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::hash_to_curve::svdw_hash_to_curve;
use crate::pluto_eris::fields::{fp::Fp, fq::Fq};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

const PLUTO_GENERATOR_X: Fq = Fq::from_raw([
    0x9ffffcd2ffffffff,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);
const PLUTO_GENERATOR_Y: Fq = Fq::from_raw([
    0x0000000000000007,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
]);

const PLUTO_A: Fq = Fq::from_raw([0, 0, 0, 0, 0, 0, 0]);
const PLUTO_B: Fq = Fq::from_raw([0x39, 0, 0, 0, 0, 0, 0]);

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

impl group::cofactor::CofactorGroup for Pluto {
    type Subgroup = Pluto;

    fn clear_cofactor(&self) -> Self {
        *self
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        1.into()
    }
}

new_curve_impl!(
    (pub),
    Pluto,
    PlutoAffine,
    false,
    Fq,
    Fp,
    (PLUTO_GENERATOR_X,PLUTO_GENERATOR_Y),
    PLUTO_A,
    PLUTO_B,
    "pluto",
    |curve_id, domain_prefix| svdw_hash_to_curve(curve_id, domain_prefix, Pluto::SVDW_Z),
);

impl Pluto {
    const SVDW_Z: Fq = Fq::ONE;
}

//TODO Update comment
// Generated using https://github.com/ConsenSys/gnark-crypto/blob/master/ecc/utils.go
// with `Fp::ZETA`
// See https://github.com/demining/Endomorphism-Secp256k1/blob/main/README.md
// to have more details about the endomorphism.

const ENDO_PARAMS_PLUTO: EndoParameters = EndoParameters {
    // gamma1: 0x2aaaaaaaaaaa955554a0aaaab2aae415b8e5c9500df52222a6a19d565
    gamma1: [
        0xdf52222a6a19d565,
        0x2aae415b8e5c9500,
        0xaaa955554a0aaaab,
        0x00000002aaaaaaaa,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ],
    // gamma2: 0x38e38e38e38e0e38e224e38e4e39d
    gamma2: [
        0xe38e224e38e4e39d,
        0x00038e38e38e38e0,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ],
    // b1: 0x60000000000030000196800005ff0065a001ae513ffffde200000001
    b1: [
        0x3ffffde200000001,
        0x05ff0065a001ae51,
        0x0000300001968000,
        0x0000000060000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ],
    // b2: 0x8000000000002000010f00000000
    b2: [
        0x2000010f00000000,
        0x0000800000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ],
};

endo7!(Pluto, Fp, ENDO_PARAMS_PLUTO);
#[test]
fn test_curve() {
    crate::tests::curve::curve_tests::<Pluto>();
}

#[test]
fn test_hash_to_curve() {
    crate::tests::curve::hash_to_curve_test::<Pluto>();
}

#[test]
fn test_serialization() {
    crate::tests::curve::random_serialization_test::<Pluto>();
    #[cfg(feature = "derive_serde")]
    crate::tests::curve::random_serde_test::<Pluto>();
}

#[test]
fn test_endo_consistency() {
    let g = Pluto::generator();
    assert_eq!(g * Fp::ZETA, g.endo());
}
