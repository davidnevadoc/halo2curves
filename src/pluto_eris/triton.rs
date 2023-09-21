use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::pluto_eris::fields::{fp::Fp, fq::Fq};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use group::cofactor::CofactorGroup;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

const TRITON_GENERATOR_X: Fq2 = Fq2 {
    // 0x13576c81faf3a13fd815d0e9bd54b845ee935948b84498b27ca972bfb93722e223c9e276a4ebe7559cfc86dd865f07d64f2b5fe6556f9066
    c0: Fq::from_raw([
        0x4f2b5fe6556f9066,
        0x9cfc86dd865f07d6,
        0x23c9e276a4ebe755,
        0x7ca972bfb93722e2,
        0xee935948b84498b2,
        0xd815d0e9bd54b845,
        0x13576c81faf3a13f,
    ]),

    //0x142164cb875db0465e5092f9380f44f555243d011699b7393029f2d201554727aeb383298fdf5847b9b3dff01bbe8d63fe7c781a8fd7bf21
    c1: Fq::from_raw([
        0xfe7c781a8fd7bf21,
        0xb9b3dff01bbe8d63,
        0xaeb383298fdf5847,
        0x3029f2d201554727,
        0x55243d011699b739,
        0x5e5092f9380f44f5,
        0x142164cb875db046,
    ]),
};
const TRITON_GENERATOR_Y: Fq2 = Fq2 {
    //0x2239f7408ead478c58e88d4df1e7418c42fdbb92e64ba85aa4dc17d7dace3f32eb471c004db774bfe78574aca67b3898cd1b78ad106ab9fe
    c0: Fq::from_raw([
        0xcd1b78ad106ab9fe,
        0xe78574aca67b3898,
        0xeb471c004db774bf,
        0xa4dc17d7dace3f32,
        0x42fdbb92e64ba85a,
        0x58e88d4df1e7418c,
        0x2239f7408ead478c,
    ]),

    // 0x1260b04d51136590dbb53dfd7caf450aeca714555bbe4f079ca65d97eb28fc9fc697b4e10bbcd9e0539ef82a731fb88ed49e3c080e6d945d
    c1: Fq::from_raw([
        0xd49e3c080e6d945d,
        0x539ef82a731fb88e,
        0xc697b4e10bbcd9e0,
        0x9ca65d97eb28fc9f,
        0xeca714555bbe4f07,
        0xdbb53dfd7caf450a,
        0x1260b04d51136590,
    ]),
};

// u + 3
const TRITON_B: Fq2 = Fq2 {
    c0: Fq::from_raw([0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
    c1: Fq::ONE,
};

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

use super::fields::fq2::Fq2;

impl group::cofactor::CofactorGroup for Triton {
    type Subgroup = Triton;

    fn clear_cofactor(&self) -> Self {
        // cofactor = 2*q - p
        //0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5d3a8a6c7be4a7d5fe91447fd6a8a7e928a00867971ffffcd300000001
        let e: [u8; 56] = [
            0x24, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x00, 0x01, 0x30, 0xe0, 0x00, 0x0d, 0x7f,
            0x70, 0xe4, 0xa8, 0x03, 0xca, 0x76, 0xf4, 0x39, 0x26, 0x6f, 0x44, 0x3f, 0x9a, 0x5d,
            0x3a, 0x8a, 0x6c, 0x7b, 0xe4, 0xa7, 0xd5, 0xfe, 0x91, 0x44, 0x7f, 0xd6, 0xa8, 0xa7,
            0xe9, 0x28, 0xa0, 0x08, 0x67, 0x97, 0x1f, 0xff, 0xfc, 0xd3, 0x00, 0x00, 0x00, 0x01,
        ];

        // self * TRITON_COFACTOR
        let mut acc = Triton::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = Triton::conditional_select(&acc, &(acc + self), bit);
        }
        acc
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        unimplemented!()
    }

    fn is_torsion_free(&self) -> Choice {
        //
        let e: [u8; 56] = [
            0x24, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x00, 0x01, 0x30, 0xe0, 0x00, 0x0d, 0x7f,
            0x70, 0xe4, 0xa8, 0x03, 0xca, 0x76, 0xf4, 0x39, 0x26, 0x6f, 0x44, 0x3f, 0x9a, 0x5c,
            0xda, 0x8a, 0x6c, 0x7b, 0xe4, 0xa7, 0xa5, 0xfe, 0x8f, 0xad, 0xff, 0xd6, 0xa2, 0xa7,
            0xe8, 0xc3, 0x00, 0x06, 0xb9, 0x45, 0x9f, 0xff, 0xfc, 0xd3, 0x00, 0x00, 0x00, 0x01,
        ];
        // self * GROUP_ORDER;
        let mut acc = Triton::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = Triton::conditional_select(&acc, &(acc + self), bit);
        }
        acc.is_identity()
    }
}

new_curve_impl!(
    (pub),
    Triton,
    TritonAffine,
    false,
    Fq2,
    Fp,
    (TRITON_GENERATOR_X,TRITON_GENERATOR_Y),
    TRITON_B,
    "triton",
);

#[test]
fn test_generator() {
    let gen = Triton::generator().clear_cofactor();
    assert!(bool::from(gen.is_on_curve()))
}

#[test]
fn test_curve() {
    crate::tests::curve::curve_tests::<Triton>();
}

#[test]
fn test_serialization() {
    crate::tests::curve::random_serialization_test::<Triton>();
    #[cfg(feature = "derive_serde")]
    crate::tests::curve::random_serde_test::<Triton>();
}

#[test]
fn test_endo_consistency() {
    let g = Triton::generator();
    assert_eq!(g * Fp::ZETA, g.endo());
}
