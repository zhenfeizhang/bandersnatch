//! This module implements non-native circuit for scalar field operations.
//!

use crate::{
    constraints::{EdwardsAffineVar, FqVar},
    *,
};
use ark_ff::{field_new, FpParameters, PrimeField, Zero};
use ark_r1cs_std::{
    alloc::AllocVar, boolean::Boolean, fields::FieldVar, groups::CurveVar,
    prelude::EqGadget, R1CSVar, ToBitsGadget,
};
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use num_bigint::{BigInt, BigUint};
// use ark_std::println;

/// The lambda parameter for decomposition.
pub(crate) const LAMBDA:Fr = field_new!(Fr, "8913659658109529928382530854484400854125314752504019737736543920008458395397");
/// Lower bits of Lambda, s.t. LAMBDA = LAMBDA_1 + 2^128 LAMBDA_2
const LAMBDA_1: Fq = field_new!(Fq, "276171084342130069873717620938461404933");
/// Higher bits of Lambda, s.t. LAMBDA = LAMBDA_1 + 2^128 LAMBDA_2
const LAMBDA_2: Fq = field_new!(Fq, "26194891433150687747794787867683243819");
/// Lower bits of r, s.t. r = r1 + 2^128 r2
const R1: Fq = field_new!(Fq, "339698375537151026184002201187337562081");
/// Higher bits of r, s.t. r = r1 + 2^128 r2
const R2: Fq = field_new!(Fq, "38523796905489664107205538631725381120");

pub fn glv_mul(
    cs: ConstraintSystemRef<Fq>,
    point_var: &EdwardsAffineVar,
    scalar_var: &[Boolean<Fq>],
) -> Result<EdwardsAffineVar, SynthesisError> {
    let k_vars = decomposition(cs, scalar_var)?;
    let endor_base_var = endomorphism(point_var)?;
    multi_scalar_mul_gadget(point_var, &k_vars[0], &endor_base_var, &k_vars[1])
}

// Here we need to implement a customized MSM algorithm, since we know that
// the high bits of Fr are restricted to be small, i.e. ~ 128 bits.
// This MSM will save us some 128 doublings.
fn multi_scalar_mul_gadget(
    base: &EdwardsAffineVar,
    scalar_1: &[Boolean<Fq>],
    endor_base: &EdwardsAffineVar,
    scalar_2: &[Boolean<Fq>],
) -> Result<EdwardsAffineVar, SynthesisError> {
    let zero = EdwardsAffineVar::zero();

    let endor_base = endor_base.clone().negate()?;

    let sum = base.clone() + endor_base.clone();

    let mut res = EdwardsAffineVar::zero();
    for i in 0..128 {
        res = res.double()?;

        let add = scalar_1[128 - i - 1].select(
            // both bits are 1 => add the sum to self
            // first bit is 1 and second bit is 0 => add the base to self
            &scalar_2[128 - i - 1].select(&sum, base)?,
            // first bit is 0 and second bit is 1 => add the endor_base to self
            // bot bits are 0 => ignore
            &scalar_2[128 - i - 1].select(&endor_base, &zero)?,
        )?;
        res += add;
    }

    Ok(res)
}

pub(crate) fn endomorphism(
    point_var: &EdwardsAffineVar,
) -> Result<EdwardsAffineVar, SynthesisError> {
    let coeff_a1_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_A1);
    let coeff_a2_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_A2);
    let coeff_a3_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_A3);
    let coeff_b1_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_B1);
    let coeff_b2_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_B2);
    let coeff_b3_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_B3);
    let coeff_c1_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_C1);
    let coeff_c2_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_C2);

    let x_var = &point_var.x;
    let y_var = &point_var.y;
    let z_var = y_var;

    let fy_var: FqVar =
        coeff_a1_var * (y_var + coeff_a2_var) * (y_var + coeff_a3_var);
    let gy_var: FqVar =
        coeff_b1_var * (y_var + coeff_b2_var) * (y_var + coeff_b3_var);
    let hy_var: FqVar = (y_var + coeff_c1_var) * (y_var + coeff_c2_var);

    let x_var = x_var * &fy_var * &hy_var;
    let y_var = &gy_var * z_var;
    let z_var = &hy_var * z_var;
    let z_var = z_var.inverse()?;

    let x_var = x_var * &z_var;
    let y_var = y_var * &z_var;

    Ok(EdwardsAffineVar::new(x_var, y_var))
}

// Input a scalar s as in bit wires,
// compute k1 and k2 s.t.
//  s = k1 - lambda * k2 mod |Fr|
// where
// * s ~ 256 bits, private input
// * lambda ~ 256 bits, public input
// * k1, k2 each ~ 128 bits, private inputs
// return the bit wires for k1 and k2
pub(crate) fn decomposition(
    cs: ConstraintSystemRef<Fq>,
    s_vars: &[Boolean<Fq>],
) -> Result<[Vec<Boolean<Fq>>; 2], SynthesisError> {
    assert_eq!(s_vars.len(), 256);

    // for an input scalar s,
    // we need to prove the following statement over ZZ
    // 
    // (0) lambda * k2 + s = t * Fr::modulus + k1
    // 
    // for some t, where
    // * k1, k2 < 2^128
    // * lambda, s, modulus are ~256 bits
    //
    // which becomes
    // (1) lambda_1 * k2 + 2^128 lambda_2 * k2 + s1 + 2^128 s2
    //        - t * r1 - t *2^128 r2 - k1 = 0
    // where
    // (2) lambda = lambda_1 + 2^128 lambda_2   <- public info
    // (3) s = s1 + 2^128 s2
    // (4) Fr::modulus = r1 + 2^128 r2          <- public info
    //
    // reorganizing (1) gives us
    // (5)          lambda_1 * k2 + s1 - t * r1 - k1
    //     + 2^128 (lambda_2 * k2 + s2 - t * r2)
    //     = 0
    //
    // Now set
    // (6) tmp = lambda_1 * k2 + s1 - t * r1 - k1
    // with
    // (7) tmp = tmp1 + 2^128 tmp2
    //
    // that is
    // (8) tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1
    //
    // i.e. tmp2 will be the carrier overflowing 2^128,
    // and on the 2^128 term, we have
    // (9) tmp2 + lambda_2 * k2 + s2  - t * r2 = 0
    //
    // the concrete statements that we need to prove (0) are
    //  (a) k1 < 2^128
    //  (b) k2 < 2^128
    //  (c) s1 < 2^128
    //  (d) s2 < 2^128
    //  (e) tmp1 < 2^128
    //  (f) tmp2 < 2^128
    //  (h) s = s1 + 2^128 s2
    //  (i) tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1
    //  (j) tmp2 + lambda_2 * k2 + s2  = t * r2
    // which can be evaluated over Fq without overflow

    // ============================================
    // step 1: build integers
    // ============================================
    // 2^128
    let two_to_128 = BigInt::from(2u64).pow(128);

    // s = s1 + 2^128 * s2
    let s = Boolean::le_bits_to_fp_var(s_vars)?.value()?;
    let s_int = fq_to_big_int!(s);
    let s_fr = Fr::from(s.into_repr());
    let s1_int = &s_int % &two_to_128;
    let s2_int = &s_int / &two_to_128;

    // lambda = lambda_1 + 2^128 lambda_2
    let lambda_int = fq_to_big_int!(LAMBDA);
    let lambda_1_int = fq_to_big_int!(LAMBDA_1);
    let lambda_2_int = fq_to_big_int!(LAMBDA_2);

    // s = k1 + lambda * k2
    let (k1, k2) = BandersnatchParameters::scalar_decomposition(&s_fr);
    let k2 = -k2;
    let k1_int = fq_to_big_int!(k1);
    let k2_int = fq_to_big_int!(k2);

    // fr_order = r1 + 2^128 r2
    let fr_order_uint: BigUint = <Fr as PrimeField>::Params::MODULUS.into();
    let fr_order_int: BigInt = fr_order_uint.into();
    let r1_int = fq_to_big_int!(R1);
    let r2_int = fq_to_big_int!(R2);

    // t = (lambda * k2 + s - k1) / fr_order
    let t_int = (&lambda_int * &k2_int + &s_int - &k1_int) / &fr_order_int;

    // tmp = tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1
    let tmp_int =
        &lambda_1_int * &k2_int + &s1_int - &t_int * &r1_int - &k1_int;
    let tmp1_int = &tmp_int % &two_to_128;
    let tmp2_int = &tmp_int / &two_to_128;

    // step 1.1 check the correctness of (0), (5), (i) and (j) in the clear
    // equation (0): lambda * k2 + s = t * Fr::modulus + k1
    assert_eq!(
        &s_int + &lambda_int * &k2_int,
        &k1_int + &t_int * &fr_order_int
    );

    // equation (5)
    //              lambda_1 * k2 + s1 - t * r1 - k1
    //     + 2^128 (lambda_2 * k2 + s2 - t * r2)
    //     = 0
    assert_eq!(
        &lambda_1_int * &k2_int + &s1_int - &t_int * &r1_int - &k1_int
            + &two_to_128
                * (&lambda_2_int * &k2_int + &s2_int - &t_int * &r2_int),
        BigInt::zero()
    );

    // equation (i): tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1
    assert_eq!(
        &tmp1_int + &two_to_128 * &tmp2_int,
        &lambda_1_int * &k2_int + &s1_int - &t_int * &r1_int - &k1_int
    );

    //  equation (j) tmp2 + lambda_2 * k2 + s2  = t * r2
    assert_eq!(
        &tmp2_int + &lambda_2_int * &k2_int + &s2_int,
        &t_int * &r2_int
    );

    // ============================================
    // step 2. build the variables
    // ============================================

    // constant variables
    let lambda_1_var = FqVar::new_constant(cs.clone(), LAMBDA_1)?;
    let lambda_2_var = FqVar::new_constant(cs.clone(), LAMBDA_2)?;

    let r1_var = FqVar::new_constant(cs.clone(), R1)?;
    let r2_var = FqVar::new_constant(cs.clone(), R2)?;

    // secret variables
    let k1_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k1_int)))?;
    let k2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k2_int)))?;

    let s1_var = convert_boolean_var_to_fp_var(&s_vars[0..128])?;
    let s2_var = convert_boolean_var_to_fp_var(&s_vars[128..])?;

    let t_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(t_int)))?;

    let tmp_var = FqVar::new_witness(cs, || Ok(int_to_fq!(tmp_int)))?;
    let tmp_boolean_vars = tmp_var.to_bits_le()?;
    let tmp2_var = convert_boolean_var_to_fp_var(&tmp_boolean_vars[128..])?;

    // ============================================
    // step 3. range proofs
    // ============================================
    //  (a) k1 < 2^128
    //  (b) k2 < 2^128
    let k1_bits_var = k1_var.to_bits_le()?;
    let k2_bits_var = k2_var.to_bits_le()?;
    enforce_bits_are_0s(&[&k1_bits_var[128..], &k2_bits_var[128..]].concat())?;

    //  (c) s1 < 2^128
    //  (d) s2 < 2^128
    // Implied by
    //  let s1_var = convert_boolean_var_to_fp_var(&s_vars[0..128])?;
    //  let s2_var = convert_boolean_var_to_fp_var(&s_vars[128..])?;

    //  (e) tmp1 < 2^128
    //  (f) tmp2 < 2^128
    //  (h) tmp = tmp1 + 2^128 tmp2
    // Implied by
    //  let tmp2_var = convert_boolean_var_to_fp_var( &tmp_boolean_vars[128..])?;

    // ============================================
    // step 4. equality proofs
    // ============================================
    //  (i) tmp =  lambda_1 * k2 + s1 - t * r1 - k1
    let right1 = lambda_1_var * (&k2_var);
    let right2 = right1 + (&s1_var);
    let right3 = (&t_var) * (&r1_var);
    let right = right2 - (&right3);
    let right = right - (&k1_var);
    right.enforce_equal(&tmp_var)?;

    //  (j) tmp2 + lambda_2 * k2 + s2 = t * r2
    let left1 = lambda_2_var * (&k2_var);
    let left = tmp2_var + (&left1);
    let left = left + (&s2_var);
    let right = t_var * (r2_var);
    right.enforce_equal(&left)?;

    // extract the output
    Ok([k1_bits_var[..128].to_vec(), k2_bits_var[..128].to_vec()])
}

#[macro_export]
macro_rules! fq_to_big_int {
    ($fq: expr) => {
        <BigInt as From<BigUint>>::from($fq.into_repr().into())
    };
}

#[macro_export]
macro_rules! int_to_fq {
    ($in: expr) => {
        Fq::from_le_bytes_mod_order(&$in.to_bytes_le().1)
    };
}

pub(crate) fn convert_boolean_var_to_fp_var(
    input: &[Boolean<Fq>],
) -> Result<FqVar, SynthesisError> {
    Boolean::le_bits_to_fp_var(input)
}


/// For a slice of bits: b_0...b_k, enforce
///     b_0^2 + b_1^2 ... + b_k^0 = 0 over Fq
/// This is slightly more efficient than Boolean::kary_or()
/// which enforces b0 | b1 | ... | bk = 0
fn enforce_bits_are_0s(input: &[Boolean<Fq>]) -> Result<(), SynthesisError> {
    let bits_fq: Vec<FqVar> = input
        .iter()
        .take(128)
        .map(|x| FqVar::from(x.clone()))
        .collect();

    let mut sum = bits_fq[0].clone() * (&bits_fq[0].clone());
    for e in bits_fq.iter().take(128).skip(1) {
        sum += &e.clone() * &e.clone();
    }
    sum.enforce_equal(&FqVar::zero())
}
