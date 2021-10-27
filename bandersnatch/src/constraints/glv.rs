//! This module implements non-native circuit for scalar field operations.
//!

use crate::{
    constraints::{EdwardsAffineVar, FqVar},
    *,
};
use ark_ec::ProjectiveCurve;
use ark_ff::{field_new, BigInteger, FpParameters, PrimeField, Zero};
use ark_r1cs_std::{
    alloc::AllocVar, boolean::Boolean, fields::FieldVar, groups::CurveVar,
    prelude::EqGadget, R1CSVar,
};
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use num_bigint::{BigInt, BigUint};

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

pub fn glv_mul_gadget(
    cs: ConstraintSystemRef<Fq>,
    point_var: &EdwardsAffineVar,
    scalar_var: &FqVar,
) -> Result<EdwardsAffineVar, SynthesisError> {
    #[cfg(debug_assertions)]
    println!(
        "number of constraints before decomposition: {}",
        cs.num_constraints()
    );

    let k_vars = scalar_decomposition_gadget(cs.clone(), scalar_var)?;

    #[cfg(debug_assertions)]
    println!(
        "number of constraints after decomposition: {}",
        cs.num_constraints()
    );

    let endor_base_var = endomorphism_gadget(cs.clone(), point_var)?;

    #[cfg(debug_assertions)]
    println!(
        "number of constraints after endomorphism: {}",
        cs.num_constraints()
    );

    multi_scalar_mul_gadget(point_var, &k_vars[0], &endor_base_var, &k_vars[1])
}

/// The circuit for 2 base scalar multiplication.
pub fn multi_scalar_mul_gadget(
    base: &EdwardsAffineVar,
    scalar_1: &[Boolean<Fq>],
    endor_base: &EdwardsAffineVar,
    scalar_2: &[Boolean<Fq>],
) -> Result<EdwardsAffineVar, SynthesisError> {
    let length = scalar_1.len();
    assert_eq!(
        length,
        scalar_2.len(),
        "the input length do not equal: {} vs {}",
        length,
        scalar_2.len()
    );

    let zero = EdwardsAffineVar::zero();

    let endor_base = endor_base.clone().negate()?;

    let sum = base.clone() + endor_base.clone();

    let mut res = EdwardsAffineVar::zero();
    for i in 0..length {
        res = res.double()?;

        let add = scalar_1[length - i - 1].select(
            // both bits are 1 => add the sum to self
            // first bit is 1 and second bit is 0 => add the base to self
            &scalar_2[length - i - 1].select(&sum, base)?,
            // first bit is 0 and second bit is 1 => add the endor_base to self
            // bot bits are 0 => ignore
            &scalar_2[length - i - 1].select(&endor_base, &zero)?,
        )?;
        res += add;
    }

    Ok(res)
}

/// The circuit for computing the point endomorphism.
pub fn endomorphism_gadget(
    cs: ConstraintSystemRef<Fq>,
    point_var: &EdwardsAffineVar,
) -> Result<EdwardsAffineVar, SynthesisError> {
    let base = point_var.value()?.into_affine();
    let endor_point =
        <BandersnatchParameters as GLVParameters>::endomorphism(&base);
    let new_x_var = FqVar::new_witness(cs.clone(), || Ok(endor_point.x))?;
    let new_y_var = FqVar::new_witness(cs.clone(), || Ok(endor_point.y))?;

    let coeff_b_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_B);
    let coeff_c_var =
        &FqVar::constant(<BandersnatchParameters as GLVParameters>::COEFF_C);

    let coeff_b_square_var = &FqVar::constant(
        <BandersnatchParameters as GLVParameters>::COEFF_B
            * <BandersnatchParameters as GLVParameters>::COEFF_B,
    );

    let x_var = &point_var.x;
    let y_var = &point_var.y;

    // xy = x * y
    let xy_var = x_var * y_var;
    // y_square = y^2
    let y_square_var = y_var * y_var;

    // f(y) = c(1-y^2)
    let f_y_var = coeff_c_var - coeff_c_var * &y_square_var;
    // g(y) = b(y^2 + b)
    let g_y_var = coeff_b_var * &y_square_var + coeff_b_square_var;
    // h(y) = y^2 - b
    let h_y_var = y_square_var - coeff_b_var;

    // x = f(y) / (xy)
    f_y_var.enforce_equal(&(&xy_var * &new_x_var))?;
    // y = g(y)/ h(y)
    g_y_var.enforce_equal(&(&h_y_var * &new_y_var))?;

    Ok(EdwardsAffineVar::new(new_x_var, new_y_var))
}

// Input a scalar s as in bit wires,
// compute k1 and k2 s.t.
//  s = k1 - lambda * k2 mod |Fr|
// where
// * s ~ 253 bits, private input
// * lambda ~ 253 bits, public input
// * k1, k2 each ~ 128 bits, private inputs
// return the bit wires for k1 and k2
pub fn scalar_decomposition_gadget(
    cs: ConstraintSystemRef<Fq>,
    s_var: &FqVar,
) -> Result<[Vec<Boolean<Fq>>; 2], SynthesisError> {
    // the oder of scalar field
    // r = 13108968793781547619861935127046491459309155893440570251786403306729687672801 < 2^253
    // q = 52435875175126190479447740508185965837690552500527637822603658699938581184513 < 2^255

    // for an input scalar s,
    // we need to prove the following statement over ZZ
    //
    // (0) lambda * k2 + s = t * Fr::modulus + k1
    //
    // for some t, where
    // * k1, k2 < sqrt{2r} < 2^127
    // * lambda, s, modulus are ~253 bits
    //
    // which becomes
    // (1) lambda_1 * k2 + 2^128 lambda_2 * k2 + s
    //        - t * r1 - t *2^128 r2 - k1 = 0
    // where
    // (2) lambda = lambda_1 + 2^128 lambda_2   <- public info
    // (3) Fr::modulus = r1 + 2^128 r2          <- public info
    // with
    //  lambda_1 and r1 < 2^128
    //  lambda_2 and r2 < 2^125
    //
    // reorganizing (1) gives us
    // (4)          lambda_1 * k2 + s - t * r1 - k1
    //     + 2^128 (lambda_2 * k2 - t * r2)
    //     = 0
    //
    // Now set
    // (5) tmp = lambda_1 * k2 + s - t * r1 - k1
    // with
    // (6) tmp = tmp1 + 2^128 tmp2
    // for tmp1 < 2^128 and tmp2 < 2^128
    //
    // that is
    // (7) tmp1 =  (lambda_1 * k2 + s - t * r1 - k1) % 2^128 = 0
    //
    // i.e. tmp2 will be the carrier overflowing 2^128,
    // and on the 2^128 term, we have
    // (8) tmp2 + lambda_2 * k2 - t * r2 = 0
    //
    // the concrete statements that we need to prove (0) are
    //  (a) k1 < 2^127
    //  (b) k2 < 2^127
    //  (c) tmp1 = 0
    //  (d) tmp2 < 2^128
    //  (e) tmp = tmp1 + 2^128 tmp2
    //  (f) tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s - t * r1 - k1
    //  (g) tmp2 + lambda_2 * k2   = t * r2
    // which can all be evaluated over Fq without overflow

    // ============================================
    // step 1: build integers
    // ============================================
    // 2^128
    let two_to_128 = BigInt::from(2u64).pow(128);

    // s
    let s = s_var.value()?;
    let s_int = fq_to_big_int!(s);
    let s_fr = Fr::from(s.into_repr());

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

    // tmp = tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s - t * r1 - k1
    let tmp_int = &lambda_1_int * &k2_int + &s_int - &t_int * &r1_int - &k1_int;
    let tmp1_int = &tmp_int % &two_to_128;
    let tmp2_int = &tmp_int / &two_to_128;

    // step 1.1 check the correctness of (0), (4), (f) and (g) in the clear
    // equation (0): lambda * k2 + s = t * Fr::modulus + k1
    assert_eq!(
        &s_int + &lambda_int * &k2_int,
        &k1_int + &t_int * &fr_order_int
    );

    // equation (4)
    //              lambda_1 * k2 + s - t * r1 - k1
    //     + 2^128 (lambda_2 * k2 - t * r2)
    //     = 0
    assert_eq!(
        &lambda_1_int * &k2_int + &s_int - &t_int * &r1_int - &k1_int
            + &two_to_128 * (&lambda_2_int * &k2_int - &t_int * &r2_int),
        BigInt::zero()
    );

    // equation (f): tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s - t * r1 - k1
    assert_eq!(
        &tmp1_int + &two_to_128 * &tmp2_int,
        &lambda_1_int * &k2_int + &s_int - &t_int * &r1_int - &k1_int
    );

    // equation (g) tmp2 + lambda_2 * k2 + s2  = t * r2
    assert_eq!(&tmp2_int + &lambda_2_int * &k2_int, &t_int * &r2_int);

    // ============================================
    // step 2. build the variables
    // ============================================

    // constant variables
    let lambda_1_var = FqVar::new_constant(cs.clone(), LAMBDA_1)?;
    let lambda_2_var = FqVar::new_constant(cs.clone(), LAMBDA_2)?;

    let r1_var = FqVar::new_constant(cs.clone(), R1)?;
    let r2_var = FqVar::new_constant(cs.clone(), R2)?;

    let two_to_128_var = FqVar::new_constant(
        cs.clone(),
        Fq::from(BigUint::from(2u64).pow(128)),
    )?;

    // secret variables
    let k1_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k1_int)))?;
    let k2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k2_int)))?;

    let t_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(t_int)))?;

    let tmp_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(tmp_int)))?;
    let tmp2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(tmp2_int)))?;

    // ============================================
    // step 3. range proofs
    // ============================================
    //  (a) k1 < 2^127
    //  (b) k2 < 2^127
    let k1_bits_vars =
        decompose_and_enforce_less_than_k_bits(cs.clone(), &k1_var, 127)?;
    let k2_bits_vars =
        decompose_and_enforce_less_than_k_bits(cs.clone(), &k2_var, 127)?;

    //  (c) tmp1 = 0        <- implied by tmp = 2^128 * tmp2
    //  (d) tmp2 < 2^128
    //  (e) tmp = tmp1 + 2^128 tmp2
    let tmp_var_rec = &tmp2_var * &two_to_128_var;
    tmp_var.enforce_equal(&tmp_var_rec)?;
    decompose_and_enforce_less_than_k_bits(cs, &tmp2_var, 128)?;

    // ============================================
    // step 4. equality proofs
    // ============================================
    //  (f) tmp =  lambda_1 * k2 + s - t * r1 - k1
    let right1 = lambda_1_var * (&k2_var);
    let right2 = right1 + s_var;
    let right3 = (&t_var) * (&r1_var);
    let right = right2 - (&right3);
    let right = right - (&k1_var);
    right.enforce_equal(&tmp_var)?;

    //  (g) tmp2 + lambda_2 * k2 = t * r2
    let left1 = lambda_2_var * (&k2_var);
    let left = tmp2_var + (&left1);
    let right = t_var * (r2_var);
    right.enforce_equal(&left)?;

    // extract the output
    Ok([k1_bits_vars, k2_bits_vars])
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

// Decomposite the elements into bit_length booleans and prove that
// input is the composition of those bit_length booleans.
// This implies that input < 2^bit_length.
// Return the bit_length boolean constrains.
// Cost: bit_length+1 constraints
fn decompose_and_enforce_less_than_k_bits(
    cs: ConstraintSystemRef<Fq>,
    input_var: &FqVar,
    bit_length: usize,
) -> Result<Vec<Boolean<Fq>>, SynthesisError> {
    if bit_length == 0 {
        panic!("invalid input bit length {}", bit_length)
    }

    let input_val = input_var.value()?;
    let input_bits = input_val.into_repr().to_bits_le();
    let input_bits_vars = input_bits
        .iter()
        .take(bit_length)
        .map(|b| Boolean::new_witness(cs.clone(), || Ok(*b)))
        .collect::<Result<Vec<_>, _>>()?;

    let mut res_var = FqVar::from(input_bits_vars[bit_length - 1].clone());
    for e in input_bits_vars.iter().rev().skip(1) {
        res_var = res_var.double()?;
        res_var = &res_var + FqVar::from(e.clone());
    }
    res_var.enforce_equal(input_var)?;

    Ok(input_bits_vars)
}
