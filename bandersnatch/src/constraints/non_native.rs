//! This module implements non-native circuit for scalar field operations.
//!

use crate::{constraints::FqVar, Fq, Fr, *};
use ark_ff::{field_new, FpParameters, PrimeField, Zero};
use ark_r1cs_std::{alloc::AllocVar, prelude::EqGadget};
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

// Constrain s = k1 - lambda * k2 mod |Fr|
// where
// * s ~ 256 bits, private input
// * lambda ~ 256 bits, public input
// * k1, k2 each ~ 128 bits, private inputs
pub(crate) fn decomposition(
    k1: Fr,
    k2: Fr,
    s: Fr,
    cs: ConstraintSystemRef<Fq>,
) -> Result<(), SynthesisError> {
    // check the correctness of the inputs
    assert_eq!(s, k1 - k2 * LAMBDA);

    // we need to prove the following statement over Fq
    // (0) lambda * k2 + s = t * Fr::modulus + k1
    //
    // where
    // * k1, k2 < 2^128
    // * t < 2^128
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
    // (6) tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1
    //
    // i.e. tmp2 will be the carrier overflowing 2^128.
    // we know that
    // (7) tmp2 + lambda_2 * k2 + s2  - t * r2 = 0
    //
    // the concrete statements that we need to prove are
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
    let two_to_128 = BigInt::from(2u64).pow(128);

    let fr_order_uint: BigUint = <Fr as PrimeField>::Params::MODULUS.into();
    let fr_order_int: BigInt = fr_order_uint.into();

    let k1_int = fq_to_big_int!(k1);
    let k2_int = fq_to_big_int!(k2);

    let lambda_int = fq_to_big_int!(LAMBDA);
    let lambda_1_int = fq_to_big_int!(LAMBDA_1);
    let lambda_2_int = fq_to_big_int!(LAMBDA_2);
    assert_eq!(lambda_int, &lambda_1_int + &lambda_2_int * &two_to_128);

    let s_int = fq_to_big_int!(s);
    let s1_int = &s_int % &two_to_128;
    let s2_int = &s_int / &two_to_128;
    assert_eq!(s_int, &s1_int + &s2_int * &two_to_128);

    let t_int = (&lambda_int * &k2_int + &s_int - &k1_int) / &fr_order_int;

    let r1_int = fq_to_big_int!(R1);
    let r2_int = fq_to_big_int!(R2);

    let tmp = &lambda_1_int * &k2_int + &s1_int - &t_int * &r1_int - &k1_int;
    let tmp1_int = &tmp % &two_to_128;
    let tmp2_int = &tmp / &two_to_128;

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
    let two_to_128_var = FqVar::new_constant(
        cs.clone(),
        Fq::from(BigUint::from(2u64).pow(128)),
    )?;
    let lambda_1_var = FqVar::new_constant(cs.clone(), LAMBDA_1)?;
    let lambda_2_var = FqVar::new_constant(cs.clone(), LAMBDA_2)?;

    let r1_var = FqVar::new_constant(cs.clone(), R1)?;
    let r2_var = FqVar::new_constant(cs.clone(), R2)?;

    // secret variables
    let k1_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k1_int)))?;
    let k2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(k2_int)))?;

    let s_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(s_int)))?;
    let s1_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(s1_int)))?;
    let s2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(s2_int)))?;

    let t_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(t_int)))?;

    let tmp1_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(tmp1_int)))?;
    let tmp2_var = FqVar::new_witness(cs.clone(), || Ok(int_to_fq!(tmp2_int)))?;

    // ============================================
    // step 3. range proofs
    // ============================================
    //  (a) k1 < 2^128
    //  (b) k2 < 2^128

    k1_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;
    k2_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;

    //  (c) s1 < 2^128
    //  (d) s2 < 2^128

    s1_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;
    s2_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;

    //  (e) tmp1 < 2^128
    //  (f) tmp2 < 2^128
    tmp1_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;
    tmp2_var.enforce_cmp(&two_to_128_var, std::cmp::Ordering::Less, false)?;

    // ============================================
    // step 4. equality proofs
    // ============================================
    //  (h) s = s1 + 2^128 s2

    let right = (&s2_var) * (&two_to_128_var);
    let right = (&right) + (&s1_var);
    right.enforce_equal(&s_var)?;

    //  (i) tmp1 + 2^128 tmp2 =  lambda_1 * k2 + s1 - t * r1 - k1

    let left = (&tmp2_var) * (&two_to_128_var);
    let left = left + (&tmp1_var);
    let right1 = lambda_1_var * (&k2_var);
    let right2 = right1 + (&s1_var);
    let right3 = (&t_var) * (&r1_var);
    let right = right2 - (&right3);
    let right = right - (&k1_var);
    right.enforce_equal(&left)?;

    //  (j) tmp2 + lambda_2 * k2 + s2 = t * r2
    let left1 = lambda_2_var * (&k2_var);
    let left = tmp2_var + (&left1);
    let left = left + (&s2_var);
    let right = t_var * (r2_var);
    right.enforce_equal(&left)?;

    Ok(())
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
