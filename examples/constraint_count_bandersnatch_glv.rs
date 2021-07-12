use std::marker::PhantomData;

use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_r1cs_std::{
    alloc::AllocVar,
    boolean::Boolean,
    eq::EqGadget,
    fields::{fp::FpVar, FieldVar},
    groups::{
        curves::{
            short_weierstrass::ProjectiveVar, twisted_edwards::AffineVar,
        },
        CurveVar,
    },
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};
use bandersnatch::{
    constraints::FqVar, EdwardsAffine, EdwardsParameters, Fq, Fr, GLVParameters,
};
// use ark_std::marker::PhantomData;

fn main() {
    // in this example we are going to argue about the statement
    //      H = xG
    // for some private H, x and G, where
    //  - G is a bandersnatch group element in the Affine form,
    //  - x is a scalar field element
    //  - H is a bandersnatch group element in the Affine form,

    let mut rng = ark_std::test_rng();

    let point_g = EdwardsAffine::rand(&mut rng);
    let x = Fr::rand(&mut rng);
    let point_h = point_g.mul(x).into_affine();
    let circuit = GroupOpCircuit {
        base: point_g,
        scalar: x,
        res: point_h,
    };

    let sanity_cs = ConstraintSystem::<Fq>::new_ref();

    circuit.generate_constraints(sanity_cs.clone()).unwrap();
    assert!(sanity_cs.is_satisfied().unwrap());
}

/// a circuit for the relation:
///   res = scalar * base
struct GroupOpCircuit {
    base: EdwardsAffine,
    scalar: Fr,
    res: EdwardsAffine,
}

impl ConstraintSynthesizer<Fq> for GroupOpCircuit {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<Fq>,
    ) -> Result<(), SynthesisError> {
        let _cs_no = cs.num_constraints();
        let base_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
            cs.clone(),
            || Ok(self.base),
        )
        .unwrap();
        #[cfg(debug_assertions)]
        println!("cs for base var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let phi_var = endomorphism_gadget(&base_var);

        #[cfg(debug_assertions)]
        println!("cs for endomorphism var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let mut scalar_bits_var = vec![];
        for e in self.scalar.into_repr().to_bits_le() {
            scalar_bits_var.push(Boolean::new_witness(cs.clone(), || Ok(e))?)
        }

        #[cfg(debug_assertions)]
        println!("cs for scalar var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let res_var_recomputed =
            base_var.scalar_mul_le(scalar_bits_var.iter())?;

        #[cfg(debug_assertions)]
        println!("cs for mul : {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let res_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
            cs.clone(),
            || Ok(self.res),
        )
        .unwrap();

        #[cfg(debug_assertions)]
        println!("cs for result var : {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        res_var.enforce_equal(&res_var_recomputed)?;

        #[cfg(debug_assertions)]
        println!("cs for equality : {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        Ok(())
    }
}

fn endomorphism_gadget(
    base_point: &AffineVar<EdwardsParameters, FqVar>,
    // cs: ConstraintSystemRef<Fq>,
) -> AffineVar<EdwardsParameters, FqVar> {
    let coeff_a1_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A1);
    let coeff_a2_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A2);
    let coeff_a3_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A3);
    let coeff_b1_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B1);
    let coeff_b2_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B2);
    let coeff_b3_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B3);
    let coeff_c1_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_C1);
    let coeff_c2_var =
        FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_C2);

    let x_var = base_point.x.clone();
    let y_var = base_point.y.clone();
    let z_var = y_var.clone();

    let fy_var: FqVar = coeff_a1_var
        * (y_var.clone() + coeff_a2_var)
        * (y_var.clone() + coeff_a3_var);
    let gy_var: FqVar = coeff_b1_var
        * (y_var.clone() + coeff_b2_var)
        * (y_var.clone() + coeff_b3_var);
    let hy_var: FqVar =
        (y_var.clone() + coeff_c1_var) * (y_var.clone() + coeff_c2_var);

    let x_var = x_var * fy_var * hy_var.clone();
    let y_var = gy_var * z_var.clone();
    let z_var = hy_var * z_var;
    let z_var = z_var.inverse().unwrap();

    let x_var = x_var * z_var.clone();
    let y_var = y_var * z_var;

    AffineVar::new(x_var, y_var)
}
