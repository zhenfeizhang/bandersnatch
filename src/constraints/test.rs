use crate::*;
use ark_bls12_381::Bls12_381;
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_groth16::{
    create_random_proof, generate_random_parameters, prepare_verifying_key,
    verify_proof,
};
use ark_r1cs_std::{
    alloc::AllocVar,
    boolean::Boolean,
    eq::EqGadget,
    fields::fp::FpVar,
    groups::{curves::twisted_edwards::AffineVar, CurveVar},
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};

#[test]
fn test_bandersnatch_group_ops() {
    // in this example we are going to argue about the statement
    //      H = xG
    // for some private H, x and G, where
    //  - G is a bandersnatch group element in the Affine form,
    //  - x is a scalar field element
    //  - H is a bandersnatch group element in the Project form,

    let mut rng = ark_std::test_rng();

    // set up the proper inputs for the circuit
    let point_g = EdwardsAffine::rand(&mut rng);
    let x = Fr::rand(&mut rng);
    let point_h = point_g.mul(x).into_affine();
    let circuit = GLVCircuit {
        base: point_g,
        scalar: x,
        res: point_h,
    };

    // circuit sanity check
    let sanity_cs = ConstraintSystem::<Fq>::new_ref();

    circuit
        .clone()
        .generate_constraints(sanity_cs.clone())
        .unwrap();
    assert!(sanity_cs.is_satisfied().unwrap());

    // generate proving key
    let pk = generate_random_parameters::<Bls12_381, _, _>(
        circuit.clone(),
        &mut rng,
    )
    .unwrap();

    // generate the proof
    let proof = create_random_proof(circuit, &pk, &mut rng).unwrap();

    // verify the proof
    let pvk = prepare_verifying_key(&pk.vk);
    assert!(verify_proof(&pvk, &proof, &[]).unwrap())
}

/// a circuit for the relation:
///   res = scalar * base
#[derive(Clone, Debug)]
struct GLVCircuit {
    base: EdwardsAffine,
    scalar: Fr,
    res: EdwardsAffine,
}

impl ConstraintSynthesizer<Fq> for GLVCircuit {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<Fq>,
    ) -> Result<(), SynthesisError> {
        let base_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
            cs.clone(),
            || Ok(self.base),
        )
        .unwrap();

        let mut scalar_bits_var = vec![];
        for e in self.scalar.into_repr().to_bits_le() {
            scalar_bits_var.push(Boolean::new_witness(cs.clone(), || Ok(e))?)
        }

        let res_var_recomputed =
            base_var.scalar_mul_le(scalar_bits_var.iter())?;

        let res_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
            cs.clone(),
            || Ok(self.res),
        )
        .unwrap();

        res_var.enforce_equal(&res_var_recomputed)
    }
}
