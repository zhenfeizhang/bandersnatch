use super::glv::*;
use crate::{constraints::*, *};
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
    R1CSVar, ToBitsGadget,
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};
use ark_std::rand::Rng;

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
        let base_var =
            AffineVar::<BandersnatchParameters, FpVar<Fq>>::new_witness(
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

        let res_var =
            AffineVar::<BandersnatchParameters, FpVar<Fq>>::new_witness(
                cs.clone(),
                || Ok(self.res),
            )
            .unwrap();

        res_var.enforce_equal(&res_var_recomputed)
    }
}

#[test]
fn test_non_native_decomposition() {
    let mut rng = ark_std::test_rng();
    for _ in 0..10 {
        let scalar: Fr = Fr::rand(&mut rng);
        let (k1, k2) = BandersnatchParameters::scalar_decomposition(&scalar);
        assert_eq!(k2 * LAMBDA + k1, scalar);

        let cs = ConstraintSystem::<Fq>::new_ref();
        let scalar_var =
            FqVar::new_witness(cs.clone(), || Ok(Fq::from(scalar.into_repr())))
                .unwrap();
        let _k_vars =
            scalar_decomposition_gadget(cs.clone(), &scalar_var).unwrap();
        println!(
            "number of constraints for decomposition: {}",
            cs.num_constraints()
        );
        assert!(cs.is_satisfied().unwrap());
    }
    // assert!(false)
}

#[test]
fn test_boolean_fp_var_conversion() {
    let mut rng = ark_std::test_rng();
    for _ in 0..10 {
        let element: Fq = Fq::rand(&mut rng);

        let cs = ConstraintSystem::<Fq>::new_ref();

        let element_var =
            FpVar::new_witness(cs.clone(), || Ok(element)).unwrap();
        let boolean_vars = element_var.to_bits_le().unwrap();

        println!(
            "number of constraints for decomposition: {}",
            cs.num_constraints()
        );
        println!(
            "number of variables for decomposition: {} {}",
            cs.num_instance_variables(),
            cs.num_witness_variables(),
        );
        assert!(cs.is_satisfied().unwrap());
        assert_eq!(element_var.value().unwrap(), element);

        let new_element_var =
            convert_boolean_var_to_fp_var(boolean_vars.as_ref()).unwrap();

        assert!(cs.is_satisfied().unwrap());
        assert_eq!(new_element_var.value().unwrap(), element);

        println!(
            "number of constraints for decomposition: {}",
            cs.num_constraints()
        );
        println!(
            "number of variables for decomposition: {} {}",
            cs.num_instance_variables(),
            cs.num_witness_variables(),
        );

        new_element_var.enforce_equal(&element_var).unwrap();

        assert!(cs.is_satisfied().unwrap());
        assert_eq!(new_element_var.value().unwrap(), element);

        println!(
            "number of constraints for decomposition: {}",
            cs.num_constraints()
        );
        println!(
            "number of variables for decomposition: {} {}\n",
            cs.num_instance_variables(),
            cs.num_witness_variables(),
        );
    }
    // assert!(false)
}

#[test]
fn test_glv_circuit() {
    let mut rng = ark_std::test_rng();
    for _ in 0..10 {
        let cs = ConstraintSystem::<Fq>::new_ref();
        let scalar: Fr = Fr::rand(&mut rng);
        let point: EdwardsAffine = rng.gen();
        let res =
            BandersnatchParameters::glv_mul(&point, &scalar).into_affine();

        // println!("{}", scalar.into_repr());

        let point_var =
            EdwardsAffineVar::new_witness(cs.clone(), || Ok(point)).unwrap();

        let scalar_var =
            FqVar::new_witness(cs.clone(), || Ok(Fq::from(scalar.into_repr())))
                .unwrap();

        println!("number of constraints before glv: {}", cs.num_constraints());

        let res_var =
            glv_mul_gadget(cs.clone(), &point_var, &scalar_var).unwrap();

        println!("number of constraints after glv: {}", cs.num_constraints());
        assert!(cs.is_satisfied().unwrap());

        // check the correctness of the result
        assert_eq!(res_var.value().unwrap(), res);
    }
    // assert!(false)
}

fn convert_boolean_var_to_fp_var(
    input: &[Boolean<Fq>],
) -> Result<FqVar, SynthesisError> {
    Boolean::le_bits_to_fp_var(input)
}
