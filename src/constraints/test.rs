// use ark_ec::AffineCurve;
// // use ark_ec::ModelParameters;
// // use ark_r1cs_std::uint8::UInt8;
// use crate::*;
// use ark_ec::ProjectiveCurve;
// use ark_ff::BigInteger;
// use ark_ff::UniformRand;
// use ark_r1cs_std::boolean::Boolean;
// use ark_r1cs_std::groups::CurveVar;
// use ark_r1cs_std::{
//     alloc::AllocVar,
//     // eq::EqGadget,
//     fields::fp::FpVar,
//     groups::curves::twisted_edwards::AffineVar,
//     // uint8::UInt8,
// };
// use ark_relations::r1cs::ConstraintSynthesizer;
// use ark_relations::r1cs::ConstraintSystemRef;
// use ark_relations::r1cs::{ConstraintSystem, SynthesisError};
// // use std::vec::Vec;
// use ark_r1cs_std::eq::EqGadget;

// #[test]
// fn test_constraints() {
//     // in this example we are going to argue about the statement
//     //      H = xG
//     // for some private H, x and G, where
//     //  - G is a bandersnatch group element in the Affine form,
//     //  - x is a scalar field element
//     //  - H is a bandersnatch group element in the Project form,

//     let mut rng = ark_std::test_rng();

//     // set up the proper inputs:
//     let point_g = EdwardsAffine::rand(&mut rng);
//     let x = Fr::rand(&mut rng);
//     let point_h = point_g.mul(x).into_affine();
//     let circuit = GroupOpCircuit {
//         base: point_g,
//         scalar: x,
//         res: point_h,
//     };

//     println!("{:?}", point_g);
//     println!("{:?}", x);
//     println!("{:?}", point_h);

//     let sanity_cs = ConstraintSystem::<Fq>::new_ref();

//     circuit.generate_constraints(sanity_cs.clone()).unwrap();
//     assert!(sanity_cs.is_satisfied().unwrap());

//     // let open = Randomness::<JubJub>(Fr::rand(&mut rng));
//     // let commit = pedersen_commit(input.as_ref(), &param, &open);

//     // // build the input of the circuit
//     // let circuit = PedersenComCircuit {
//     // 	param: param.clone(),
//     // 	input: input.to_vec(),
//     // 	open,
//     // 	commit,
//     // };

//     // // check the circuit is satisfied
//     // assert!(sanity_check(circuit.clone()));

//     // // generate the CRS for Groth16
//     // let zk_param = groth_param_gen(param);

//     // // generate the proof
//     // let proof = groth_proof_gen(&zk_param, circuit, &[0u8; 32]);

//     // // verify the proof
//     // assert!(groth_verify(&zk_param, &proof, &commit));

//     assert!(false)
// }

// /// a circuit for the relation:
// ///   res = scalar * base
// struct GroupOpCircuit {
//     base: EdwardsAffine,
//     scalar: Fr,
//     res: EdwardsAffine,
// }

// impl ConstraintSynthesizer<Fq> for GroupOpCircuit {
//     fn generate_constraints(
//         self,
//         cs: ConstraintSystemRef<Fq>,
//     ) -> Result<(), SynthesisError> {
//         let base_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
//             cs.clone(),
//             || Ok(self.base),
//         )
//         .unwrap();

//         let mut scalar_bits_var = vec![];
//         for e in self.scalar.0.to_bits_le() {
//             scalar_bits_var.push(Boolean::new_witness(cs.clone(), || Ok(e))?)
//         }
//         println!("{:?}", scalar_bits_var);
//         let res_var_recomp = base_var.scalar_mul_le(scalar_bits_var.iter())?;

//         let res_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
//             cs.clone(),
//             || Ok(self.res),
//         )
//         .unwrap();

//         res_var.enforce_equal(&res_var_recomp)
//     }
// }
