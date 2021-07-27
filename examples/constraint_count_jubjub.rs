use ark_ec::ProjectiveCurve;
use ark_ed_on_bls12_381::{Fq, Fr, JubjubParameters, SWProjective};
use ark_ff::{BigInteger, PrimeField, UniformRand};
use ark_r1cs_std::{
    alloc::AllocVar,
    boolean::Boolean,
    eq::EqGadget,
    fields::fp::FpVar,
    groups::{curves::short_weierstrass::ProjectiveVar, CurveVar},
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};

fn main() {
    // in this example we are going to argue about the statement
    //      H = xG
    // for some private H, x and G, where
    //  - G is a Jubjub group element in the Affine form,
    //  - x is a scalar field element
    //  - H is a Jubjub group element in the Affine form,

    let mut rng = ark_std::test_rng();

    let point_g = SWProjective::rand(&mut rng);
    let x = Fr::rand(&mut rng);
    let point_h = point_g.mul(x.into_repr());
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
    base: SWProjective,
    scalar: Fr,
    res: SWProjective,
}

impl ConstraintSynthesizer<Fq> for GroupOpCircuit {
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<Fq>,
    ) -> Result<(), SynthesisError> {
        let _cs_no = cs.num_constraints();
        let base_var =
            ProjectiveVar::<JubjubParameters, FpVar<Fq>>::new_witness(
                cs.clone(),
                || Ok(self.base),
            )
            .unwrap();

        #[cfg(debug_assertions)]
        println!("cs for base var: {}", cs.num_constraints() - _cs_no);
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

        let res_var =
            ProjectiveVar::<JubjubParameters, FpVar<Fq>>::new_witness(
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
