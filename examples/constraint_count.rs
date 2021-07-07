use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ed_on_bls12_381::{EdwardsAffine, EdwardsParameters, Fq, Fr};
use ark_ff::{BigInteger, UniformRand};
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

fn main() {
    // in this example we are going to argue about the statement
    //      H = xG
    // for some private H, x and G, where
    //  - G is a bandersnatch group element in the Affine form,
    //  - x is a scalar field element
    //  - H is a bandersnatch group element in the Project form,

    let mut rng = ark_std::test_rng();

    let point_g = EdwardsAffine::rand(&mut rng);
    let x = Fr::rand(&mut rng);
    let point_h = point_g.mul(x).into_affine();
    let circuit = GroupOpCircuit {
        base: point_g,
        scalar: x,
        res: point_h,
    };

    // println!("{:?}", point_g);
    // println!("{:?}", x);
    // println!("{:?}", point_h);

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

        let mut scalar_bits_var = vec![];
        for e in self.scalar.0.to_bits_le() {
            scalar_bits_var.push(Boolean::new_witness(cs.clone(), || Ok(e))?)
        }

        #[cfg(debug_assertions)]
        println!("cs for scalar var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let res_var_recomputed =
            base_var.scalar_mul_le(scalar_bits_var.iter().rev())?;

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
