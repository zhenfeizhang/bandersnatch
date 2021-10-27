use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{PrimeField, UniformRand};
use ark_r1cs_std::{
    alloc::AllocVar, eq::EqGadget, fields::fp::FpVar,
    groups::curves::twisted_edwards::AffineVar,
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};
use bandersnatch::{
    constraints::{
        endomorphism_gadget, multi_scalar_mul_gadget,
        scalar_decomposition_gadget,
    },
    BandersnatchParameters, EdwardsAffine, Fq, Fr,
};

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
        let base_var =
            AffineVar::<BandersnatchParameters, FpVar<Fq>>::new_witness(
                cs.clone(),
                || Ok(self.base),
            )
            .unwrap();

        let scalar_var = FpVar::<Fq>::new_witness(cs.clone(), || {
            Ok(Fq::from(self.scalar.into_repr()))
        })?;

        #[cfg(debug_assertions)]
        println!("cs for setup: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let phi_var = endomorphism_gadget(cs.clone(), &base_var)?;

        // // sanity check
        // let phi = <EdwardsParameters as GLVParameters>::endomorphism(&self.base);
        // let phi_var_2 = AffineVar::new_witness(cs.clone(), ||Ok(phi)).unwrap();
        // phi_var.enforce_equal(&phi_var_2).unwrap();

        #[cfg(debug_assertions)]
        println!("cs for endomorphism var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let (k, k2_sign) =
            scalar_decomposition_gadget(cs.clone(), &scalar_var)?;

        #[cfg(debug_assertions)]
        println!(
            "cs for scalar decomposition var: {}",
            cs.num_constraints() - _cs_no
        );
        let _cs_no = cs.num_constraints();

        let res_var_recomputed = multi_scalar_mul_gadget(
            &base_var, &k[0], &phi_var, &k[1], &k2_sign,
        )?;

        #[cfg(debug_assertions)]
        println!("cs for msm: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let res_var =
            AffineVar::<BandersnatchParameters, FpVar<Fq>>::new_witness(
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

        #[cfg(debug_assertions)]
        println!("total constraints : {}", cs.num_constraints());

        Ok(())
    }
}
