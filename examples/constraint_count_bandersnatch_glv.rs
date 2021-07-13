use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::{BigInteger, FpParameters, PrimeField, UniformRand};
use ark_r1cs_std::{
    alloc::AllocVar,
    boolean::Boolean,
    eq::EqGadget,
    fields::{fp::FpVar, FieldVar},
    groups::{curves::twisted_edwards::AffineVar, CurveVar},
    uint128, R1CSVar,
};
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef,
    SynthesisError,
};
use bandersnatch::{
    constraints::FqVar, EdwardsAffine, EdwardsParameters, Fq, Fr, FrParameters,
    GLVParameters,
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
        let base_var = AffineVar::<EdwardsParameters, FpVar<Fq>>::new_witness(
            cs.clone(),
            || Ok(self.base),
        )
        .unwrap();
        #[cfg(debug_assertions)]
        println!("cs for base var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let phi_var = endomorphism_gadget(&base_var);

        // // sanity check
        // let phi = <EdwardsParameters as GLVParameters>::endomorphism(&self.base);
        // let phi_var_2 = AffineVar::new_witness(cs.clone(), ||Ok(phi)).unwrap();
        // phi_var.enforce_equal(&phi_var_2).unwrap();

        #[cfg(debug_assertions)]
        println!("cs for endomorphism var: {}", cs.num_constraints() - _cs_no);
        let _cs_no = cs.num_constraints();

        let (k1, k2) = scalar_decomposition_gadget(&self.scalar, cs.clone());

        #[cfg(debug_assertions)]
        println!(
            "cs for scalar decomposition var: {}",
            cs.num_constraints() - _cs_no
        );
        let _cs_no = cs.num_constraints();

        let res_var_recomputed = multi_scalar_mul_gadget(
            &base_var, &k1.0, &k1.1, &phi_var, &k2.0, &k2.1,
        );

        #[cfg(debug_assertions)]
        println!("cs for msm: {}", cs.num_constraints() - _cs_no);
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

/// Mapping a point G to phi(G):= lambda G where phi is the endomorphism
fn endomorphism_gadget(
    base_point: &AffineVar<EdwardsParameters, FqVar>,
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

/// Decompose a scalar s into k1, k2, s.t. s = k1 + lambda k2 mod r
/// for lambda = 8913659658109529928382530854484400854125314752504019737736543920008458395397
/// via a Babai's nearest plane algorithm.
fn scalar_decomposition_gadget(
    scalar: &Fr,
    cs: ConstraintSystemRef<Fq>,
) -> (
    (Vec<Boolean<Fq>>, Boolean<Fq>),
    (Vec<Boolean<Fq>>, Boolean<Fq>),
) {
    let (mut k1, mut k2) =
        <EdwardsParameters as GLVParameters>::scalar_decomposition(scalar);

    // todo: we need to prove that
    //  k1 + lambda k2 = s mod r
    // this seems to require a lot of constraints

    // println!("{}, {}, {}", scalar, k1, k2);

    let r_over_2: Fr =
        <FrParameters as FpParameters>::MODULUS_MINUS_ONE_DIV_TWO.into();

    let k1_is_pos = if k1 > r_over_2 {
        k1 = -k1;
        false
    } else {
        true
    };

    let k2_is_pos = if k2 > r_over_2 {
        k2 = -k2;
        false
    } else {
        true
    };

    // println!("{:?} {:?}",  k1.into_repr(), k2.into_repr());
    let k1 = k1.into_repr().to_bits_le();
    let k2 = k2.into_repr().to_bits_le();

    let len = std::cmp::max(get_bits(&k1), get_bits(&k2));
    // println!("{} {:?} {:?}", len, k1, k2);
    let mut k1_var = vec![];
    let mut k2_var = vec![];
    for i in 0..len as usize {
        k1_var.push(Boolean::new_witness(cs.clone(), || Ok(k1[i])).unwrap());
        k2_var.push(Boolean::new_witness(cs.clone(), || Ok(k2[i])).unwrap());
    }
    let k1_is_pos = Boolean::new_witness(cs.clone(), || Ok(k1_is_pos)).unwrap();
    let k2_is_pos = Boolean::new_witness(cs.clone(), || Ok(k2_is_pos)).unwrap();

    ((k1_var, k1_is_pos), (k2_var, k2_is_pos))
}

/// return the highest non-zero bits of a bit string.
fn get_bits(a: &[bool]) -> u16 {
    let mut res = 256;
    for e in a.iter().rev() {
        if !e {
            res -= 1;
        } else {
            return res;
        }
    }
    res
}

// Here we need to implement a customized MSM algorithm, since we know that
// the high bits of Fr are restricted to be small, i.e. ~ 128 bits.
// This MSM will save us some 128 doublings.
fn multi_scalar_mul_gadget(
    base: &AffineVar<EdwardsParameters, FqVar>,
    scalar_1: &[Boolean<Fq>],
    scalar_1_is_pos: &Boolean<Fq>,
    endor_base: &AffineVar<EdwardsParameters, FqVar>,
    scalar_2: &[Boolean<Fq>],
    scalar_2_is_pos: &Boolean<Fq>,
) -> AffineVar<EdwardsParameters, FqVar> {
    let zero = AffineVar::<EdwardsParameters, FqVar>::zero();

    let base = scalar_1_is_pos
        .select(base, &base.clone().negate().unwrap())
        .unwrap();

    let endor_base = scalar_2_is_pos
        .select(endor_base, &endor_base.clone().negate().unwrap())
        .unwrap();

    let sum = base.clone() + endor_base.clone();
    let len = scalar_1.len();
    let mut res = AffineVar::<EdwardsParameters, FqVar>::zero();
    for i in 0..len {
        res = res.double().unwrap();

        let add = scalar_1[len - i - 1]
            .select(
                // both bits are 1 => add the sum to self
                // first bit is 1 and second bit is 0 => add the base to self
                &scalar_2[len - i - 1].select(&sum, &base).unwrap(),
                // first bit is 0 and second bit is 1 => add the endor_base to self
                // bot bits are 0 => ignore
                &scalar_2[len - i - 1].select(&endor_base, &zero).unwrap(),
            )
            .unwrap();
        res += add;
    }

    res
}
