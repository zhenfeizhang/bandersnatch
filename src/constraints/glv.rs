use ark_ec::{
    short_weierstrass_jacobian::GroupProjective, ProjectiveCurve,
    SWModelParameters,
};
use ark_ff::Field;
use ark_r1cs_std::{
    groups::curves::short_weierstrass::AffineVar,
    prelude::{
        AllocVar, CondSelectGadget, CurveVar, EqGadget, FieldOpsBounds,
        FieldVar, GroupOpsBounds,
    },
    R1CSVar, ToBitsGadget, ToBytesGadget,
};
use ark_std::{
    fmt::Debug,
    ops::{AddAssign, SubAssign},
};
use crate::{constraints::FqVar, GLVParameters, EdwardsParameters};


/// This GLVVar trait inherit the CurveVar trait
pub trait GLVVar<C: ProjectiveCurve, ConstraintF: Field>:
    'static
    + Sized
    + Clone
    + Debug
    + R1CSVar<ConstraintF, Value = C>
    + ToBitsGadget<ConstraintF>
    + ToBytesGadget<ConstraintF>
    + EqGadget<ConstraintF>
    + CondSelectGadget<ConstraintF>
    + AllocVar<C, ConstraintF>
    + AllocVar<C::Affine, ConstraintF>
    + for<'a> GroupOpsBounds<'a, C, Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + AddAssign<C>
    + SubAssign<C>
    + AddAssign<Self>
    + SubAssign<Self>
    + CurveVar<C, ConstraintF>
{
    /// Mapping a point G to phi(G):= lambda G where phi is the endomorphism
    fn endomorphism(&self) -> Self;

    /// Decompose a scalar s into k1, k2, s.t. s = k1 + lambda k2
    /// via a Babai's nearest plane algorithm.
    fn scalar_decomposition(
        scalar: C::ScalarField,
    ) -> (C::ScalarField, C::ScalarField);

    /// perform GLV multiplication
    fn glv_mul(&self, scalar: &C::ScalarField) -> Self;
}

impl<P, F> GLVVar<GroupProjective<P>, <P::BaseField as Field>::BasePrimeField>
    for AffineVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn endomorphism(&self) -> Self {
        let coeff_a1_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A1);
        let coeff_a2_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A2);
        let coeff_a3_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_A3);
        let coeff_b1_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B1);
        let coeff_b2_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B2);
        let coeff_b3_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_B3);
        let coeff_c1_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_C1);
        let coeff_c1_var = FqVar::constant(<EdwardsParameters as GLVParameters>::COEFF_C2);

        let x_var: FqVar = self.x.clone();
        let y_var = self.y.clone().into();
        let z_var = y_var;

        let fy = coeff_a1_var * (y_var + coeff_a2_var) * (y_var + coeff_a3_var);


        todo!()
    }

    fn scalar_decomposition(
        _scalar: P::ScalarField,
    ) -> (P::ScalarField, P::ScalarField) {
        todo!()
    }

    /// perform GLV multiplication
    fn glv_mul(&self, _scalar: &P::ScalarField) -> Self {
        todo!()
    }
}

#[test]
fn test_glv_var() {}