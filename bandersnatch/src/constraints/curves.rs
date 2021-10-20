use crate::{constraints::FqVar, *};
use ark_r1cs_std::groups::curves::{
    short_weierstrass::ProjectiveVar, twisted_edwards::AffineVar,
};

/// A variable that is the R1CS equivalent of `crate::EdwardsAffine`.
pub type EdwardsAffineVar = AffineVar<BandersnatchParameters, FqVar>;

/// A variable that is the R1CS equivalent of `crate::EdwardsProjective`
pub type SWProjectiveVar = ProjectiveVar<BandersnatchParameters, FqVar>;

#[test]
fn test() {
    ark_curve_constraint_tests::curves::te_test::<_, EdwardsAffineVar>()
        .unwrap();
    ark_curve_constraint_tests::curves::sw_test::<
        BandersnatchParameters,
        SWProjectiveVar,
    >()
    .unwrap();
    ark_curve_constraint_tests::curves::group_test::<
        EdwardsProjective,
        Fq,
        EdwardsAffineVar,
    >()
    .unwrap();
}
