use crate::{
    BandersnatchParameters, EdwardsAffine, EdwardsProjective, Fq, Fr,
    FrParameters,
};
use ark_ec::{AffineCurve, ModelParameters, ProjectiveCurve};
use ark_ff::{
    field_new, prelude::*, BigInteger, BigInteger256, FpParameters, One,
};
use ark_std::{cmp::max, vec, vec::Vec, Zero};
use num_bigint::BigUint;

/// The GLV parameters that are useful to compute the endomorphism
/// and scalar decomposition.
pub trait GLVParameters: Send + Sync + 'static + ModelParameters {
    type CurveAffine;
    type CurveProjective;

    // phi(P) = lambda*P for all P
    // constants that are used to calculate phi(P)
    const COEFF_A1: Self::BaseField;
    const COEFF_A2: Self::BaseField;
    const COEFF_A3: Self::BaseField;
    const COEFF_B1: Self::BaseField;
    const COEFF_B2: Self::BaseField;
    const COEFF_B3: Self::BaseField;
    const COEFF_C1: Self::BaseField;
    const COEFF_C2: Self::BaseField;

    // constants that are used to perform scalar decomposition
    // This is a matrix which is practically the LLL reduced bases
    const COEFF_N11: Self::ScalarField;
    const COEFF_N12: Self::ScalarField;
    const COEFF_N21: Self::ScalarField;
    const COEFF_N22: Self::ScalarField;

    /// mapping a point G to phi(G):= lambda G where psi is the endomorphism
    fn endomorphism(base: &Self::CurveAffine) -> Self::CurveAffine;

    /// decompose a scalar s into k1, k2, s.t. s = k1 + lambda k2
    fn scalar_decomposition(
        k: &Self::ScalarField,
    ) -> (Self::ScalarField, Self::ScalarField);

    /// perform GLV multiplication
    fn glv_mul(
        base: &Self::CurveAffine,
        scalar: &Self::ScalarField,
    ) -> Self::CurveProjective;
}

impl GLVParameters for BandersnatchParameters {
    type CurveAffine = crate::EdwardsAffine;
    type CurveProjective = crate::EdwardsProjective;

    // phi(P) = lambda*P for all P
    // constants that are used to calculate phi(P)
    const COEFF_A1: Self::BaseField = field_new!(
        Fq,
        "16179988757916560824577558193084210236647645729299773892093730683504906651604"
    );

    const COEFF_A2: Self::BaseField = field_new!(
        Fq,
        "37446463827641770816307242315180085052603635617490163568005256780843403514036"
    );

    const COEFF_A3: Self::BaseField = field_new!(
        Fq,
        "14989411347484419663140498193005880785086916883037474254598401919095177670477"
    );

    const COEFF_B1: Self::BaseField = field_new!(
        Fq,
        "37446463827641770816307242315180085052603635617490163568005256780843403514036"
    );

    const COEFF_B2: Self::BaseField = field_new!(
        Fq,
        "36553259151239542273674161596529768046449890757310263666255995151154432137034"
    );

    const COEFF_B3: Self::BaseField = field_new!(
        Fq,
        "15882616023886648205773578911656197791240661743217374156347663548784149047479"
    );

    const COEFF_C1: Self::BaseField = field_new!(
        Fq,
        "42910309089382041158038545419309140955400939872179826051492616687477682993077"
    );

    const COEFF_C2: Self::BaseField = field_new!(
        Fq,
        "9525566085744149321409195088876824882289612628347811771111042012460898191436"
    );

    // constants that are used to perform scalar decomposition
    // This is a matrix which is practically the LLL reduced bases
    // N = Matrix(
    // [[113482231691339203864511368254957623327,
    // 10741319382058138887739339959866629956],
    // [21482638764116277775478679919733259912,
    // -113482231691339203864511368254957623327]])

    const COEFF_N11: Self::ScalarField =
        field_new!(Fr, "113482231691339203864511368254957623327");

    const COEFF_N12: Self::ScalarField =
        field_new!(Fr, "10741319382058138887739339959866629956");

    const COEFF_N21: Self::ScalarField =
        field_new!(Fr, "21482638764116277775478679919733259912");

    const COEFF_N22: Self::ScalarField =
        field_new!(Fr, "-113482231691339203864511368254957623327");

    /// Mapping a point G to phi(G):= lambda G where phi is the endomorphism
    fn endomorphism(base: &Self::CurveAffine) -> Self::CurveAffine {
        let mut x = base.x;
        let mut y = base.y;
        let mut z = y;

        // z = y;
        let fy = Self::COEFF_A1 * (y + Self::COEFF_A2) * (y + Self::COEFF_A3);
        let gy = Self::COEFF_B1 * (y + Self::COEFF_B2) * (y + Self::COEFF_B3);
        let hy = (y + Self::COEFF_C1) * (y + Self::COEFF_C2);

        x = x * fy * hy;
        y = gy * z;
        z = hy * z;

        Self::CurveProjective::new(x, y, Fq::one(), z).into_affine()
    }

    /// Decompose a scalar s into k1, k2, s.t. s = k1 + lambda k2
    /// via a Babai's nearest plane algorithm.
    fn scalar_decomposition(
        scalar: &Self::ScalarField,
    ) -> (Self::ScalarField, Self::ScalarField) {
        let tmp: BigInteger256 = (*scalar).into();
        let scalar_z: BigUint = tmp.into();

        let tmp: BigInteger256 = Self::COEFF_N11.into();
        let n11: BigUint = tmp.into();

        let tmp: BigInteger256 = Self::COEFF_N12.into();
        let n12: BigUint = tmp.into();

        let r: BigUint = <FrParameters as FpParameters>::MODULUS.into();

        // beta = vector([n,0]) * self.curve.N_inv
        let beta_1 = scalar_z.clone() * n11;
        let beta_2 = scalar_z * n12;

        let beta_1 = beta_1 / r.clone();
        let beta_2 = beta_2 / r;

        // b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
        let beta_1 = Fr::from(beta_1);
        let beta_2 = Fr::from(beta_2);
        let b1 = beta_1 * Self::COEFF_N11 + beta_2 * Self::COEFF_N21;
        let b2 = beta_1 * Self::COEFF_N12 + beta_2 * Self::COEFF_N22;

        let k1 = (*scalar) - b1;
        let k2 = -b2;
        (k1, k2)
    }

    /// perform GLV multiplication
    fn glv_mul(
        base: &Self::CurveAffine,
        scalar: &Self::ScalarField,
    ) -> Self::CurveProjective {
        let psi_base = Self::endomorphism(&base);
        let (k1, k2) = Self::scalar_decomposition(scalar);
        two_scalar_mul(&base, &k1, &psi_base, &k2)
    }
}

// Here we need to implement a customized MSM algorithm, since we know that
// the high bits of Fr are restricted to be small, i.e. ~ 128 bits.
// This MSM will save us some 128 doublings.
pub fn two_scalar_mul(
    base: &crate::EdwardsAffine,
    scalar_1: &Fr,
    endor_base: &crate::EdwardsAffine,
    scalar_2: &Fr,
) -> crate::EdwardsProjective {
    let mut b1 = (*base).into_projective();
    let mut s1 = *scalar_1;
    let mut b2 = (*endor_base).into_projective();
    let mut s2 = *scalar_2;

    let r_over_2: Fr =
        <FrParameters as FpParameters>::MODULUS_MINUS_ONE_DIV_TWO.into();

    if s1 > r_over_2 {
        b1 = -b1;
        s1 = -s1;
    }
    if s2 > r_over_2 {
        b2 = -b2;
        s2 = -s2;
    }
    let s1: BigInteger256 = s1.into();
    let s2: BigInteger256 = s2.into();

    let b1b2 = b1 + b2;

    let s1_bits = s1.to_bits_le();
    let s2_bits = s2.to_bits_le();
    let s1_len = get_bits(&s1_bits);
    let s2_len = get_bits(&s2_bits);
    let len = max(s1_len, s2_len) as usize;

    let mut res = crate::EdwardsProjective::zero();
    for i in 0..len {
        res = res.double();
        if s1_bits[len - i - 1] && !s2_bits[len - i - 1] {
            res += b1
        }
        if !s1_bits[len - i - 1] && s2_bits[len - i - 1] {
            res += b2
        }
        if s1_bits[len - i - 1] && s2_bits[len - i - 1] {
            res += b1b2
        }
    }
    res
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

/// The result of this function is only approximately `ln(a)`
/// [`Explanation of usage`]
///
/// [`Explanation of usage`]: https://github.com/scipr-lab/zexe/issues/79#issue-556220473
fn ln_without_floats(a: usize) -> usize {
    // log2(a) * ln(2)
    (ark_std::log2(a) * 69 / 100) as usize
}

pub fn multi_scalar_mul_with_glv(
    bases: &[EdwardsAffine],
    scalars: &[Fr],
) -> EdwardsProjective {
    let mut bases_ext = Vec::new();
    let mut scalars_ext = Vec::new();

    let r_over_2: Fr =
        <FrParameters as FpParameters>::MODULUS_MINUS_ONE_DIV_TWO.into();

    for i in 0..bases.len() {
        let phi =
            <BandersnatchParameters as GLVParameters>::endomorphism(&bases[i]);
        let (k1, k2) =
            <BandersnatchParameters as GLVParameters>::scalar_decomposition(
                &scalars[i],
            );
        if k1 > r_over_2 {
            bases_ext.push(-bases[i]);
            scalars_ext.push((-k1).into_repr());
        } else {
            bases_ext.push(bases[i]);
            scalars_ext.push(k1.into_repr());
        }

        if k2 > r_over_2 {
            bases_ext.push(-phi);
            scalars_ext.push((-k2).into_repr());
        } else {
            bases_ext.push(phi);
            scalars_ext.push(k2.into_repr());
        }
    }
    let mut len = 0;
    for e in scalars_ext.iter() {
        let l = get_bits(&e.to_bits_le());
        // println!("{}", l);
        if l > len {
            len = l
        }
    }
    // println!("{} {}", bases.len(), len);
    multi_scalar_mul(&bases_ext, &scalars_ext, len as u32)
}

pub fn multi_scalar_mul<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
    len: u32,
) -> G::Projective {
    let size = ark_std::cmp::min(bases.len(), scalars.len());
    let scalars = &scalars[..size];
    let bases = &bases[..size];
    let scalars_and_bases_iter =
        scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero());

    let c = if size < 32 {
        3
    } else {
        ln_without_floats(size) + 2
    };

    let num_bits = len;
    // <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
    let fr_one = G::ScalarField::one().into_repr();

    let zero = G::Projective::zero();
    let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

    // Each window is of size `c`.
    // We divide up the bits 0..num_bits into windows of size `c`, and
    // in parallel process each such window.
    let window_sums: Vec<_> = ark_std::cfg_into_iter!(window_starts)
        .map(|w_start| {
            let mut res = zero;
            // We don't need the "zero" bucket, so we only have 2^c - 1 buckets.
            let mut buckets = vec![zero; (1 << c) - 1];
            // This clone is cheap, because the iterator contains just a
            // pointer and an index into the original vectors.
            scalars_and_bases_iter.clone().for_each(|(&scalar, base)| {
                if scalar == fr_one {
                    // We only process unit scalars once in the first window.
                    if w_start == 0 {
                        res.add_assign_mixed(base);
                    }
                } else {
                    let mut scalar = scalar;

                    // We right-shift by w_start, thus getting rid of the
                    // lower bits.
                    scalar.divn(w_start as u32);

                    // We mod the remaining bits by 2^{window size}, thus taking `c` bits.
                    let scalar = scalar.as_ref()[0] % (1 << c);

                    // If the scalar is non-zero, we update the corresponding
                    // bucket.
                    // (Recall that `buckets` doesn't have a zero bucket.)
                    if scalar != 0 {
                        buckets[(scalar - 1) as usize].add_assign_mixed(base);
                    }
                }
            });

            // Compute sum_{i in 0..num_buckets} (sum_{j in i..num_buckets} bucket[j])
            // This is computed below for b buckets, using 2b curve additions.
            //
            // We could first normalize `buckets` and then use mixed-addition
            // here, but that's slower for the kinds of groups we care about
            // (Short Weierstrass curves and Twisted Edwards curves).
            // In the case of Short Weierstrass curves,
            // mixed addition saves ~4 field multiplications per addition.
            // However normalization (with the inversion batched) takes ~6
            // field multiplications per element,
            // hence batch normalization is a slowdown.

            // `running_sum` = sum_{j in i..num_buckets} bucket[j],
            // where we iterate backward from i = num_buckets to 0.
            let mut running_sum = G::Projective::zero();
            buckets.into_iter().rev().for_each(|b| {
                running_sum += &b;
                res += &running_sum;
            });
            res
        })
        .collect();

    // We store the sum for the lowest window.
    let lowest = *window_sums.first().unwrap();

    // We're traversing windows from high to low.
    lowest
        + &window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, sum_i| {
                total += sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            })
}
