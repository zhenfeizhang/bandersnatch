#[macro_use]
extern crate criterion;

use ark_ec::{
    msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::{fields::PrimeField, BigInteger256};
use ark_std::{
    ops::MulAssign,
    rand::{RngCore, SeedableRng},
    UniformRand,
};
use bandersnatch::{
    BandersnatchParameters, EdwardsAffine as BandersnatchAffine, FixedBaseGLV,
    FixedBaseMul, Fr, GLVParameters,
};
use criterion::Criterion;
use rand_chacha::ChaCha20Rng;

criterion_group!(
    bandersnatch_bench,
    bench_endomorphism,
    bench_decomposition,
    bench_2sm,
    bench_bandersnatch,
    bench_jubjub,
    bench_ed_on_bls_12_377,
    bench_msm,
    bench_bls12_381_g1,
    bench_bls12_381_g2,
    bench_bls12_381_pairing,
    bench_bls12_377_g1,
    bench_bls12_377_g2,
    bench_bls12_377_pairing,
    bench_blst_g1,
    bench_blst_g2,
);

criterion_main!(bandersnatch_bench);

fn bench_decomposition(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("micro benchmark");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("decomposition");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        b.iter(|| BandersnatchParameters::scalar_decomposition(&r))
    });
    bench_group.finish();
}

fn bench_endomorphism(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("micro benchmark");

    let base_point = bandersnatch::EdwardsAffine::prime_subgroup_generator();

    let bench_str = format!("endomorphism");
    bench_group.bench_function(bench_str, move |b| {
        b.iter(|| BandersnatchParameters::endomorphism(&base_point))
    });
    bench_group.finish();
}

fn bench_2sm(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("micro benchmark");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let base_point = bandersnatch::EdwardsAffine::prime_subgroup_generator();
    let psi_point = BandersnatchParameters::endomorphism(&base_point);

    let bench_str = format!("two-scalar-mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        let (r1, r2) = BandersnatchParameters::scalar_decomposition(&r);

        b.iter(|| {
            bandersnatch::two_scalar_mul(&base_point, &r1, &psi_point, &r2)
        })
    });
    bench_group.finish();
}

fn bench_msm(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("micro benchmark");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    for d in [
        2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
    ]
    .iter()
    {
        let base: Vec<BandersnatchAffine> = (0..*d)
            .map(|_| BandersnatchAffine::rand(&mut rng))
            .collect();
        let base_clone = base.clone();

        let scalar: Vec<Fr> = (0..*d).map(|_| Fr::rand(&mut rng)).collect();
        let scalar_repr: Vec<BigInteger256> =
            scalar.iter().map(|x| (*x).into_repr()).collect();

        let bench_str = format!("ark-multi-scalar-mul for dim {}", d);
        bench_group.bench_function(bench_str, move |b| {
            b.iter(|| {
                let _ = VariableBaseMSM::multi_scalar_mul(
                    &base_clone,
                    &scalar_repr,
                );
            });
        });

        let bench_str = format!("glv-multi-scalar-mul for dim {}", d);
        bench_group.bench_function(bench_str, move |b| {
            b.iter(|| {
                let _ = bandersnatch::multi_scalar_mul_with_glv(&base, &scalar);
            });
        });
    }
    bench_group.finish();
}

fn bench_bandersnatch(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("bandersnatch curve");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let mut base_point =
        bandersnatch::EdwardsAffine::prime_subgroup_generator();

    let base_vector: Vec<bandersnatch::EdwardsProjective> =
        <BandersnatchParameters as FixedBaseMul>::preprocessing(
            base_point.into_projective(),
        );

    let base_ext_vector: Vec<[bandersnatch::EdwardsProjective; 4]> =
        <BandersnatchParameters as FixedBaseGLV>::preprocessing(
            base_point.into_projective(),
        );

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("variable base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        b.iter(|| {
            base_point.mul_assign(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        b.iter(|| {
            <BandersnatchParameters as FixedBaseMul>::fixed_base_mul(
                &base_vector,
                r,
            )
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("glv mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = BandersnatchParameters::glv_mul(&base_point, &r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("fix base mul with GLV");
    bench_group.bench_function(bench_str, move |b| {
        let r = bandersnatch::Fr::rand(&mut rng);
        b.iter(|| {
            <BandersnatchParameters as FixedBaseGLV>::fixed_base_mul(
                &base_ext_vector,
                r,
            )
        })
    });

    bench_group.finish();
}

fn bench_jubjub(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("JubJub curve");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let mut base_point =
        ark_ed_on_bls12_381::EdwardsAffine::prime_subgroup_generator();
    let mut random_point =
        ark_ed_on_bls12_381::EdwardsAffine::from_random_bytes(bytes.as_ref())
            .unwrap();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_ed_on_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            random_point.mul_assign(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_ed_on_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            base_point.mul_assign(r);
        })
    });
    bench_group.finish();
}

fn bench_ed_on_bls_12_377(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("ed_on_bls_12_377 curve");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let mut base_point =
        ark_ed_on_bls12_377::EdwardsAffine::prime_subgroup_generator();
    let mut random_point =
        ark_ed_on_bls12_377::EdwardsAffine::from_random_bytes(bytes.as_ref())
            .unwrap();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_ed_on_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            random_point.mul_assign(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_ed_on_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            base_point.mul_assign(r);
        })
    });
    bench_group.finish();
}

fn bench_bls12_381_g1(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-381 curve G1");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let base_point = ark_bls12_381::G1Affine::prime_subgroup_generator();

    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = base_point.mul(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let r = ark_bls12_381::Fr::rand(&mut rng);
    let random_point = base_point.mul(r).into_affine();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = random_point.mul(r);
        })
    });

    bench_group.finish();
}

fn bench_bls12_381_g2(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-381 curve G2");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let base_point = ark_bls12_381::G2Affine::prime_subgroup_generator();

    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = base_point.mul(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let r = ark_bls12_381::Fr::rand(&mut rng);
    let random_point = base_point.mul(r).into_affine();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_381::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = random_point.mul(r);
        })
    });

    bench_group.finish();
}

fn bench_bls12_381_pairing(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-381 curve pairing");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    let mut g1 =
        ark_bls12_381::G1Affine::prime_subgroup_generator().into_projective();
    let mut g2 =
        ark_bls12_381::G2Affine::prime_subgroup_generator().into_projective();

    let bench_str = format!("random base pairing");
    bench_group.bench_function(bench_str, move |b| {
        let r1 = ark_bls12_381::Fr::rand(&mut rng);
        let r2 = ark_bls12_381::Fr::rand(&mut rng);

        g1.mul_assign(r1);
        g2.mul_assign(r2);

        b.iter(|| {
            let _ = ark_bls12_381::Bls12_381::pairing(g1, g2);
        })
    });

    bench_group.finish();
}

fn bench_bls12_377_g1(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-377 curve G1");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let base_point = ark_bls12_377::G1Affine::prime_subgroup_generator();

    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = base_point.mul(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let r = ark_bls12_377::Fr::rand(&mut rng);
    let random_point = base_point.mul(r).into_affine();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = random_point.mul(r);
        })
    });

    bench_group.finish();
}

fn bench_bls12_377_g2(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-377 curve G2");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut bytes = [0u8; 32];
    rng.fill_bytes(&mut bytes);

    let base_point = ark_bls12_377::G2Affine::prime_subgroup_generator();

    let bench_str = format!("fix base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = base_point.mul(r);
        })
    });

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let r = ark_bls12_377::Fr::rand(&mut rng);
    let random_point = base_point.mul(r).into_affine();

    let bench_str = format!("random base mul");
    bench_group.bench_function(bench_str, move |b| {
        let r = ark_bls12_377::Fr::rand(&mut rng);
        b.iter(|| {
            let _ = random_point.mul(r);
        })
    });

    bench_group.finish();
}

fn bench_bls12_377_pairing(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("BLS12-377 curve pairing");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    let mut g1 =
        ark_bls12_377::G1Affine::prime_subgroup_generator().into_projective();
    let mut g2 =
        ark_bls12_377::G2Affine::prime_subgroup_generator().into_projective();

    let bench_str = format!("random base pairing");
    bench_group.bench_function(bench_str, move |b| {
        let r1 = ark_bls12_377::Fr::rand(&mut rng);
        let r2 = ark_bls12_377::Fr::rand(&mut rng);

        g1.mul_assign(r1);
        g2.mul_assign(r2);

        b.iter(|| {
            let _ = ark_bls12_377::Bls12_377::pairing(g1, g2);
        })
    });

    bench_group.finish();
}

fn bench_blst_g1(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("blst/bls12-381 g1");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
    let mut ikm = [0u8; 32];
    rng.fill_bytes(&mut ikm);

    let sk = blst::min_pk::SecretKey::key_gen(&ikm, &[]).unwrap();

    let bench_str = format!("fixed base mul");
    bench_group.bench_function(bench_str, move |b| {
        b.iter(|| {
            let _ = sk.sk_to_pk();
        })
    });

    // the following code is cleaner but there is some memory overflow
    // let pk = sk.sk_to_pk();
    // let mut g = unsafe { blst::blst_p1_generator() };
    // let mut res = unsafe { blst::blst_p1_generator() as *mut blst_p1 };
    // // let mut f = g;
    // // let mut res: *mut blst_p1 = h;
    // // let mut f = unsafe { *(::std::ptr::null::<blst_p1>()) };
    // // let mut res = unsafe { *(::std::ptr::null::<blst_p1>()) };

    // let r = ark_bls12_381::Fr::rand(&mut rng);
    // let mut r_bytes = vec![0u8; 32];
    // // r.serialize(&mut r_bytes).unwrap();
    // let r_pt = r_bytes.as_ptr();

    // unsafe { blst::blst_p1_mult(res, g, r_pt, 0) }

    // let bench_str = format!("fixed base mul");
    // bench_group.bench_function(bench_str, move |b| {
    //     // b.iter(|| unsafe { blst::blst_p1_mult(&mut res, &f, r_pt, 256) })
    // });

    bench_group.finish();
}

fn bench_blst_g2(c: &mut Criterion) {
    let mut bench_group = c.benchmark_group("blst/bls12-381 g2");

    let mut rng = ChaCha20Rng::from_seed([0u8; 32]);

    let mut ikm = [0u8; 32];
    rng.fill_bytes(&mut ikm);

    let sk = blst::min_sig::SecretKey::key_gen(&ikm, &[]).unwrap();

    let bench_str = format!("fixed base mul");
    bench_group.bench_function(bench_str, move |b| {
        b.iter(|| {
            let _ = sk.sk_to_pk();
        })
    });

    bench_group.finish();
}
