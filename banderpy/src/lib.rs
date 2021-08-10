#![allow(clippy::manual_strip)]

use ark_ec::AffineCurve;
use ark_ff::{BigInteger256, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use bandersnatch::{
    multi_scalar_mul, BandersnatchParameters, EdwardsAffine, Fr, GLVParameters,
};
use cpython::{py_fn, py_module_initializer, PyResult, Python};

py_module_initializer!(libbanderpy, |py, m| {
    m.add(py, "__doc__", "This module is implemented in Rust.")?;

    // ========================
    // generator
    // ========================
    m.add(py, "get_generator_rust", py_fn!(py, get_generator()))?;

    // ========================
    // random elements
    // ========================
    // generate a random point
    m.add(py, "random_point_rust", py_fn!(py, random_point()))?;

    // generate a random scalar
    m.add(py, "random_scalar_rust", py_fn!(py, random_scalar()))?;

    // ========================
    // arithmetics
    // ========================
    // point add
    m.add(py, "add_rust", py_fn!(py, add(a: Vec<u8>, b: Vec<u8>)))?;

    // point double
    m.add(py, "double_rust", py_fn!(py, double(a: Vec<u8>)))?;

    // group multiplication
    m.add(
        py,
        "mul_rust",
        py_fn!(py, mul(point: Vec<u8>, scalar: Vec<u8>)),
    )?;

    // glv multiplication
    m.add(
        py,
        "glv_rust",
        py_fn!(py, glv(point: Vec<u8>, scalar: Vec<u8>)),
    )?;

    // multi-scalar multiplication
    m.add(
        py,
        "msm_rust",
        py_fn!(py, msm(points: Vec<Vec<u8>>, scalars: Vec<Vec<u8>>)),
    )?;

    // ========================
    // serdes
    // ========================
    // serialize a point
    m.add(
        py,
        "point_to_string_rust",
        py_fn!(py, point_to_string(obj: Vec<u8>)),
    )?;

    // serialize a scalar
    m.add(
        py,
        "scalar_to_string_rust",
        py_fn!(py, scalar_to_string(obj: Vec<u8>)),
    )?;

    // serialize a point
    m.add(
        py,
        "point_serialize_rust",
        py_fn!(py, point_serialize(obj: Vec<u8>)),
    )?;

    // serialize a scalar
    m.add(
        py,
        "scalar_serialize_rust",
        py_fn!(py, scalar_serialize(obj: Vec<u8>)),
    )?;

    Ok(())
});

// ========================
// generator
// ========================
fn get_generator(_: Python) -> PyResult<Vec<u8>> {
    let mut buf = Vec::new();
    EdwardsAffine::prime_subgroup_generator()
        .serialize_unchecked(&mut buf)
        .expect("serialization error");

    Ok(buf)
}

// ========================
// random elements
// ========================
fn random_point(_: Python) -> PyResult<Vec<u8>> {
    let mut buf = Vec::new();
    let mut rng = rand::thread_rng();
    EdwardsAffine::rand(&mut rng)
        .serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}

fn random_scalar(_: Python) -> PyResult<Vec<u8>> {
    let mut buf = Vec::new();
    let mut rng = rand::thread_rng();
    Fr::rand(&mut rng)
        .serialize(&mut buf)
        .expect("serialization error");
    Ok(buf)
}
// ========================
// arithmetics
// ========================
fn add(_: Python, a: Vec<u8>, b: Vec<u8>) -> PyResult<Vec<u8>> {
    let a = EdwardsAffine::deserialize_unchecked(&a[..])
        .expect("deserialization error");
    let b = EdwardsAffine::deserialize_unchecked(&b[..])
        .expect("deserialization error");
    let c = a + b;
    let mut buf = Vec::new();
    c.serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}
fn double(_: Python, a: Vec<u8>) -> PyResult<Vec<u8>> {
    use ark_ec::group::Group;
    let a = EdwardsAffine::deserialize_unchecked(&a[..])
        .expect("deserialization error");

    let c = a.double();
    let mut buf = Vec::new();
    c.serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}

fn mul(_: Python, point: Vec<u8>, scalar: Vec<u8>) -> PyResult<Vec<u8>> {
    let p = EdwardsAffine::deserialize_unchecked(&point[..])
        .expect("deserialization error");
    let s = Fr::deserialize(&scalar[..]).expect("deserialization error");
    let r = p.mul(s);
    let mut buf = Vec::new();
    r.serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}

fn glv(_: Python, point: Vec<u8>, scalar: Vec<u8>) -> PyResult<Vec<u8>> {
    let p = EdwardsAffine::deserialize_unchecked(&point[..])
        .expect("deserialization error");
    let s = Fr::deserialize(&scalar[..]).expect("deserialization error");
    let r = <BandersnatchParameters as GLVParameters>::glv_mul(&p, &s);
    let mut buf = Vec::new();
    r.serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}

fn msm(
    _: Python,
    points: Vec<Vec<u8>>,
    scalar: Vec<Vec<u8>>,
) -> PyResult<Vec<u8>> {
    let points: Vec<EdwardsAffine> = points
        .iter()
        .map(|x| {
            EdwardsAffine::deserialize_unchecked(&x[..])
                .expect("deserialization error")
        })
        .collect();
    let scalars: Vec<BigInteger256> = scalar
        .iter()
        .map(|x| {
            Fr::deserialize(&x[..])
                .expect("deserialization error")
                .into_repr()
        })
        .collect();
    let res = multi_scalar_mul(&points, &scalars, 256);

    let mut buf = Vec::new();
    res.serialize_unchecked(&mut buf)
        .expect("serialization error");
    Ok(buf)
}

// ========================
// serdes
// ========================
fn point_to_string(_: Python, obj: Vec<u8>) -> PyResult<String> {
    Ok(EdwardsAffine::deserialize_unchecked(&obj[..])
        .expect("deserialization error")
        .to_string())
}

fn scalar_to_string(_: Python, obj: Vec<u8>) -> PyResult<String> {
    Ok(Fr::deserialize(&obj[..])
        .expect("deserialization error")
        .to_string())
}

fn point_serialize(_: Python, obj: Vec<u8>) -> PyResult<Vec<u8>> {
    let p = EdwardsAffine::deserialize_unchecked(&obj[..])
        .expect("deserialization error");
    let mut buf = Vec::new();
    p.serialize(&mut buf).expect("serialization error");

    Ok(buf)
}

fn scalar_serialize(_: Python, obj: Vec<u8>) -> PyResult<Vec<u8>> {
    let p = Fr::deserialize(&obj[..]).expect("deserialization error");
    let mut buf = Vec::new();
    p.serialize(&mut buf).expect("serialization error");

    Ok(buf)
}
