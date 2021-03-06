bandersnatch-rust
------
![GitHub branch checks state](https://img.shields.io/github/checks-status/zhenfeizhang/bandersnatch/main)
![docs.rs](https://img.shields.io/docsrs/bandersnatch/0.1.1)
![Crates.io (version)](https://img.shields.io/crates/dv/bandersnatch/0.1.1)
![GitHub](https://img.shields.io/github/license/zhenfeizhang/bandersnatch)

This is a reference implementation of [Bandersnatch curve](https://ethresear.ch/t/introducing-bandersnatch-a-fast-elliptic-curve-built-over-the-bls12-381-scalar-field/9957) using [Arkwork](https://github.com/arkworks-rs/curves)'s framework in Rust.
The spec of the curve is available 
[here](https://github.com/asanso/Bandersnatch/blob/main/README.md).
There is also a Python reference implementation [here](https://github.com/asanso/Bandersnatch/),
and a python wrapper of this library [banderpy](https://github.com/zhenfeizhang/bandersnatch/tree/main/banderpy).

# Logistics

- This code is released under MIT license.
- This code is not audited and may contain severe security flaws. Use at your own risk.
- Version 0.1.1.
- This repo is upstreamed to Arkworks [curve](https://github.com/arkworks-rs/curves/) crate.

# Change log

__0.1.1__: use a zcash style generator
__0.1.0__: release

# Howto

## API docs

```
cargo doc --open
```

## Benchmarks

```
cargo bench
```

## Examples
Counting the number of constraints in group operations
```
cargo run --example constraint_count_bandersnatch
cargo run --example constraint_count_jubjub
cargo run --example constraint_count_bandersnatch_glv
```