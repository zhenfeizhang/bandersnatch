bandersnatch-rust
------


This is a reference implementation of [Bendersnatch curve](https://ethresear.ch/t/introducing-bandersnatch-a-fast-elliptic-curve-built-over-the-bls12-381-scalar-field/9957) using [Arkwork](https://github.com/arkworks-rs/curves)'s framework in Rust.
The spec of the curve is available 
[here](https://github.com/asanso/Bandersnatch/blob/main/README.md).
There was also a Python reference implementation [here](https://github.com/asanso/Bandersnatch/).


# Logistics

- This code is released under MIT license.
- This code is not audited and may contain severe security flaws. Use at your own risk.
- Version 0.1.0.
- This repo is upstreamed to Arkworks [curve](https://github.com/arkworks-rs/curves/) crate.
- Todos:
    - [x] GVL multiplication
    - [ ] R1CS for GLV multiplication
    - [ ] Documentation coverage
    - [ ] Test coverage
    - [x] Update benchmark data for Benchmark section.

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
```