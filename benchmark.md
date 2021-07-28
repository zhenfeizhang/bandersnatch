benchmarking cost for curve operation of Arkwork's implementation
------
AMD 5900x; ubuntu 20.04; arkworks 0.3.0; rust 1.52.1

# GLV micro benchmarks

|  ops | cost |
|:---|---| 
| | |
| endomorphism | 2.15 us |
| scalar decomposition| 0.75 us|
| multi-scalar-mult | 41.5 us|

# MSM with(out) GLV

| dim | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
|:---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| msm | 312 us | 355 us | 406 us | 545 us | 1.07 ms | 1.77 ms | 2.43 ms | 3.99 ms | 7.01 ms | 11.3 ms | 19.8 ms | 36.6 ms | 60.54 ms | 110.5 ms | 
| msm with GLV | 178 us | 227 us | 304 us | 599 us | 1.07 ms | 1.58 ms | 2.79 ms | 5.08 ms | 8.13 ms | 15.8 ms | 29.7 ms | 54.3 ms | 102.1 ms | 197.8 ms |


# Edward curves ops

|   | bandersnatch (with GLV)| bandersnatch (without GLV)| Jubjub | ed_on_bls12_377|
|:---|---| --- | ---|---|
| fix base mul | 44 us | 80 us | 75 us  | 73 us |

# G1 ops

|   |  bls12-381 | bls12_377 | blst-bls12-381 |
|:---|---| --- | --- |
| fix base mul | 123 us  | 121 us | 87 us |

# G2 ops

|   |  bls12-381 | bls12_377 | blst-bls12-381 |
|:---|---| --- | --- | 
| fix base mul | 372 us  | 432 us | 163 us | 

