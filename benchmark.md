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

# \# constraints

| | \# constraints for 1 group op|
|:---|---| 
| Jubjub |3325|
| Bandersantch without GLV |  3325 |
| Bandersantch with GLV| 1468 |

# Edward curves ops

|   | bandersnatch (without GLV)| bandersnatch (without GLV)| Jubjub | ed_on_bls12_377|
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

