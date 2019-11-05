# sageCurves
## Overview
A sage script to compute all parameters needed to add a new elliptic curve to [libff](https://github.com/EYBlockchain/zk-swap-libff). The script support M- and D- type curves of embedding degree `k=12` constructed as tower fields `Fq_12` over `Fq_6` over `Fq_2`. Families include Barreto-Naehrig (BN), Barreto-Lynn-Scott (BLS12) and Brezing-Weng (BW12).

## Requirements
sageCurves requires a working [SageMath](http://www.sagemath.org) installation, and has been tested on SageMath version 8.8.

## Tutorial
```python
sage -python curve_parameters_libff.py path/to/curve_file.txt > path/to/output_dir/curve_name_init.cpp
```
`curve_file.txt` should be written as follow: 
```
type:curve_name
a
b
u
```
Where: 
+ `type` specifies the curve family (`bn`, `bls12`, `bw12`)
+ `curve_name` e.g.: `alt_bn128` or `bls12_381`
+ `a` curve coefficient as in short Weierstrass equation: `y^2=x^3+a*x+b`
+ `b` curve coefficient as in short Weierstrass equation: `y^2=x^3+a*x+b`
+ `u` from the polynomials `q(u)` and `r(u)` for `bn`, `bls12` and `bw12` families

### Examples
```python
sage -python curve_parameters_libff.py Curves/bls12_377.txt > ./bls12_377_init.cpp
```
This repo was used to add 8 curves to libff [here](https://github.com/EYBlockchain/zk-swap-libff/tree/ey/libff/algebra/curves)

