# sageCurves
## Overview
A sage script to compute all parameters needed to add a new elliptic curve to [libff](https://github.com/EYBlockchain/zk-swap-libff), [bellman](https://github.com/zkcrypto/bellman), [zexe](https://github.com/scipr-lab/zexe) or [py_ecc](https://github.com/ethereum/py_ecc) libraries (wip: works for *libff* for now). The script support M- and D- type curves of embedding degree `k=12` constructed as tower fields `Fq_12` over `Fq_6` over `Fq_2`. Families include Barreto-Naehrig (BN), Barreto-Lynn-Scott (BLS12) and Brezing-Weng (BW12).

## Requirements
sageCurves requires:
+ a working [SageMath](http://www.sagemath.org) installation (tested on SageMath version 8.8)
+ a working [pystache](https://github.com/defunkt/pystache) installation inside `sage` (Download the package and `sage -pip /path/to/package/`)

## Tutorial
```
sage curve12_parameters.py -i <curve-file> -o <output-directory> -l <target-library> 
```
```
options:
    -h, --help: prints this help message
    -i, --infile: input file (supported families: bn, bls12 and bw12) 
    -o, --outdir: output directory
    -l, --lib: libff (default), bellman, zexe, py_ecc
```
`<curve-file>` should be written as follows:
```
type:curve_name
u
```
Where: 
+ `type` specifies the curve family (`bn`, `bls12`, `bw12`)
+ `curve_name` e.g.: `alt_bn128` or `bls12_381`
+ `u` from the polynomials `q(u)` and `r(u)` for `bn`, `bls12` and `bw12` families

### Examples
```python
sage curve_parameters_libff.py -i eg/bls12_377.txt -o ./bls12_377_init.cpp -l libff
```
The output files generated under `<outdir>` are:
+ `bls12_377_init.cpp`
+ `bls12_377_init.hpp`
+ `bls12_377_g1.hpp`
+ `bls12_377_g1.cpp`
+ `bls12_377_g2.hpp`
+ `bls12_377_g2.cpp`
+ `bls12_377_pp.hpp`
+ `bls12_377_pp.cpp`
+ `bls12_377_pairing.hpp`
+ `bls12_377_pairing.cpp`
 
