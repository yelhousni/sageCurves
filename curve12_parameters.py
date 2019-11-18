# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

"""
Computes parameters necessary for the implementation
of an elliptic curve from the families BN, BLS12 and BW12

For these families: 
    * embedding degree k=12
    * complex multiplication discriminant D=-3
"""

from utils import *
import sys, getopt

def help():
    print('')
    print('usage: sage' + sys.argv[0] + ' -i <input-file> -o <output-directory> -l <target-library>')
    print('e.g. : sage' + sys.argv[0] + ' -i ./Curves/alt_bn128.txt -o alt_bn128/ -l libff')
    print('\n options:')
    print('\t-h, --help: prints this help message')
    print('\t-i, --infile: input file (supported families: bn, bls12 and bw12)')
    print('\t-o, --outdir: output directory')
    print('\t-l, --lib: libff (default), bellman, zexe, py_ecc')
    print('')
	
def main(argv):
    curveFile = ''
    outdir = ''
    lib = 'libff'
    try:
        opts, args = getopt.getopt(argv,"hi:o:l:",["infile=","outdir=","lib="])
    except getopt.GetoptError:
       sys.exit(help())
    for opt, arg in opts:
       if opt == '-h':
          sys.exit(help())
       elif opt in ("-i", "--infile"):
          curveFile = arg
       elif opt in ("-o", "--outdir"):
          outdir = arg
       elif opt in ("-l", "--lib"):
          lib = arg
    try:
        with open(curveFile, 'r') as f:
            lines = f.readlines()
            curve_family = lines[0].split(":")[0]
            curve_name = lines[0].split(":")[1][:-1]
            poly_u = Integer(lines[1][:-1])

        # Curve
        modulus_q = get_q(curve_family, poly_u)
        frobenius_trace = get_t(curve_family, poly_u)
        modulus_r = get_r(curve_family, poly_u)
        embedding_degree = Integer(12)
        CM_discriminant = Integer(-3)
        elliptic_curve = small_B_twist(make_curve(modulus_q,frobenius_trace,modulus_r,embedding_degree,CM_discriminant))
        coeff_a = elliptic_curve.a4()
        coeff_b = elliptic_curve.a6()
            
        # Fr parameters
        R2_64_r, R3_64_r, inv_64_r, R2_32_r, R3_32_r, inv_32_r, num_bits_r, euler_r, s_r, t_r, t_minus_1_over_2_r, multiplicative_generator_r, root_of_unity_r, nqr_r, nqr_to_t_r = parameters_Fp(modulus_r)

        # Fq parameters
        R2_64_q, R3_64_q, inv_64_q, R2_32_q, R3_32_q, inv_32_q, num_bits_q, euler_q, s_q, t_q, t_minus_1_over_2_q, multiplicative_generator_q, root_of_unity_q, nqr_q, nqr_to_t_q = parameters_Fp(modulus_q)

        # Fq2 parameters 
        euler_q2, s_q2, t_q2, t_minus_1_over_2_q2, non_residue_q2, nqr_q2, nqr_to_t_q2, frobenius_q2 = parameters_Fp2(modulus_q)

        # Fq6 parameters
        non_residue_q6, frobenius_q6_1, frobenius_q6_2, mul_by_q, twist_coeff_b, twist_type, G1_one, G2_one = parameters_Fp6(modulus_q, non_residue_q2, coeff_b, modulus_r)

        # Fq12 parameters
        non_residue_q12, frobenius_q12 = parameters_Fp12(modulus_q, non_residue_q6)

        # pairing parameters
        ate_loop_count, ate_loop_count_bool, final_exponent, x_bool = pairing_parameters(curve_family, poly_u, modulus_q, embedding_degree, modulus_r)

        # print to libff
        fill_in_libff_templates(outdir, curve_name, curve_family, modulus_r, R2_64_r, R3_64_r, inv_64_r, R2_32_r, R3_32_r, inv_32_r, num_bits_r, euler_r, s_r, t_r, t_minus_1_over_2_r, multiplicative_generator_r, root_of_unity_r, nqr_r, nqr_to_t_r, modulus_q, R2_64_q, R3_64_q, inv_64_q, R2_32_q, R3_32_q, inv_32_q, num_bits_q, euler_q, s_q, t_q, t_minus_1_over_2_q, multiplicative_generator_q, root_of_unity_q, nqr_q, nqr_to_t_q, euler_q2, s_q2, t_q2, t_minus_1_over_2_q2, non_residue_q2, nqr_q2, nqr_to_t_q2, frobenius_q2, non_residue_q6, frobenius_q6_1, frobenius_q6_2, mul_by_q, coeff_b, twist_coeff_b, twist_type, G1_one, G2_one, non_residue_q12, frobenius_q12, ate_loop_count, ate_loop_count_bool, final_exponent, poly_u, x_bool)

    except Exception as e:
        sys.exit(e)
    return 0

if __name__ == '__main__':
    main(sys.argv[1:])
