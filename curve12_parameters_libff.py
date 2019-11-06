# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

"""
Computes parameters necessary for libff implementation
of an elliptic curve from the families BN, BLS12 and BW12

For these families: 
    * embedding degree k=12
    * complex multiplication discriminant D=-3
"""

from sage.all_cmdline import *   # import sage library
from utils import *

def help():
    print('')
    print('usage: ' + sys.argv[0] + ' <curve_file>')
    print('e.g. : ' + sys.argv[0] + ' ./curves/alt_bn128.txt')
    print('\n supported families:')
    print('\tBN')
    print('\tBLS12')
    print('\tBW12')
    print('')
	
def main():
    if (len(sys.argv) != 2): sys.exit(help())
    try:
        curveFile = sys.argv[1]
        with open(curveFile, 'r') as f:
            lines = f.readlines()
            family = lines[0].split(":")[0]
            curve = lines[0].split(":")[1][:-1]
            x = Integer(lines[1])

        # Curve
        q = get_q(family, x)
        t = get_t(family, x)
        r = get_r(family, x)
        k = Integer(12)
        D = Integer(-3)
        E = make_curve(q,t,r,k,D)
        coeff_a = E.a4()
        coeff_b = E.a6()
            
        # header
        print_header(curve)

        # Fr parameters
        Fr_params = parameters_Fp(r)
        print_Fp_parameters(curve, 'Fr', Fr_params)

        # Fq parameters
        Fq_params = parameters_Fp(q)
        print_Fp_parameters(curve, 'Fq', Fq_params)

        # Fq2 parameters 
        Fq2_params = parameters_Fp2(q)
        print_Fp2_parameters(curve, Fq2_params)

        # Fq6 parameters
        Fq6_params = parameters_Fp6(q, Fq2_params[4], coeff_b, r)
        print_Fp6_parameters(curve, Fq6_params)

        # Fq12 parameters
        Fq12_params = parameters_Fp12(q, Fq6_params[0])
        print_Fp12_parameters(curve, Fq12_params)

        # curve and its twist
        coeff_b_twist = Fq6_params[4]
        twist_type = Fq6_params[5]
        mul_by_q_x = Fq6_params[3][0]
        mul_by_q_y = Fq6_params[3][1]
        print('\t\n/* choice of short Weierstrass curve and its twist */\n')
        print('\t' + curve + '_coeff_b = ' + curve + '_Fq("' + str(coeff_b) + '");')
        print('\t' + curve + '_twist_coeff_b = ' + curve + '_Fq2(' + curve + '_Fq("' + str(coeff_b_twist.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0  if len(coeff_b_twist.polynomial().list()) == 1  else coeff_b_twist.polynomial().list()[1]) + '"); // ' + twist_type + '-type curve')
        print('\t' + curve + '_twist_mul_by_c0 = ' + curve + "coeff_b * " + curve + "_Fq2::non_residue;") 
        print('\t' + curve + '_twist_mul_by_c1 = ' + curve + "coeff_b * " + curve + "_Fq2::non_residue;") 
        print('\t' + curve + '_twist_mul_by_q_X = ' + curve + '_Fq2(' + curve + '_Fq("' + str(mul_by_q_x.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0  if len(mul_by_q_x.polynomial().list()) == 1  else mul_by_q_x.polynomial().list()[1]) + '");')
        print('\t' + curve + '_twist_mul_by_q_Y = ' + curve + '_Fq2(' + curve + '_Fq("' + str(mul_by_q_y.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0  if len(mul_by_q_y.polynomial().list()) == 1  else mul_by_q_y.polynomial().list()[1]) + '");')

        # choice of group G1
        G1 = Fq6_params[6]
        print('\t\n/* choice of group G1 */\n')
        print('\t' + curve + '_G1::G1_zero = ' + curve + '_G1(' + curve + '_Fq::zero(),' + curve + '_Fq::one(),' + curve + '_Fq::zero());')
        print('\t' + curve + '_G1::G1_one = ' + curve + '_G1(' + curve + '_Fq("' + str(G1[0]) + '"),' + curve + ' _Fq("' + str(G1[1]) + '"),' + curve + '_Fq::zero());')

        # choice of group G2
        G2 = Fq6_params[7] 
        print('\t\n/* choice of group G2 */\n')
        print('\t' + curve + '_G2::G2_zero = ' + curve + '_G2(' + curve + '_Fq2::zero(),' + curve + '_Fq2::one(),' + curve + '_Fq2::zero());')
        print('\t' + curve + '_G2::G2_one = ' + curve + '_G2(' + curve + '_Fq2(' + curve + '_Fq("' + str(G2[0].polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0  if len(G2[0].polynomial().list()) == 1  else G2[0].polynomial().list()[1]) + '")),' + curve + '_Fq2(' + curve + '_Fq("' + str(G2[1].polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0  if len(G2[1].polynomial().list()) == 1  else G2[1].polynomial().list()[1]) + '")),' + curve + 'Fq2::one()));')

        # pairing parameters
        if (family == 'bn' or family == 'bw12'):
            ate_loop_count = 6*x+2 
        elif (family == 'bls12'):
            ate_loop_count = x
        else: 
            raise Exception("family unknown or not supported")
        if (ate_loop_count < 0 ): 
            ate_loop_is_neg = 'true'; 
        else:
            ate_loop_is_neg = 'false'; 
        final_exponent = (q**k-1)/r 
        if (x<0 ):
            final_exponent_is_z_neg = 'true'
        else:
            final_exponent_is_z_neg = 'false'
        print('\t\n/* choice of pairing */\n')
        print('\t' + curve + '_ate_loop_count = bigint_q("' + str(abs(ate_loop_count)) + '");')
        print('\t' + curve + '_ate_is_loop_count_neg = bigint_q("' + ate_loop_is_neg + '");')
        print('\t' + curve + '_final_exponent = bigint<12*' + curve + 'q_limbs>("' + str(final_exponent) + '");')
        print('\t' + curve + '_final_exponent_z = bigint_q("' + str(abs(x)) + '");')
        print('\t' + curve + '_final_exponent_is_z_neg = bigint_q("' + final_exponent_is_z_neg + '");')

    except Exception as e:
        sys.exit(0 )
    return 0 

if __name__ == "__main__":
    main()

