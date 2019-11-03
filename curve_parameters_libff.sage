# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

"""
Computes BN curve parameters necessary for libff implementation
"""

import sys

def next_power_of_2(x):  
     return 1 if x == 0 else int(log(2**(x - 1).nbits())/log(2))

def round_up(num, word):
     return num if num % word == 0 else num + word - num % word

def montgomery(p):
    size64 = round_up(next_power_of_2(p), 64)
    R64 = 2^size64 
    R2_64 = mod(R64^2, p)
    R3_64 = mod(R64^3, p)
    inv64 = hex(-1/p % 2^64)

    size32 = round_up(next_power_of_2(p), 32)
    R32 = 2^size32 
    R2_32 = mod(R32^2, p)
    R3_32 = mod(R32^3, p)
    inv32 = hex(-1/p % 2^32)

    return R2_64, R3_64, inv64, R2_32, R3_32, inv32  

def twoAdicity(num):
    if num % 2 != 0: return 1
    factor = 0
    while num % 2 == 0:
        num /= 2
        factor += 1
    return factor

def parameters_Fp(modulus):
    Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4 = montgomery(modulus) 
    num_bits = modulus.nbits()
    euler = int((modulus-1)/2)
    s = twoAdicity(modulus-1)
    assert( (modulus-1)%2^s == 0)
    t = int((modulus-1)/2^s)
    t_minus_1_over_2 = int((t-1)/2)
    multiplicative_generator = primitive_root(modulus)
    root_of_unity = power_mod(multiplicative_generator, t, modulus)
    " not necessiraly the least, can be chosen manually "
    nqr = least_quadratic_nonresidue(modulus)
    nqr_to_t = power_mod(nqr, t, modulus)
    return modulus, Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t

def parameters_Fp2(modulus):
    euler = int((modulus^2-1)/2)
    s = twoAdicity(modulus^2-1)
    t = int((modulus^2-1)/2^s)
    t_minus_1_over_2 = int((t-1)/2)
    " not necessiraly the least, can be chosen manually "
    non_residue = least_quadratic_nonresidue(modulus)
    Fq = GF(modulus)
    P.<z> = Fq[]
    assert((z^2-non_residue).is_irreducible())
    Fq2.<u> = GF(modulus^2, modulus=z^2-non_residue)
    PP.<zz> = Fq2[]
    " find quadratic non-residue "
    " not necessiraly the least, can be chosen manually "
    nqr = u
    while(not (zz^2-nqr).is_irreducible()):
        nqr += 1
    nqr_to_t = nqr^t
    Frobenius_coeff_C1 = lift(Mod(non_residue, modulus)^((modulus^0-1)/2)), lift(Mod(non_residue, modulus)^((modulus^1-1)/2))
    return euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1[0], Frobenius_coeff_C1[1]

def parameters_Fp6(modulus, non_residue, coeff_b, r):
    Fq = GF(modulus)
    P.<z> = Fq[]
    Fq2.<u> = GF(modulus^2, modulus=z^2-non_residue)
    PP.<zz> = Fq2[]
    nqr = u
    while(not (zz^3-nqr).is_irreducible() or not (zz^6-nqr).is_irreducible()):
        nqr += 1
    Frobenius_coeffs_c1 = []
    Frobenius_coeffs_c2 = []
    mul_by_q = []
    for i in range(6):
        Frobenius_coeffs_c1.append(nqr^int((modulus^i-1)/3))
    for i in range(6):
        Frobenius_coeffs_c2.append(nqr^int((2*modulus^i-2)/3))
    "compute mul_by_q for later"
    mul_by_q.append(Frobenius_coeffs_c1[2])
    mul_by_q.append(nqr^int((modulus^i-1)/2))
    "compute coeff_b_twist for later"
    coeff_b_twist = coeff_b / nqr
    Et = EllipticCurve(Fq2, [0, coeff_b_twist])
    if (Et.order() % r != 0): 
        coeff_b_twist = coeff_b * nqr
        twist_type = 'M'
        Et = EllipticCurve(Fq2, [0, coeff_b_twist])
    else:
        twist_type = 'D'
    assert(Et.order() % r == 0)
    "compute G1 generator for later"
    E = EllipticCurve(Fq, [0, coeff_b])
    G1_one = E.random_point() # BN is prime order. TODO: find the simplest generator 
    trace = modulus+1-E.order()
    "compute G2 generator for later"
    Et = EllipticCurve(Fq2, [0, coeff_b_twist])
    G2_one = Et.random_point()
    return nqr, Frobenius_coeffs_c1, Frobenius_coeffs_c2, mul_by_q, coeff_b_twist, twist_type, G1_one, G2_one, trace

def parameters_Fp12(modulus, nqr):
    Frob12_coeff_C1 = []
    for i in range(12):
        Frob12_coeff_C1.append(nqr^int((modulus^i-1)/6))
    return nqr, Frob12_coeff_C1

def tab_print(string):
    print('\t' + string)

def print_header(curve_name):
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_g1.hpp>')
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_g2.hpp>')
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_init.hpp>\n')
    print('namespace libff {\n')
    print('bigint<' + curve_name + '_r_limbs> ' + curve_name + '_modulus_r;')
    print('bigint<' + curve_name + '_q_limbs> ' + curve_name + '_modulus_q;\n')
    print(curve_name + '_Fq ' + curve_name + '_coeff_b;')
    print(curve_name + '_Fq2 ' + curve_name + '_twist;')
    print(curve_name + '_Fq2 ' + curve_name + '_twist_coeff_b;')
    print(curve_name + '_Fq ' + curve_name + '_twist_mul_by_b_c0;')
    print(curve_name + '_Fq ' + curve_name + '_twist_mul_by_b_c1;')
    print(curve_name + '_Fq2 ' + curve_name + '_twist_mul_by_q_X;')
    print(curve_name + '_Fq2 ' + curve_name + '_twist_mul_by_q_Y;\n')
    print('bigint<' + curve_name + '_q_limbs> ' + curve_name + '_ate_loop_count;')
    print('bool ' + curve_name + '_ate_is_loop_count_neg;')
    print('bigint<12*' + curve_name + '_q_limbs> ' + curve_name + '_final_exponent;') # assumes k=12
    print('bigint<' + curve_name + '_q_limbs> ' + curve_name + '_final_exponent_z;')
    print('bool ' + curve_name + '_final_exponent_is_z_neg;\n')
    print('void init_' + curve_name + '_params()')
    print('{')
    tab_print('typedef bigint<' + curve_name + '_r_limbs> bigint_r;')
    tab_print('typedef bigint<' + curve_name + '_q_limbs> bigint_q;\n')
    tab_print('assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4);')

def print_Fp_parameters(curve_name, field_name, parameters):
    tab_print('\n /* ' + curve_name + ' ' + field_name + ' parameters */\n')
    tab_print(curve_name + '_modulus_' + field_name[-1] + ' = bigint_' + field_name[-1] + '("' + str(parameters[0]) + '");')  
    tab_print('assert(' + curve_name + '_' + field_name + '::modulus_is_valid());')
    tab_print('if (sizeof(mp_limb_t) == 8)')
    tab_print('{');
    tab_print('\t' + curve_name + '_' + field_name + '::Rsquared = bigint_' + field_name[-1] + '("' + str(parameters[1]) + '");')  
    tab_print('\t' + curve_name + '_' + field_name + '::Rcubed = bigint_' + field_name[-1] + '("' + str(parameters[2]) + '");')  
    tab_print('\t' + curve_name + '_' + field_name + '::inv = 0x' + str(parameters[3]) + ';')  
    tab_print('}');
    tab_print('if (sizeof(mp_limb_t) == 4)')
    tab_print('{');
    tab_print('\t' + curve_name + '_' + field_name + '::Rsquared = bigint_' + field_name[-1] + '("' + str(parameters[4]) + '");')  
    tab_print('\t' + curve_name + '_' + field_name + '::Rcubed = bigint_' + field_name[-1] + '("' + str(parameters[5]) + '");')  
    tab_print('\t' + curve_name + '_' + field_name + '::inv = 0x' + str(parameters[6]) + ';')  
    tab_print('}');
    tab_print(curve_name + '_' + field_name + '::num_bits = ' + str(parameters[7]) + ';')  
    tab_print(curve_name + '_' + field_name + '::euler = bigint_' + field_name[-1] + '("' + str(parameters[8]) + '");')  
    tab_print(curve_name + '_' + field_name + '::s = ' + str(parameters[9]) + '; ')  
    tab_print(curve_name + '_' + field_name + '::t = bigint_' + field_name[-1] + '("' + str(parameters[10]) + '");')  
    tab_print(curve_name + '_' + field_name + '::t_minus_1_over_2 = bigint_' + field_name[-1] + '("' + str(parameters[11]) + '");')  
    tab_print(curve_name + '_' + field_name + '::multiplicative_generator = ' + curve_name + '_' + field_name + '("' + str(parameters[12]) + '");')  
    tab_print(curve_name + '_' + field_name + '::root_of_unity = ' + curve_name + '_' + field_name + '("' + str(parameters[13]) + '");')  
    tab_print(curve_name + '_' + field_name + '::nqr = ' + curve_name + '_' + field_name + '("' + str(parameters[14]) + '"); ')  
    tab_print(curve_name + '_' + field_name + '::nqr_to_t = ' + curve_name + '_' + field_name + '("' + str(parameters[15]) + '");')

def print_Fp2_parameters(curve_name, parameters):
    tab_print('\n /* ' + curve_name + ' Fq2' + ' parameters */\n')
    tab_print(curve_name + '_Fq2::euler = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[0]) + '");')
    tab_print(curve_name + '_Fq2::s = ' + str(parameters[1]) + ';')
    tab_print(curve_name + '_Fq2::t = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[2]) + '");')
    tab_print(curve_name + '_Fq2::t_minus_1_over_2 = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[3]) + '");')
    tab_print(curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq("' + str(parameters[4]) + '");')
    tab_print(curve_name + '_Fq2::nqr = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[5].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[5].polynomial().list())==1 else parameters[5].polynomial().list()[1]) + '"));')
    tab_print(curve_name + '_Fq2::nqr_to_t = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[6].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[6].polynomial().list())==1 else parameters[6].polynomial().list()[1]) + '"));')
    tab_print(curve_name + '_Fq2::Frobenius_coeffs_c1[0] = ' + curve_name + '_Fq("' + str(parameters[7]) + '");')
    tab_print(curve_name + '_Fq2::Frobenius_coeffs_c1[1] = ' + curve_name + '_Fq("' + str(parameters[8]) + '");')

def print_Fp6_parameters(curve_name, parameters):
    tab_print('\n /* ' + curve_name + ' Fq6' + ' parameters */\n')
    tab_print(curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[0].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[0].polynomial().list())==1 else parameters[0].polynomial().list()[1]) + '"));')
    for i in range(6):
        tab_print(curve_name + '_Fq2::Frobenius_coeffs_c1[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[1][i].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[1][i].polynomial().list())==1 else parameters[1][i].polynomial().list()[1]) + '"));')
    for i in range(6):
        tab_print(curve_name + '_Fq2::Frobenius_coeffs_c2[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[2][i].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[2][i].polynomial().list())==1 else parameters[2][i].polynomial().list()[1]) + '"));')

def print_Fp12_parameters(curve_name, parameters):
    tab_print('\n /* ' + curve_name + ' Fq12' + ' parameters */\n')
    tab_print(curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[0].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[0].polynomial().list())==1 else parameters[0].polynomial().list()[1]) + '"));')
    for i in range(12):
        tab_print(curve_name + '_Fq2::Frobenius_coeffs_c1[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[1][i].polynomial().list()[0]) + '"),' + curve_name + '_Fq("' + str(0 if len(parameters[1][i].polynomial().list())==1 else parameters[1][i].polynomial().list()[1]) + '"));')
   
def help():
    print('')
    print('usage: ' + sys.argv[0] + ' <curve_file>')
    print('e.g. : ' + sys.argv[0] + ' ./curves/alt_bn128.txt')
    print('\n supported curves:')
    tab_print('BN')
    tab_print('BLS12')
    tab_print('BW12')
    tab_print('SW12')
    print('\n TODO:')
    tab_print('SW6')
    tab_print('MNT6')
    tab_print('MNT4')
    print('')
	
def main():
    if (len(sys.argv) != 2): sys.exit(help())
    try:
        curveFile = sys.argv[1]
        with open(curveFile, 'r') as f:
            lines = f.readlines()
            family = lines[0].split(":")[0]
            curve = lines[0].split(":")[1][:-1]
            coeff_a = Integer(lines[1])
            coeff_b = Integer(lines[2]) 
            r = Integer(lines[3])
            q = Integer(lines[4])
            x =  Integer(lines[5])

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
        tab_print('\n/* choice of short Weierstrass curve and its twist */\n')
        tab_print(curve + '_coeff_b = ' + curve + '_Fq("' + str(coeff_b) + '");')
        tab_print(curve + '_twist_coeff_b = ' + curve + '_Fq2(' + curve + '_Fq("' + str(coeff_b_twist.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0 if len(coeff_b_twist.polynomial().list()) == 1 else coeff_b_twist.polynomial().list()[1]) + '"); // ' + twist_type + '-type curve')
        tab_print(curve + '_twist_mul_by_c0 = ' + curve + "coeff_b * " + curve + "_Fq2::non_residue;") 
        tab_print(curve + '_twist_mul_by_c1 = ' + curve + "coeff_b * " + curve + "_Fq2::non_residue;") 
        tab_print(curve + '_twist_mul_by_q_X = ' + curve + '_Fq2(' + curve + '_Fq("' + str(mul_by_q_x.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0 if len(mul_by_q_x.polynomial().list()) == 1 else mul_by_q_x.polynomial().list()[1]) + '");')
        tab_print(curve + '_twist_mul_by_q_Y = ' + curve + '_Fq2(' + curve + '_Fq("' + str(mul_by_q_y.polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0 if len(mul_by_q_y.polynomial().list()) == 1 else mul_by_q_y.polynomial().list()[1]) + '");')

        # choice of group G1
        G1 = Fq6_params[6]
        tab_print('\n/* choice of group G1 */\n')
        tab_print(curve + '_G1::G1_zero = ' + curve + '_G1(' + curve + '_Fq::zero(),' + curve + '_Fq::one(),' + curve + '_Fq::zero());')
        tab_print(curve + '_G1::G1_one = ' + curve + '_G1(' + curve + '_Fq("' + str(G1[0]) + '"),' + curve + ' _Fq("' + str(G1[1]) + '"),' + curve + '_Fq::zero());')

        # choice of group G2
        G2 = Fq6_params[7] 
        tab_print('\n/* choice of group G2 */\n')
        tab_print(curve + '_G2::G2_zero = ' + curve + '_G2(' + curve + '_Fq2::zero(),' + curve + '_Fq2::one(),' + curve + '_Fq2::zero());')
        tab_print(curve + '_G2::G2_one = ' + curve + '_G2(' + curve + '_Fq2(' + curve + '_Fq("' + str(G2[0].polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0 if len(G2[0].polynomial().list()) == 1 else G2[0].polynomial().list()[1]) + '")),' + curve + '_Fq2(' + curve + '_Fq("' + str(G2[1].polynomial().list()[0]) + '"),' + curve + '_Fq("' + str(0 if len(G2[1].polynomial().list()) == 1 else G2[1].polynomial().list()[1]) + '")),' + curve + 'Fq2::one()));')

        # pairing parameters
        emb_degree = 12
        trace = Fq6_params[8]
        if (family == 'bn' or family == 'bw'):
            ate_loop_count = 6*x+2
        elif (family == 'bls'):
            ate_loop_count = x
        else: 
            ate_loop_count = trace-1
        if (ate_loop_count < 0): 
            ate_loop_is_neg = 'true'; 
        else:
            ate_loop_is_neg = 'false'; 
        final_exponent = (q^emb_degree-1)/r 
        if (x<0):
            final_exponent_is_z_neg = 'true'
        else:
            final_exponent_is_z_neg = 'false'
        tab_print('\n/* choice of pairing */\n')
        tab_print(curve + '_ate_loop_count = bigint_q("' + str(abs(ate_loop_count)) + '");')
        tab_print(curve + '_ate_is_loop_count_neg = bigint_q("' + ate_loop_is_neg + '");')
        tab_print(curve + '_final_exponent = bigint<12*' + curve + 'q_limbs>("' + str(final_exponent) + '");')
        tab_print(curve + '_final_exponent_z = bigint_q("' + str(abs(x)) + '");')
        tab_print(curve + '_final_exponent_is_z_neg = bigint_q("' + final_exponent_is_z_neg + '");')

    except Exception as e:
        sys.exit(0)
    return 0

if __name__ == "__main__":
    main()
