# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

from sage.all_cmdline import *   # import sage library
import sys

def q_bn(u):
    return Integer(36*u**4+36*u**3+24*u**2+6*u+1)

def r_bn(u):
    return Integer(36*u**4+36*u**3+18*u**2+6*u+1)

def q_bls12(u):
    return Integer((u-1)**2 * ((u**4-u**2+1)/3) + u)

def r_bls12(u):
    return Integer(u**4-u**2+1)

def q_bw12(u):
    return Integer(1728*x**6+2160*x**5+1548*x**4+756*x**3+240*x**2+54*x+7) 

def r_bw12(u):
    return Integer(36*x**4+36*x**3+18*x**2+6*x+1) 

def next_power_of_2(x):  
     return 1  if x == 0  else int(log(2 **(x - 1 ).nbits())/log(2 ))

def round_up(num, word):
     return num if num % word == 0  else num + word - num % word

def montgomery(p):
    size64 = round_up(next_power_of_2(p), 64 )
    R64 = 2 **size64 
    R2_64 = mod(R64**2 , p)
    R3_64 = mod(R64**3 , p)
    inv64 = hex(-1 /p % 2 **64 )

    size32 = round_up(next_power_of_2(p), 32 )
    R32 = 2 **size32 
    R2_32 = mod(R32**2 , p)
    R3_32 = mod(R32**3 , p)
    inv32 = hex(-1 /p % 2 **32 )

    return R2_64, R3_64, inv64, R2_32, R3_32, inv32  

def twoAdicity(num):
    if num % 2  != 0 : return 1 
    factor = 0 
    while num % 2  == 0 :
        num /= 2 
        factor += 1 
    return factor

def parameters_Fp(modulus):
    Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4 = montgomery(modulus) 
    num_bits = modulus.nbits()
    euler = int((modulus-1 )/2 )
    s = twoAdicity(modulus-1 )
    assert( (modulus-1 )%2 **s == 0 )
    t = int((modulus-1 )/2 **s)
    t_minus_1_over_2 = int((t-1 )/2 )
    multiplicative_generator = primitive_root(modulus)
    root_of_unity = power_mod(multiplicative_generator, t, modulus)
    " not necessiraly the least, can be chosen manually "
    nqr = least_quadratic_nonresidue(modulus)
    nqr_to_t = power_mod(nqr, t, modulus)
    return modulus, Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t

def parameters_Fp2(modulus):
    euler = int((modulus**2 -1 )/2 )
    s = twoAdicity(modulus**2 -1 )
    t = int((modulus**2 -1 )/2 **s)
    t_minus_1_over_2 = int((t-1 )/2 )
    " not necessiraly the least, can be chosen manually "
    non_residue = least_quadratic_nonresidue(modulus)
    Fq = GF(modulus)
    P = Fq['z']; (z,) = P._first_ngens(1)
    assert((z**2 -non_residue).is_irreducible())
    Fq2 = GF(modulus**2 , modulus=z**2 -non_residue, names=('u',)); (u,) = Fq2._first_ngens(1)
    PP = Fq2['zz']; (zz,) = PP._first_ngens(1)
    " find quadratic non-residue "
    " not necessiraly the least, can be chosen manually "
    nqr = u
    while(not (zz**2 -nqr).is_irreducible()):
        nqr += 1 
    nqr_to_t = nqr**t
    Frobenius_coeff_C1 = lift(Mod(non_residue, modulus)**((modulus**0 -1 )/2 )), lift(Mod(non_residue, modulus)**((modulus**1 -1 )/2 ))
    return euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1[0 ], Frobenius_coeff_C1[1 ]

def parameters_Fp6(modulus, non_residue, coeff_b, r):
    Fq = GF(modulus)
    P = Fq['z']; (z,) = P._first_ngens(1)
    Fq2 = GF(modulus**2 , modulus=z**2 -non_residue, names=('u',)); (u,) = Fq2._first_ngens(1)
    PP = Fq2['zz']; (zz,) = PP._first_ngens(1)
    nqr = u
    while(not (zz**3 -nqr).is_irreducible() or not (zz**6 -nqr).is_irreducible()):
        nqr += 1 
    Frobenius_coeffs_c1 = []
    Frobenius_coeffs_c2 = []
    mul_by_q = []
    for i in range(6 ):
        Frobenius_coeffs_c1.append(nqr**int((modulus**i-1 )/3 ))
    for i in range(6 ):
        Frobenius_coeffs_c2.append(nqr**int((2 *modulus**i-2 )/3 ))
    "compute mul_by_q for later"
    mul_by_q.append(Frobenius_coeffs_c1[2 ])
    mul_by_q.append(nqr**int((modulus**i-1 )/2 ))
    "compute coeff_b_twist for later"
    coeff_b_twist = coeff_b / nqr
    Et = EllipticCurve(Fq2, [0 , coeff_b_twist])
    if (Et.order() % r != 0 ): 
        coeff_b_twist = coeff_b * nqr
        twist_type = 'M'
        Et = EllipticCurve(Fq2, [0 , coeff_b_twist])
    else:
        twist_type = 'D'
    assert(Et.order() % r == 0 )
    "compute G1 generator for later"
    E = EllipticCurve(Fq, [0 , coeff_b])
    G1_one = E.random_point() # BN is prime order. TODO: find the simplest generator 
    trace = modulus+1 -E.order()
    "compute G2 generator for later"
    Et = EllipticCurve(Fq2, [0 , coeff_b_twist])
    G2_one = Et.random_point()
    return nqr, Frobenius_coeffs_c1, Frobenius_coeffs_c2, mul_by_q, coeff_b_twist, twist_type, G1_one, G2_one, trace

def parameters_Fp12(modulus, nqr):
    Frob12_coeff_C1 = []
    for i in range(12 ):
        Frob12_coeff_C1.append(nqr**int((modulus**i-1 )/6 ))
    return nqr, Frob12_coeff_C1

def print_header(curve_name):
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_g1.hpp>')
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_g2.hpp>')
    print('#include <libff/algebra/curves/' + curve_name + '/' + curve_name + '_init.hpp>\n')
    print('namespace libff {\n')
    print('\tbigint<' + curve_name + '_r_limbs> ' + curve_name + '_modulus_r;')
    print('\tbigint<' + curve_name + '_q_limbs> ' + curve_name + '_modulus_q;\n')
    print('\t' + curve_name + '_Fq ' + curve_name + '_coeff_b;')
    print('\t' + curve_name + '_Fq2 ' + curve_name + '_twist;')
    print('\t' + curve_name + '_Fq2 ' + curve_name + '_twist_coeff_b;')
    print('\t' + curve_name + '_Fq ' + curve_name + '_twist_mul_by_b_c0;')
    print('\t' + curve_name + '_Fq ' + curve_name + '_twist_mul_by_b_c1;')
    print('\t' + curve_name + '_Fq2 ' + curve_name + '_twist_mul_by_q_X;')
    print('\t' + curve_name + '_Fq2 ' + curve_name + '_twist_mul_by_q_Y;\n')
    print('\tbigint<' + curve_name + '_q_limbs> ' + curve_name + '_ate_loop_count;')
    print('\tbool ' + curve_name + '_ate_is_loop_count_neg;')
    print('\tbigint<12*' + curve_name + '_q_limbs> ' + curve_name + '_final_exponent;') # assumes k=12
    print('\tbigint<' + curve_name + '_q_limbs> ' + curve_name + '_final_exponent_z;')
    print('\tbool ' + curve_name + '_final_exponent_is_z_neg;\n')
    print('\tvoid init_' + curve_name + '_params()')
    print('\t{')
    print('\ttypedef bigint<' + curve_name + '_r_limbs> bigint_r;')
    print('\ttypedef bigint<' + curve_name + '_q_limbs> bigint_q;\n')
    print('\tassert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4);')

def print_Fp_parameters(curve_name, field_name, parameters):
    print('\t\n /* ' + curve_name + ' ' + field_name + ' parameters */\n')
    print('\t' + curve_name + '_modulus_' + field_name[-1 ] + ' = bigint_' + field_name[-1 ] + '("' + str(parameters[0 ]) + '");')  
    print('\tassert(' + curve_name + '_' + field_name + '::modulus_is_valid());')
    print('\tif (sizeof(mp_limb_t) == 8)')
    print('\t{');
    print('\t\t' + curve_name + '_' + field_name + '::Rsquared = bigint_' + field_name[-1 ] + '("' + str(parameters[1 ]) + '");')  
    print('\t\t' + curve_name + '_' + field_name + '::Rcubed = bigint_' + field_name[-1 ] + '("' + str(parameters[2 ]) + '");')  
    print('\t\t' + curve_name + '_' + field_name + '::inv = 0x' + str(parameters[3 ]) + ';')  
    print('\t}');
    print('\tif (sizeof(mp_limb_t) == 4)')
    print('\t{');
    print('\t' + curve_name + '_' + field_name + '::Rsquared = bigint_' + field_name[-1 ] + '("' + str(parameters[4 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::Rcubed = bigint_' + field_name[-1 ] + '("' + str(parameters[5 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::inv = 0x' + str(parameters[6 ]) + ';')  
    print('\t}');
    print('\t' + curve_name + '_' + field_name + '::num_bits = ' + str(parameters[7 ]) + ';')  
    print('\t' + curve_name + '_' + field_name + '::euler = bigint_' + field_name[-1 ] + '("' + str(parameters[8 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::s = ' + str(parameters[9 ]) + '; ')  
    print('\t' + curve_name + '_' + field_name + '::t = bigint_' + field_name[-1 ] + '("' + str(parameters[10 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::t_minus_1_over_2 = bigint_' + field_name[-1 ] + '("' + str(parameters[11 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::multiplicative_generator = ' + curve_name + '_' + field_name + '("' + str(parameters[12 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::root_of_unity = ' + curve_name + '_' + field_name + '("' + str(parameters[13 ]) + '");')  
    print('\t' + curve_name + '_' + field_name + '::nqr = ' + curve_name + '_' + field_name + '("' + str(parameters[14 ]) + '"); ')  
    print('\t' + curve_name + '_' + field_name + '::nqr_to_t = ' + curve_name + '_' + field_name + '("' + str(parameters[15 ]) + '");')

def print_Fp2_parameters(curve_name, parameters):
    print('\t\n /* ' + curve_name + ' Fq2' + ' parameters */\n')
    print('\t'  +curve_name + '_Fq2::euler = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[0 ]) + '");')
    print('\t'  +curve_name + '_Fq2::s = ' + str(parameters[1 ]) + ';')
    print('\t'  +curve_name + '_Fq2::t = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[2 ]) + '");')
    print('\t'  +curve_name + '_Fq2::t_minus_1_over_2 = bigint<2*' + curve_name + '_q_limbs>("' + str(parameters[3 ]) + '");')
    print('\t'  +curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq("' + str(parameters[4 ]) + '");')
    print('\t'  +curve_name + '_Fq2::nqr = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[5 ].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[5 ].polynomial().list())==1  else parameters[5 ].polynomial().list()[1 ]) + '"));')
    print('\t'  +curve_name + '_Fq2::nqr_to_t = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[6 ].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[6 ].polynomial().list())==1  else parameters[6 ].polynomial().list()[1 ]) + '"));')
    print('\t'  +curve_name + '_Fq2::Frobenius_coeffs_c1[0] = ' + curve_name + '_Fq("' + str(parameters[7 ]) + '");')
    print('\t'  +curve_name + '_Fq2::Frobenius_coeffs_c1[1] = ' + curve_name + '_Fq("' + str(parameters[8 ]) + '");')

def print_Fp6_parameters(curve_name, parameters):
    print('\t\n /* ' + curve_name + ' Fq6' + ' parameters */\n')
    print('\t' + curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[0 ].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[0 ].polynomial().list())==1  else parameters[0 ].polynomial().list()[1 ]) + '"));')
    for i in range(6 ):
        print('\t' + curve_name + '_Fq2::Frobenius_coeffs_c1[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[1 ][i].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[1 ][i].polynomial().list())==1  else parameters[1 ][i].polynomial().list()[1 ]) + '"));')
    for i in range(6 ):
        print('\t' + curve_name + '_Fq2::Frobenius_coeffs_c2[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[2 ][i].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[2 ][i].polynomial().list())==1  else parameters[2 ][i].polynomial().list()[1 ]) + '"));')

def print_Fp12_parameters(curve_name, parameters):
    print('\t\n /* ' + curve_name + ' Fq12' + ' parameters */\n')
    print('\t' + curve_name + '_Fq2::non_residue = ' + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[0 ].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[0 ].polynomial().list())==1  else parameters[0 ].polynomial().list()[1 ]) + '"));')
    for i in range(12 ):
        print('\t' + curve_name + '_Fq2::Frobenius_coeffs_c1[' + str(i) + '] = '  + curve_name + '_Fq2(' + curve_name + '_Fq("' + str(parameters[1 ][i].polynomial().list()[0 ]) + '"),' + curve_name + '_Fq("' + str(0  if len(parameters[1 ][i].polynomial().list())==1  else parameters[1 ][i].polynomial().list()[1 ]) + '"));')
