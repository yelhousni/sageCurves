# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

from sage.all_cmdline import *   # import sage library
import sys

# Curve utils
def make_curve(q,t,r,k,D,debug=False):
    """
    credits to https://github.com/scipr-lab/ecfactory

    Description:
    
        Finds the curve equation for the elliptic curve (q,t,r,k,D) using the Complex Multiplication method
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        E - elliptic curve over F_q with trace t,
            a subgroup of order r with embedding degree k,
            and fundamental discriminant D
    
    """
    assert is_valid_curve(q,t,r,k,D), 'Invalid input. No curve exists.' # check inputs
    if debug:
        print('Tested input')
    poly = hilbert_class_polynomial(D) # compute hilbert class polynomial
    if debug:
        print('Computed Hilbert class polynomial')
    check = False
    j_inv = poly.any_root(GF(q)) # find j-invariant    
    orig_curve = EllipticCurve(GF(q), j=j_inv) # make a curve
    E = orig_curve
    check = test_curve(q,t,r,k,D,E) # see if this is the right curve
    twist = False
    if not check: # not the right curve, use quadratic twist
        E = E.quadratic_twist()
        check = test_curve(q,t,r,k,D,E)
        if check:
            twist = True
        else: # twist didnt work => j = 0 or 1728
            if j_inv == 0: # for j = 0, use sextic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 6:
                    E = orig_curve.sextic_twist(power_mod(prim,i,q))
                    i+=1
            elif j_inv == 1728: # for j = 1728, use quartic twists
                prim = primitive_root(q)
                i = 1
                while t != E.trace_of_frobenius() and i < 4:
                    E = orig_curve.quartic_twist(power_mod(prim,i,q))
                    i+=1
            else: # twist didnt work and j != 0, 1728. this should never happen, so write input to a file for debugging
                print('Error. Quadratic twist failed to find the correct curve with j != 0, 1728. Logging output to debug.txt') # this line should never be reached'
                f = open('debug.txt', 'w')
                f.write('Twist: ' + str(twist) + '\n')
                f.write('q: ' + str(q) + '\n')
                f.write('t: ' + str(t) + '\n')
                f.write('r: ' + str(r) + '\n')
                f.write('k: ' + str(k) + '\n')
                f.write('D: ' + str(D) + '\n')
                f.write('E: ' + str(E) + '\n')
                f.write('orig_curve: ' + str(orig_curve))
                f.close()
                return False
            check = test_curve(q,t,r,k,D,E)
            twist = True
    if not check: # didnt find a curve. this should never happen, so write input to a file for debugging
        print('Error. Failed to find curve. Logging output to debug.txt')
        f = open('debug.txt', 'w')
        f.write('Twist: ' + str(twist) + '\n')
        f.write('q: ' + str(q) + '\n')
        f.write('t: ' + str(t) + '\n')
        f.write('r: ' + str(r) + '\n')
        f.write('k: ' + str(k) + '\n')
        f.write('D: ' + str(D) + '\n')
        f.write('E: ' + str(E) + '\n')
        f.write('orig_curve: ' + str(orig_curve))
        f.close()
        return False
    return E

def small_B_twist(E):
    """
    Description:
        
        Finds a curve isogenous to E that has small B in the curve equation y^2 = x^3 + A*x + B
    
    Input:
    
        E - elliptic curve
    
    Output:
    
        E' - elliptic curve isogenous to E that has small B in the curve equation y^2 = x^3 + A*x + B
    
    """
    b = E.ainvs()[4]
    q = E.base_field().order()
    b = power_mod(Integer(b), -1, q)
    d = 0
    s = Mod(1,q)
    bool = True
    while bool:
        try:
            d = (s*b)
            d = d.nth_root(3)
            d = Integer(d)
            bool = False
        except ValueError as e:
            s+=1 
            pass
    ainvs = [i for i in E.ainvs()]
    ainvs[3] *= d**2
    ainvs[4] *= d**3
    return EllipticCurve(E.base_field(), ainvs)

def test_curve(q,t,r,k,D,E): 
    """
    Description:
    
       Tests that E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    Input:
    
        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant
    
    Output:
    
        bool - true iff E is an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    
    """    
    bool = True
    bool = bool and (power_mod(q, k, r) == 1) #q^k -1 ==0 mod r
    bool = bool and (E.trace_of_frobenius() == t)
    bool = bool and (kronecker((t*t-4*q) * Integer(D).inverse_mod(q),q) == 1)
    bool = bool and (E.cardinality() == q+1-t)
    bool = bool and (E.cardinality() % r ==0)
    return bool


def is_valid_curve(q,t,r,k,D):
    """
    Description:

        Tests that (q,t,r,k,D) is a valid elliptic curve

    Input:

        q - size of prime field
        t - trace of Frobenius
        r - size of prime order subgroup
        k - embedding degree
        D - (negative) fundamental discriminant

    Output:

        bool - true iff there exists an elliptic curve over F_q with trace t, a subgroup of order r with embedding degree k, and fundamental discriminant D
    """
    if q == 0 or t == 0 or r == 0 or k == 0 or D == 0:
        return False
    if not is_prime(q):
        return False
    if not is_prime(r):
        return False
    if not fundamental_discriminant(D) == D:
        return False
    if D % 4 == 0: #check CM equation
        if not is_square(4*(t*t - 4*q)//D):
            return False
    if D % 4 == 1:
        if not is_square((t*t - 4*q)//D):
            return False
    if not (q+1-t) % r == 0: #check r | #E(F_q)
        return False
    if not power_mod(q,k,r) == 1: #check embedding degree is k
        return False
    return True

def get_q(family, u):
    """
    Description:

        Computes the field size of a BN/BLS12/BW12 elliptic curve

    Input:

        family - pairing-friendly family (bn, bls12 or bw12) 
        u - curve integer

    Output:

        Integer - field size
    """
    if (family == 'bn'):
       return Integer(36*u**4+36*u**3+24*u**2+6*u+1)
    elif (family == 'bls12'):
        return Integer((u-1)**2 * ((u**4-u**2+1)/3) + u)
    elif (family == 'bw12'):
        return Integer(1728*u**6+2160*u**5+1548*u**4+756*u**3+240*u**2+54*u+7) 
    else:
        raise Exception('family not supported')

def get_r(family, u):
    """
    Description:

        Computes the subgroup order of a BN/BLS12/BW12 elliptic curve

    Input:

        family - pairing-friendly family (bn, bls12 or bw12) 
        u - curve integer

    Output:

        Integer - subgroup order
    """
    if (family == 'bn'):
        return Integer(36*u**4+36*u**3+18*u**2+6*u+1)
    elif (family == 'bls12'):
        return Integer(u**4-u**2+1)
    elif (family == 'bw12'):
        return Integer(36*u**4+36*u**3+18*u**2+6*u+1) 
    else:
        raise Exception('family not supported')

def get_t(family, u):
    """
    Description:

        Computes the Frebonius trace of a BN/BLS12/BW12 elliptic curve

    Input:

        family - pairing-friendly family (bn, bls12 or bw12) 
        u - curve integer

    Output:

        Integer - Frobenius trace
    """
    if (family == 'bn'):
        return Integer(6*u**2+1)
    elif (family == 'bls12'):
        return Integer(u+1)
    elif (family == 'bw12'):
        return Integer(-6*u**2+1)
    else:
        raise Exception('family not supported')

# Field utils
def next_power_of_2(x):  
     return 1 if x == 0 else int(log(2 **(x - 1).nbits())/log(2))

def round_up(num, word):
     return num if num % word == 0  else num + word - num % word

def montgomery(p):
    """
    Description:

        Computes parameters needed for Montgomery arithmetic

    Input:

        p - modulus

    Output:

        2x Integer - Montgomery constant squared for 64-arch and 32-arch
        2x Integer - Montgomery constant cubed for 64-arch and 32-arch
        2x hex - inverse of minus modulus for 64-arch and 32-arch
    """ 
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
    """
    Description:

        Computes the 2-adicity order of an integer

    Input:

        num - integer

    Output:

        Integer - 2-adicity order 
    """ 
    if num % 2  != 0 : return 1 
    factor = 0 
    while num % 2  == 0 :
        num /= 2 
        factor += 1 
    return factor

def parameters_Fp(modulus):
    """
    Description:

        Computes parameters of finite field F_modulus

    Input:

        mmodulus - integer

    Output:
        modulus, Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t
    """ 
    Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4 = montgomery(modulus) 
    num_bits = modulus.nbits()
    euler = int((modulus-1)/2)
    s = twoAdicity(modulus-1)
    assert( (modulus-1)%2 **s == 0 )
    t = int((modulus-1)/2 **s)
    t_minus_1_over_2 = int((t-1)/2)
    multiplicative_generator = primitive_root(modulus)
    root_of_unity = power_mod(multiplicative_generator, t, modulus)
    nqr = least_quadratic_nonresidue(modulus)
    nqr_to_t = power_mod(nqr, t, modulus)
    return modulus, Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t

def parameters_Fp2(modulus):
    """
    Description:

        Computes parameters of a quadratic extension finite field F_(modulus^2)
        F_(modulus^2)[u] = F_modulus / u^2 - non_residue

    Input:

        mmodulus - integer

    Output:
        euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1[0], Frobenius_coeff_C1[1]
    """ 
    euler = int((modulus**2-1)/2)
    s = twoAdicity(modulus**2-1)
    t = int((modulus**2-1)/2**s)
    t_minus_1_over_2 = int((t-1)/2)
    Fq = GF(modulus)
    if (modulus % 4 == 3):
        non_residue = Fq(-1)
    else:
        non_residue = least_quadratic_nonresidue(modulus)
    P = Fq['z']; (z,) = P._first_ngens(1)
    assert((z**2 -non_residue).is_irreducible())
    Fq2 = GF(modulus**2, modulus=z**2-non_residue, names=('u',)); (u,) = Fq2._first_ngens(1)
    PP = Fq2['zz']; (zz,) = PP._first_ngens(1)
    nqr = u
    while(not (zz**2 -nqr).is_irreducible()):
        nqr += 1 
    nqr_to_t = nqr**t
    Frobenius_coeff_C1 = lift(Mod(non_residue, modulus)**((modulus**0 -1 )/2 )), lift(Mod(non_residue, modulus)**((modulus**1 -1 )/2 ))
    return euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1[0 ], Frobenius_coeff_C1[1 ]

def parameters_Fp6(modulus, non_residue, coeff_b, r):
    """
    Description:

        + Computes parameters of a sextic extension finite field F_(modulus^6)
          F_(modulus^6)[v] = F_(modulus^2) / v^3 - non_residue

        + Computes generators of subgroups G1 and G2

    Input:

        mmodulus - integer

    Output:
        nqr, Frobenius_coeffs_c1, Frobenius_coeffs_c2, mul_by_q, coeff_b_twist, twist_type, G1_one, G2_one
    """ 
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
    "compute G2 generator for later"
    Et = EllipticCurve(Fq2, [0 , coeff_b_twist])
    G2_one = Et.random_point()
    return nqr, Frobenius_coeffs_c1, Frobenius_coeffs_c2, mul_by_q, coeff_b_twist, twist_type, G1_one, G2_one

def parameters_Fp12(modulus, nqr):
    """
    Description:

        Computes parameters of a dodedic':)' extension finite field F_(modulus^12)
        F_(modulus^12)[u] = F_(modulus^6) / u^2 - non_residue

    Input:

        mmodulus - integer

    Output:
        euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1[0], Frobenius_coeff_C1[1]
    """ 

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
