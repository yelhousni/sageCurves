# -*- coding: utf-8 -*-

# Author: Youssef El Housni
# youssef.el.housni@fr.ey.com / youssef.housni21@gmail.com 

import pystache
from sage.all_cmdline import *   # import sage library
import sys

#TODO: for BN 'friendly' curves write coeff_b as c^4+d^6 or c^6+4d^4 and most parameters follow easily

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


# Field utils
def linearPoly2tab(poly):
    coeffs_list = poly.polynomial().list()
    c0 = coeffs_list[0]
    c1 = 0  if len(coeffs_list)==1  else coeffs_list[1]
    return c0, c1

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
        Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t
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
    return Rsq8, Rcu8, inv8, Rsq4, Rcu4, inv4, num_bits, euler, s, t, t_minus_1_over_2, multiplicative_generator, root_of_unity, nqr, nqr_to_t

def parameters_Fp2(modulus):
    """
    Description:

        Computes parameters of a quadratic extension finite field F_(modulus^2)
        F_(modulus^2)[u] = F_modulus / u^2 - non_residue

    Input:

        mmodulus - integer

    Output:
        euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1
    """ 
    euler = int((modulus**2-1)/2)
    s = twoAdicity(modulus**2-1)
    t = int((modulus**2-1)/2**s)
    t_minus_1_over_2 = int((t-1)/2)
    Fq = GF(modulus)
    if (modulus % 4 == 3): 
        non_residue = Fq(-1) # Legendre symbol is clearly -1
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
    return euler, s, t, t_minus_1_over_2, non_residue, nqr, nqr_to_t, Frobenius_coeff_C1

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
    for i in range(6):
        Frobenius_coeffs_c1.append(nqr**int((modulus**i-1)/3))
    for i in range(6):
        Frobenius_coeffs_c2.append(nqr**int((2*modulus**i-2)/3))
    "compute mul_by_q for later"
    mul_by_q.append(Frobenius_coeffs_c1[2])
    mul_by_q.append(nqr**int((modulus**i-1)/2))
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
    # to change to fit chosen generators in standard curves
    # for bls12-381, take: https://github.com/zkcrypto/bls12_381/blob/master/src/notes/design.rs
    E = EllipticCurve(Fq, [0 , coeff_b])
    cofactor = Integer(E.order() / r)
    if ((1+coeff_b).is_square()): # it is the case for BN
        G1_one = cofactor * E(1, sqrt(1+coeff_b))
    else:
        G1_one = cofactor * E.random_point() 
    "compute G2 generator for later"
    # to change to fit chosen generators in standard curves
    # for bls12-381, take: https://github.com/zkcrypto/bls12_381/blob/master/src/notes/design.rs
    cofactor_t = Integer(Et.order()/r)
    G2_one = cofactor_t * Et.random_point()
    return nqr, Frobenius_coeffs_c1, Frobenius_coeffs_c2, mul_by_q, coeff_b_twist, twist_type, G1_one, G2_one

def parameters_Fp12(modulus, nqr):
    """
    Description:

        Computes parameters of a dodedic':)' extension finite field F_(modulus^12)
        F_(modulus^12)[u] = F_(modulus^6) / u^2 - non_residue

    Input:

        mmodulus - integer

    Output:
        non_residue, frobenius_coeffs
    """ 

    Frob12_coeff_C1 = []
    for i in range(12 ):
        Frob12_coeff_C1.append(nqr**int((modulus**i-1 )/6 ))
    return nqr, Frob12_coeff_C1

def pairing_parameters(family, x, q, k, r):
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
    return ate_loop_count, ate_loop_is_neg, final_exponent, final_exponent_is_z_neg

def print_to_libff(curve_name, curve_family, modulus_r, R2_64_r, R3_64_r, inv_64_r, R2_32_r, R3_32_r, inv_32_r, num_bits_r, euler_r, s_r, t_r, t_minus_1_over_2_r, multiplicative_generator_r, root_of_unity_r, nqr_r, nqr_to_t_r, R2_64_q, R3_64_q, inv_64_q, modulus_q, R2_32_q, R3_32_q, inv_32_q, num_bits_q, euler_q, s_q, t_q, t_minus_1_over_2_q, multiplicative_generator_q, root_of_unity_q, nqr_q, nqr_to_t_q, euler_q2, s_q2, t_q2, t_minus_1_over_2_q2, non_residue_q2, nqr_q2, nqr_to_t_q2, frobenius_q2, non_residue_q6, frobenius_q6_1, frobenius_q6_2, mul_by_q, twist_coeff_b, twist_type, G1_one, G2_one, non_residue_q12, frobenius_q12, ate_loop_count, ate_loop_count_bool, final_exponent, x_bool):
    renderer = pystache.Renderer()
    print(renderer.render_path('./templates/init.mustache', 
        {'curve_name': curve_name},
        {'curve_family': curve_family},
        {'modulus_r': modulus_r},
        {'R2_64_r': R2_64_r},
        {'R3_64_r': R3_64_r},
        {'inv_64_r': inv_64_r},
        {'R2_32_r': R2_32_r},
        {'R3_32_r': R3_32_r},
        {'inv_32_r': inv_32_r},
        {'num_bits_r': num_bits_r},
        {'euler_r': euler_r},
        {'s_r': s_r},
        {'t_r': t_r},
        {'t_minus_1_over_2_r': t_minus_1_over_2_r},
        {'multiplicative_generator_r': multiplicative_generator_r},
        {'root_of_unity_r': root_of_unity_r},
        {'nqr_r': nqr_r},
        {'nqr_to_t_r':  nqr_to_t_r},
        {'modulus_q': modulus_q},
        {'R2_64_q': R2_64_q},
        {'R3_64_q': R3_64_q},
        {'inv_64_q': inv_64_q},
        {'R2_32_q': R2_32_q},
        {'R3_32_q': R3_32_q},
        {'inv_32_q': inv_32_q},
        {'num_bits_q': num_bits_q},
        {'euler_q': euler_q},
        {'s_q': s_q},
        {'t_q': t_q},
        {'t_minus_1_over_2_q': t_minus_1_over_2_q},
        {'multiplicative_generator_q': multiplicative_generator_q},
        {'root_of_unity_q': root_of_unity_q},
        {'nqr_q': nqr_q},
        {'nqr_to_t_q':  nqr_to_t_q},
        {'euler_q2': euler_q2},
        {'s_q2': s_q2},
        {'t_q2': t_q2},
        {'t_minus_1_over_2_q2': t_minus_1_over_2_q2}))
    """
        {'non_residue_q2': non_residue_q2},
        {'nqr_q2': nqr_q2},
        {'nqr_to_t_q2': nqr_to_t_q2},
        {'frobenius_q2': frobenius_q2},
        {'non_residue_q6': non_residue_q6},
        {'frobenius_q6_1': frobenius_q6_1},
        {'frobenius_q6_2': frobenius_q6_2},
        {'mul_by_q': mul_by_q},
        {'twist_coeff_b': twist_coeff_b},
        {'twist_type': twist_type},
        {'G1_one': G1_one},
        {'G2_one': G2_one},
        {'non_residue_q12': non_residue_q12},
        {'frobenius_q12': frobenius_q12},
        {'ate_loop_count': ate_loop_count},
        {'ate_loop_count_bool': ate_loop_count_bool},
        {'final_exponent': final_exponent},
        {'x_bool': x_bool}))
    """
