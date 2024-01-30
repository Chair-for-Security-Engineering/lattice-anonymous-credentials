# -*- coding: utf-8 -*-
"""
    Brief: Functions and estimators for SIS and LWE hardness
"""

#-- Import --#
from math import sqrt, log2, floor, ceil, pi, e
import sys, os
estimator_path = os.path.realpath('') + '/lattice-estimator'
if estimator_path not in sys.path:
    sys.path.insert(0, estimator_path)
from estimator import *
#-- End Import --#

def rhf_from_bkz(b: int):
    """
    [Chen13] heuristics to determine the Root Hermite Factor from the
    BKZ blocksize
    - input: 
        (int)   b       -- BKZ block size
    - output: 
        (flt)   delta_0 -- Corresponding Root Hermite Factor
    """
    delta_0 = float( (b/(2 * pi * e) * (pi * b)**(1/b))**(1/(2 * (b-1))) )
    return delta_0

def find_bkz_blocksize(delta_0: float):
    """
    Finds the target BKZ blocksize for a given Root Hermite Factor
    - input: 
        (flt)   delta_0 -- Root Hermite Factor
    - output: 
        (int)   b       -- Corresponding BKZ block size
    """
    b = 36
    while rhf_from_bkz(b + 100) > delta_0:
        b += 100
    while rhf_from_bkz(b + 10) > delta_0:
        b += 10
    while rhf_from_bkz(b + 2) > delta_0:
        b += 2
    while rhf_from_bkz(b + 1) > delta_0:
        b += 1
    return b + 1

def core_svp_cost(b: int, it=1, cost_model='sieving'):
    """
    Returns the classical and quantum security estimates for `it` calls
    to an SVP oracle in the corresponding cost model based on the BKZ 
    blocksize. 
    - input: 
        (int)   b           -- BKZ block size
        (int)   it          -- Number of calls to SVP oracle [Optional]
        (str)   cost_model  -- SVP oracle modeling [Optional]
    - output: 
        (int)   csec        -- Classical bit security
        (int)   qsec        -- Quantum bit security
    """
    if cost_model == 'sieving':
        csec, qsec = floor(0.292 * b + log2(it)), floor(0.257 * b + log2(it))
    elif cost_model == 'realistic_sieving':
        csec, qsec = floor(0.292 * b + log2(it) + 16.4), floor(0.257 * b + log2(it) + 16.4)
    elif cost_model == 'pessimistic':
        csec, qsec = floor(0.2075 * b + log2(it)), floor(0.2075 * b + log2(it))
    elif cost_model == 'minspace':
        csec, qsec = floor(0.368 * b + log2(it)), floor(0.2975 * b + log2(it))
    elif cost_model == 'enumeration_APS':
        exponent = 0.187 * b * log2(b) - 1.019 * b + 16.1
        csec, qsec = floor(exponent + log2(it)), floor(exponent / 2 + log2(it))
    elif cost_model == 'enumeration_HPS':
        exponent = 0.000784 * b**2 + 0.366 * b - 0.9
        csec, qsec = floor(exponent + log2(it)), floor(exponent / 2 + log2(it))
    else:
        raise NotImplementedError('Cost Model not available. Select among: \'sieving\', \'realistic_sieving\', \'pessimistic\', \'minspace\', \'enumeration_APS\', \'enumeration_HPS\'')
    return csec, qsec

def estimate_SIS(n: int, m: int, q: int, beta: float, betainf=-1, cost_model='sieving', it=1):
    """
    Estimate the Core-SVP hardness of SIS with dimension n, samples m, modulus q, L2 bound beta,
    Linf bound betainf. That is given A in Z_q^{n x m}, find x in Z_q^m such that Ax = 0 mod q 
    verifying the norm bounds.
    - input: 
        (int)   n           -- Number of rows in A
        (int)   m           -- Number of columns in A
        (int)   q           -- Modulus
        (float) beta        -- Euclidean norm bound of the solution
        (float) betainf     -- Infinity norm bound of the solution
        (str)   cost_model  -- SVP oracle modeling [Optional]
        (int)   it          -- Number of calls to SVP oracle [Optional]
    - output: 
        (int)   opt_subdim  -- Optimal Number of columns of A to keep
        (float) delta_0     -- Root Hermite Factor
        (int)   b           -- BKZ block size
        (int)   csec        -- Classical bit security
        (int)   qsec        -- Quantum bit security
    """
    if betainf != -1 and betainf > beta:
        betainf = beta
    # Trivial solution
    if betainf > q:
        raise ValueError('[WARNING] Infinity norm bound larger than modulus: trivial solution')

    # Finding optimal subdimension to attack SIS (L2 norm)
    opt_subdim = round(2 * n * log2(q) / log2(beta))
    if opt_subdim < n:
        opt_subdim = n 
    elif opt_subdim > m:
        opt_subdim = m 

    # Finding associated Root Hermite Factor
    delta_0 = float( (beta * q**(-n/opt_subdim))**(1/opt_subdim) )
    # Finding corresponding BKZ blocksize
    if delta_0 < 1.00057456569590:
        print('[!] Root Hermite Factor smaller than 1.000574: SIS security over 1024 bits in pessimistic Core-SVP model')
        b = 4935
        csec, qsec = core_svp_cost(b, it=1, cost_model='pessimistic')
    else:
        b = find_bkz_blocksize(delta_0)
        csec, qsec = core_svp_cost(b, it=it, cost_model=cost_model)
    return opt_subdim, delta_0, b, csec, qsec

def estimate_ISIS(n: int, m: int, q: int, beta: float, cost_model='sieving', it=1):
    """
    Estimate the Core-SVP hardness of ISIS with dimension n, samples m, modulus q, L2 bound beta. 
    That is given A in Z_q^{n x m} and u in Z_q^n, find x in Z_q^m such that Ax = u mod q 
    verifying the norm bound.
    - input: 
        (int)   n           -- Number of rows in A
        (int)   m           -- Number of columns in A
        (int)   q           -- Modulus
        (float) beta        -- Euclidean norm bound of the solution
        (str)   cost_model  -- SVP oracle modeling [Optional]
        (int)   it          -- Number of calls to SVP oracle [Optional]
    - output: 
        (int)   remove_k    -- 
        (float) delta_0     -- Root Hermite Factor
        (int)   b           -- BKZ block size
        (int)   csec        -- Classical bit security
        (int)   qsec        -- Quantum bit security
    """
    if beta < q:
        subdim, delta_0, b, csec, qsec = estimate_SIS(n=n, m=m, q=q, beta=beta, cost_model=cost_model, it=it)
        return m - subdim, delta_0, b, csec, qsec
    elif beta > q * sqrt(n / 12):
        raise ValueError('[WARNING] Euclidean bound too large: trivial solution')
    
    # Finding optimal number of rows to discard
    remove_k = round(m - 2 * n * log2(q) / log2(beta))
    if remove_k > m - n:
        remove_k = m - n 

    # Finding associated Root Hermite Factor
    delta_0 = float( beta ** (1/(m - remove_k)) / q ** (n / (m - remove_k)**2) ) 
    # Finding corresponding BKZ blocksize
    if delta_0 < 1.00057456569590:
        print('[!] Root Hermite Factor smaller than 1.000574: ISIS security over 1024 bits in pessimistic Core-SVP model')
        b = 4935
        csec, qsec = core_svp_cost(b, it=1, cost_model='pessimistic')
    else:
        b = find_bkz_blocksize(delta_0)
        csec, qsec = core_svp_cost(b, it=it, cost_model=cost_model)
    return remove_k, delta_0, b, csec, qsec

def estimate_LWE(n: int, q: int, Xs, Xe, m: int, it=1, cost_model='sieving'):
    """
    Estimate the Core-SVP hardness of LWE with dimension n, samples m, modulus q, secret 
    distribution Xs, error distribution Xe. That is given A in Z_q^{m x n} and b in Z_q^m,
    find (s,e) drawn from Xs^n and Xe^m such that b = As + e mod q.
    - input: 
        (int)   n           -- LWE dimension
        (int)   q           -- Modulus
        (ND)    Xs          -- Secret distribution
        (ND)    Xe          -- Error distribution
        (int)   m           -- Number of LWE samples
        (int)   it          -- Number of calls to SVP oracle [Optional]
        (str)   cost_model  -- SVP oracle modeling [Optional]
    - output: 
        (int)   b           -- BKZ block size
        (float) delta_0     -- Root Hermite Factor
        (int)   csec        -- Classical bit security
        (int)   qsec        -- Quantum bit security
    """
    # Running estimator
    results = LWE.estimate.rough(LWE.Parameters(n=n, q=q, Xs=Xs, Xe=Xe, m=m))

    b = min([results[alg]['beta'] for alg in results.keys()])
    delta_0 = rhf_from_bkz(b)
    csec, qsec = core_svp_cost(b, it=it, cost_model=cost_model)
    return b, delta_0, csec, qsec