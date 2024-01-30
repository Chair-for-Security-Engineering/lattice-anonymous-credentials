# -*- coding: utf-8 -*-
"""
    Brief: Main parameters classes
"""
#-- Import --#
from math import sqrt, exp, log, log2, floor, ceil, pi, e
from scipy.special import comb
from scipy.optimize import root
from sympy import isprime, prevprime
from estimate_SIS_LWE import *
#-- End Import --#
#-- Global Parameters --#
COST_MODEL      = 'realistic_sieving'
QUANTUM         = False
LOG2_EPS        = -40
SPECTRAL_SLACK  = 6
#-- End Global Parameters --#

def c_star(dim:int, sec:float):
    """
    Find the Gaussian tailcut rate c so that the upper bound c*s*sqrt(dim) is
    verified with probability at least 1 - 2^(-sec)
    - input:
        (int)   dim     -- Dimension
        (int)   sec     -- Security parameter
    - output:
        (flt)   c_star  -- Tailcut rate
    """
    f = lambda c: sec + dim * (log2(c * sqrt(2*pi)) + (1/2 - pi*c**2)*log2(e))
    return root(f, 1)['x'][0]

class SEP_Parameters:
    """
    Main class containing all the parameters for the 
    signature scheme. 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m:int, q_start:int):
        """
        Computing all the signature parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m               -- Message dimension (without secret key)
            (int)   q_start         -- Starting search modulus (q largest adequately splitted prime below q_start)
        - output:
            SEP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   m_s     -- User secret dimension (2d)
                    (int)   m       -- Secret (2d) + Message dimension
                    (int)   q       -- Modulus
                    (int)   k       -- Gadget vector dimension
                    (int)   b       -- Gadget base
                    (flt)   s_G     -- Gadget sampling width
                    (flt)   s_1     -- Top preimage sampling width
                    (flt)   s_2     -- Bottom preimage sampling width
                    (int)   w       -- Hamming weight of tags
                    (int)   kappa   -- Number of modulus splitting factor
                    (int)   Q       -- Maximal number of queries
                    (flt)   r       -- Smoothing of Z
                    (flt)   r_ndk   -- Smoothing of Z^(ndk)
                    (flt)   r_ndk2  -- Smoothing of Z^(nd(2+k))
                    (flt)   alpha   -- Rejection sampling slack
                    (flt)   M       -- Rejection sampling repetition rate
                    (flt)   B_1     -- Verification bound for v_1
                    (flt)   B_1_s   -- Verification bound for v_1 squared
                    (flt)   B_1_prime -- Verification bound for v_1 - r (hiding case)
                    (flt)   B_1_prime_s -- Verification bound for v_1 - r (hiding case) squared
                    (flt)   B_2     -- Verification bound for v_2
                    (flt)   B_2_s   -- Verification bound for v_2 squared
                    (flt)   B_3     -- Verification bound for v_3
                    (flt)   B_3_s   -- Verification bound for v_3 squared
                [+] Required M-SIS security (for x = I or x = II)
                    (int)   req_msis_coresvp_x  -- Required M-SIS CoreSVP hardness
                [+] Key Recovery Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical CoreSVP hardness
                    (int)   mlwe_coresvp_q  -- Quantum CoreSVP hardness
                [+] Forgery Security (for x = I or x = II)
                    (int)   msis_beta_inf_x -- M-SIS Infinity norm bound
                    (flt)   msis_beta_x     -- M-SIS Euclidean norm bound
                    (int)   msis_subdim_x   -- Optimal M-SIS subdimension
                    (flt)   msis_rhf_x      -- Root Hermite Factor
                    (int)   msis_coresvp_c_x-- Classical CoreSVP hardness
                    (int)   msis_coresvp_q_x-- Quantum CoreSVP hardness
                [+] Sizes
                    (int)   pk_bitsize  -- Bitsize of the public key (B = AR)
                    (int)   sk_bitsize  -- Bitsize of the secret key (R)
                    (int)   sig_bitsize -- Bitsize of the signature (tag + v_{1,2} + v_2 + v_3)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R = Z[X]/<X^n + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # User secret key dimensionality (2*d)
        self.m_s = 2 * self.d

        # Message dimensionality including user secret key (2*d)
        self.m = 2 * self.d + m

        # Maximal number of signature queries (hardcoded)
        self.Q = 2 ** 32

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = 4

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q = prevprime(q_start + 1)
        while (q%(4*self.kappa) != 2*self.kappa + 1) or (q < (2*sqrt(self.kappa)) ** self.kappa):
            q = prevprime(q)
        self.q = q

        # Finding minimum Hamming weight for the tag space
        w = 1
        while comb(self.n, w) < self.Q:
            w += 1
        self.w = w

        # Gadget base for the gadget vector g = [1 b b^2 ... b^(k-1)]
        self.b = ceil(q ** (1/5))

        # Gadget dimension
        self.k = ceil(log(self.q) / log(self.b))

        # Smoothing parameters
        self.r = sqrt((log(2) - LOG2_EPS * log(2)) / pi)
        self.r_ndk = sqrt((log(2 * self.n * self.d * self.k) - LOG2_EPS * log(2)) / pi)
        self.r_ndk2 = sqrt((log(2 * self.n * self.d * (2 + self.k)) - LOG2_EPS * log(2)) / pi)

        # Temporary variables: 
        #   - bound on spectral norm of R
        #   - bound on euclidean norm of U*m
        norm_R = 7/10 * (sqrt(2 * self.d * self.n) + sqrt(self.k * self.d * self.n) + SPECTRAL_SLACK)
        self.bound_sk = norm_R
        self.square_bound_sk = norm_R ** 2
        norm_Um = sqrt(self.d * self.n) * sqrt(self.m * self.n)

        # Computing Gaussian widths for elliptic sampler
        self.s_G = self.r_ndk * sqrt(self.b ** 2 + 1)
        s_MP12 = sqrt(2 * self.s_G ** 4 / (self.s_G ** 2 - 1)) * norm_R
        self.s_1 = max( sqrt(pi / log(2)) * (norm_Um + sqrt(2 * self.d * self.n)), s_MP12)
        self.s_2 = sqrt(2 * self.s_G ** 2 + self.r_ndk2 ** 2)

        # Computing rejection sampling parameters for security proof
        self.alpha = self.s_1 / (norm_Um + sqrt(2 * self.d * self.n))
        self.M = exp(pi / self.alpha ** 2)

        ### Security

        # Computing hardness of M-LWE_{n,d,d,q,U(S_1)} for key recovery
        Xs = ND.CenteredBinomial(1)
        Xe = ND.CenteredBinomial(1)
        res = estimate_LWE(n = self.n * self.d,
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.d,
                            cost_model = COST_MODEL)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing required M-SIS security for type I+II forgeries
        tag_space_size = comb(self.n, self.w)
        # M-LWE loss
        mlwe_hardness_bound = 2 ** (- (self.mlwe_coresvp_q if QUANTUM else self.mlwe_coresvp_c))
        Ek = self.k * mlwe_hardness_bound
        Em = self.m * mlwe_hardness_bound
        # Trapdoor switching loss
        a = 2 * self.sec / (2 * self.sec - 1)
        o = 2 * self.sec
        eps = 2 ** LOG2_EPS
        delta = ((1+eps)/(1-eps))**(12*self.d*(self.n - 1)+5) * ((1+eps/(self.n*self.d*self.k))/(1-eps/(self.n*self.d*self.k)))**(2*self.n*self.d*self.k)
        loss_TS = (1 + o*(o-1)/(2*(2-delta)**(o+1)) * (delta-1)**2)**(self.Q/o)

        # Computing loss through hybrid argument
        loss = lambda Adv: ((((Adv - Ek)*1/loss_TS) ** a - 2*Ek)*1/loss_TS) ** a - Ek # Equation (2)
        advantage_I = 2 ** (-self.sec)
        for _ in range(self.d):
            advantage_I = loss(advantage_I)
        advantage_I = 1/(2*(tag_space_size - self.Q)) * advantage_I
        if advantage_I < 0:
            raise ValueError('ERROR in M-SIS required security: M-LWE hardness not sufficient for Type I')
        else:
            self.req_msis_coresvp_I = ceil( - log2(advantage_I) )
        # Computing loss through hybrid argument
        advantage_II = 1/(4*self.M) * (1-eps)/(1+eps) * 2 ** (-self.sec)
        for _ in range(self.d):
            advantage_II = loss(advantage_II)
        advantage_II = 1/(4*self.Q) * advantage_II
        if advantage_II < 0:
            raise ValueError('ERROR in M-SIS required security: M-LWE hardness not sufficient for Type II')
        else:
            self.req_msis_coresvp_II = ceil( - log2(advantage_II) )

        # Computing M-SIS security for type I+II forgeries
        # Square verification bounds
        self.B_1_s = floor(c_star(2 * self.d * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (2 * self.d * self.n))
        self.B_1 = sqrt(self.B_1_s)
        self.B_1_prime_s = floor((self.B_1 + sqrt(2 * self.d * self.n)) ** 2)
        self.B_1_prime = sqrt(self.B_1_prime_s)
        self.B_2_s = floor(c_star(self.k * self.d * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.k * self.d * self.n))
        self.B_2 = sqrt(self.B_2_s)
        self.B_3_s = floor(c_star(self.k * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.k * self.n))
        self.B_3 = sqrt(self.B_3_s)

        self.msis_beta_inf_I = self.s_1 * log2(self.sec) + self.k * self.d * self.n * self.s_2 * log2(self.sec) + 1
        self.msis_beta_I = sqrt(floor( (self.B_1_prime + sqrt(self.d * self.n) * self.B_2) ** 2 + self.B_3 ** 2 + self.m * self.n + 1 ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k + self.m + 1),
                        q = self.q,
                        beta = self.msis_beta_I,
                        betainf = self.msis_beta_inf_I,
                        cost_model = COST_MODEL)
        self.msis_subdim_I = res[0]
        self.msis_rhf_I = res[1]
        self.msis_bkz_I = res[2]
        self.msis_coresvp_c_I = res[3]
        self.msis_coresvp_q_I = res[4]

        self.msis_beta_inf_II = 2 * self.s_1 * log2(self.sec) + self.k * self.d * self.n * self.s_2 * log2(self.sec) + 1
        self.msis_beta_II = sqrt(floor( (2 * self.B_1_prime + 2 * sqrt(self.d * self.n) * self.B_2 + norm_Um) ** 2 + 4 * self.B_3 ** 2 ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k),
                        q = self.q,
                        beta = self.msis_beta_II,
                        betainf = self.msis_beta_inf_II,
                        cost_model = COST_MODEL)
        self.msis_subdim_II = res[0]
        self.msis_rhf_II = res[1]
        self.msis_bkz_II = res[2]
        self.msis_coresvp_c_II = res[3]
        self.msis_coresvp_q_II = res[4]

        ### Efficiency

        # Size of public key: |pk| = |B| (+ seed)
        self.pk_bitsize = self.d * self.k * self.d * self.n * ceil(log2(self.q)) + 256

        # Size of secret key: |sk| = |R| (perturbation sampling material not included)
        self.sk_bitsize = 2 * self.d * self.k * self.d * self.n * ceil(log2(3))

        # Size of signature: |sig| = |tag| + |v_{1,2}| + |v_2| + |v_3|
        self.tag_bitsize    = self.n
        self.v_12_bitsize   = ceil(self.n * self.d * (1/2 + log2(self.s_1)))
        self.v_2_bitsize    = ceil(self.k * self.d * self.n * (1/2 + log2(self.s_2)))
        self.v_3_bitsize    = ceil(self.k * self.n * (1/2 + log2(self.s_2)))
        
        self.sig_bitsize    = self.tag_bitsize + self.v_12_bitsize + self.v_2_bitsize + self.v_3_bitsize       
    
    def __repr__(self):
        """
        Printing a SEP_Parameters object
        """
        tmp = '\n[+] Signature Scheme Parameters\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Ring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('User secret dimension (2d)', 'm_s', self.m_s)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Overall message dimension (secret+message)', 'm', self.m)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget dimension', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget base', 'b', self.b)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gadget sampling width', 's_G', self.s_G)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Top preimage sampling width', 's_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Bottom preimage sampling width', 's_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Hamming weights of tags', 'w', self.w)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack (security proof)', 'α', self.alpha)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate (security proof)', 'M', self.M)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Maximal number of queries', 'Q', floor(log2(self.Q)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Smoothing loss', 'ε', LOG2_EPS)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z', 'η(1)', self.r)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(ndk)', 'η(ndk)', self.r_ndk)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(nd(2+k))', 'η(nd(2+k))', self.r_ndk2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 1', 'B_1', self.B_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 1 (hiding)', 'B_1\'', self.B_1_prime)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 2', 'B_2', self.B_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 3', 'B_3', self.B_3)
        tmp += '\n[+] Required M-SIS hardness for {:3d} bits of signature security\n'.format(self.sec)
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required M-SIS hardness for type I forgeries', 'λ_I', self.req_msis_coresvp_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required M-SIS hardness for type II forgeries', 'λ_II', self.req_msis_coresvp_II)
        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Type-I Forgeries (M-SIS)') + '\n'
        # tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Infinity norm bound on M-SIS solution', 'β_oo (I)', self.msis_beta_inf_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (I)', self.msis_beta_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (I)', self.msis_subdim_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (I)', self.msis_rhf_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (I)', self.msis_bkz_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec (I)', self.msis_coresvp_c_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec (I)', self.msis_coresvp_q_I)
        tmp += '{:-^100s}'.format('Type-II Forgeries (M-SIS)') + '\n'
        # tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Infinity norm bound on M-SIS solution', 'β_oo (II)', self.msis_beta_inf_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (II)', self.msis_beta_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (II)', self.msis_subdim_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (II)', self.msis_rhf_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (II)', self.msis_bkz_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec (II)', self.msis_coresvp_c_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec (II)', self.msis_coresvp_q_II)
        tmp += '{:-^100s}'.format('Key Recovery (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_coresvp_q)
        tmp += '\n[+] Signature Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Public key size (KB)', '|pk|', self.pk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Secret key size (KB)', '|sk|', self.sk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Tag size (B)', '|tag|', self.tag_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_{1,2} size (B)', '|v_{1,2}|', self.v_12_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_2 size (B)', '|v_2|', self.v_2_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_3 size (B)', '|v_3|', self.v_3_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Overall Signature size (B)', '|sig|', self.sig_bitsize / 2 ** 3.)
        return tmp


class Issue_ZKP_Parameters:
    """
    Main class containing all the parameters for the 
    zero-knowledge proof system (Show protocol). 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_1_start:int, n_attr:int, sig):
        """
        Computing all the zero-knowledge argument parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m_2             -- Commitment randomness dimension
            (int)   q_1_start       -- Starting search modulus (q_1 largest adequately splitted prime below q_1_start)
            (int)   n_attr          -- Number of disclosed attributes
            (SEP_Parameters) sig    -- Signature parameters
        - output:
            Show_ZKP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q_1     -- Modulus factor
                    (int)   q_2     -- Second modulus factor (signature modulus)
                    (int)   q_min   -- Smallest modulus factor
                    (int)   q       -- Modulus
                    (int)   l       -- Parameter for soundness amplification
                    (int)   n_attr  -- Number of disclosed attributes
                    (int)   m_1     -- Dimension of witness
                    (int)   m_2     -- Dimension of ABDLOP commitment randomness
                    (int)   xi_s2   -- Infinity norm of ABDLOP commitment randomness
                    (int)   rho     -- Infinity norm of challenges
                    (int)   eta     -- Manhattan-like norm of challenges
                    (int)   challenge_space_size -- Size of challenge space
                    For x = 1, 2, 3
                        (flt)   M_x         -- Rejection sampling repetition rate for z_x
                        (flt)   log_eps_x   -- Rejection sampling (log) loss for z_x
                        (flt)   alpha_x     -- Rejection sampling slack for z_x
                        (flt)   s_x         -- Gaussian width for y_x
                    (flt)   bound_witness -- Euclidean norm bound on the witness
                    (flt)   B_3     -- Euclidean norm bound for approximate range proofs
                [+] Zero-Knowledge Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical CoreSVP hardness
                    (int)   mlwe_coresvp_q  -- Quantum CoreSVP hardness
                [+] Soundness Security
                    (flt)   msis_beta       -- M-SIS Euclidean norm bound
                    (int)   msis_subdim     -- Optimal M-SIS subdimension
                    (flt)   msis_rhf        -- Root Hermite Factor
                    (int)   msis_coresvp_c  -- Classical CoreSVP hardness
                    (int)   msis_coresvp_q  -- Quantum CoreSVP hardness
                    (flt)   soundness_error -- Soundness error
                [+] Sizes
                    (int)   incompressible_bitsize  -- Bitsize of incompressible elements (in R_q)
                    (int)   compressible_bitsize    -- Bitsize of compressible elements (Gaussians)
                    (int)   challenge_bitsize       -- Bitsize of challenge (c)
                    (int)   proof_bitsize           -- Bitsize of proof (π)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R'' = Z[X]/<X^n' + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Number of disclosed attributes
        self.n_attr = n_attr

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = sig.kappa

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q_1 = prevprime(q_1_start + 1)
        while (q_1%(4*self.kappa) != 2*self.kappa + 1) or (q_1 < (2*sqrt(self.kappa)) ** self.kappa):
            q_1 = prevprime(q_1)
        self.q_1 = q_1
        self.q_2 = sig.q
        self.q = self.q_1 * sig.q
        self.q_min = min(self.q_1, sig.q)

        # Repetition for soundness amplification
        self.l = ceil(self.sec / log2(self.q_min))

        # Infinity norm bound on the challenges
        self.rho = ceil(1/2 * (2 ** (2*(self.sec + 1)/self.n) - 1))

        # Manhattan-like norm bound on the challenges (hardcoded)
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]

        # Size of challenge space
        self.challenge_space_size = (2 * self.rho + 1) ** (self.n // 2) / 2

        # Subring gap
        self.k = sig.n // self.n

        # Witness dimension
        self.m_1 = self.k * (2 * sig.d + sig.m - self.n_attr) # m includes secret and message

        # Commitment randomness dimension and infinity norm bound (hardcoded)
        self.m_2 = m_2
        self.xi_s2 = 1

        # Bound on Euclidean norm of the witness
        self.bound_witness = sqrt(self.n * self.m_1)

        # Rejection sampling parameters (hardcoded)
        self.M_1 = 2
        self.M_2 = 2
        self.M_3 = 2
        self.log_eps_1 = -130
        self.log_eps_2 = -130
        self.log_eps_3 = -130
        self.alpha_1 = sqrt(pi)/log(self.M_1) * (sqrt(-self.log_eps_1 * log(2) + log(self.M_1)) + sqrt(-self.log_eps_1 * log(2)))
        self.alpha_2 = sqrt(pi)/log(self.M_2) * (sqrt(-self.log_eps_2 * log(2) + log(self.M_2)) + sqrt(-self.log_eps_2 * log(2)))
        self.alpha_3 = sqrt(pi)/log(self.M_2) * (sqrt(-self.log_eps_3 * log(2) + log(self.M_3)) + sqrt(-self.log_eps_3 * log(2)))

        # Gaussian widths
        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * self.bound_witness

        # Checking approximate range proofs bounds
        self.B_256_s = floor(c_star(256, self.sec + 3) ** 2 * self.s_3 ** 2 * 256)
        self.B_256 = sqrt(self.B_256_s)
        assert self.q > self.bound_witness ** 2
        assert self.q > 2 * 41/sqrt(26) * self.n * self.m_1 * self.B_256
        assert self.q > 2 * self.B_256_s / 13 - self.B_256

        # Square Verification bounds
        self.B_1_s = floor(c_star(self.m_1 * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (self.m_1 * self.n))
        self.B_2_s = floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))

        ### Security

        # Computing hardness of M-LWE for zero-knowledge
        # Xs = ND.Uniform(-self.xi_s2,self.xi_s2)
        # Xe = ND.Uniform(-self.xi_s2,self.xi_s2)
        Xs = ND.CenteredBinomial(self.xi_s2)
        Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n = self.n * (self.m_2 - (self.d + floor(256/self.n) + self.l + 1)),
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.m_2,
                            cost_model = COST_MODEL)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for soundness
        self.msis_beta = 8 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (self.m_1 + self.m_2),
                        q = self.q,
                        beta = self.msis_beta,
                        cost_model = COST_MODEL)
        self.msis_subdim = res[0]
        self.msis_rhf = res[1]
        self.msis_bkz = res[2]
        self.msis_coresvp_c = res[3]
        self.msis_coresvp_q = res[4]

        self.soundness_error = self.q_min ** (-self.l) + self.q_min ** (-self.n/self.kappa) + 2 / self.challenge_space_size + \
                                    2 ** (-(self.msis_coresvp_q if QUANTUM else self.msis_coresvp_c))

        ### Efficiency

        # CRS size
        ajtai_crs_size = self.d * (self.m_1 + self.m_2) * (self.n * ceil(log2(self.q)))
        bdlop_crs_size = (256 / self.n + self.l + 1) * self.m_2 * (self.n * ceil(log2(self.q)))
        self.crs_size  = ajtai_crs_size + bdlop_crs_size

        # Size of incompressible elements (those uniform in R_q)
        self.incompressible_bitsize = self.n * (self.d + 256/self.n + 2 * self.l + 1) * ceil(log2(self.q))

        # Size of compressible elements (Gaussians): 
        self.compressible_bitsize = ceil(self.n * self.m_1 * (1/2 + log2(self.s_1))) + \
                                    ceil(self.n * self.m_2 * (1/2 + log2(self.s_2))) + \
                                    ceil(256 * (1/2 + log2(self.s_3)))

        # Size of challenge
        self.challenge_bitsize = self.n * ceil(log2(2 * self.rho + 1))

        # Size of proof
        self.proof_bitsize = self.incompressible_bitsize + self.compressible_bitsize + self.challenge_bitsize
    
    def __repr__(self):
        """
        Printing a Show_ZKP_Parameters object
        """
        tmp = '\n[+] Zero-Knowledge Proof Parameters\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring gap factor (n_signature / n_proof)', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus factor', 'q_1', self.q_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Smallest modulus factor', 'q_min', self.q_min)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Parameter for soundness amplification', 'l', self.l)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Number of disclosed attributes', '|I|', self.n_attr)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of witness', 'm_1', self.m_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of ABDLOP commitment randomness', 'm_2', self.m_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of commitment randomness', 'ξ', self.xi_s2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of challenges', 'ρ', self.rho)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Manhattan-like norm of challenges', 'η', self.eta)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Size of challenge space', '|C|', floor(log2(self.challenge_space_size)))        
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 1', 'M_1', self.M_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 2', 'M_2', self.M_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3', 'M_3', self.M_3)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 1', 'ε_1', self.log_eps_1)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 2', 'ε_2', self.log_eps_2)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3', 'ε_3', self.log_eps_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 1 (security proof)', 'α_1', self.alpha_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 2 (security proof)', 'α_2', self.alpha_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (security proof)', 'α_3', self.alpha_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_1', 'σ_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_2', 'σ_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_3', 'σ_3', self.s_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness', 'B', self.bound_witness)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Norm bound for approximate range proof', 'B_256', self.B_256)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 1: B^2', 'Cond. 1', ceil(self.bound_witness ** 2))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 2: 16*n*m_1*B_256', 'Cond. 2', ceil(2 * 41/sqrt(26) * self.n * self.m_1 * self.B_256))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 3: 2*B_256^2 / 13 - B_256', 'Cond. 3', ceil(2 * self.B_256 ** 2 / 13 - self.B_256))
        tmp += '| {:60s} | {:^10s} | 2^{:<18.5f} |\n'.format('Soundness error', 'δ_s', log2(self.soundness_error))

        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Soundness (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β', self.msis_beta)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim',self.msis_subdim)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.msis_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.msis_bkz)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec', self.msis_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.msis_coresvp_q)
        tmp += '{:-^100s}'.format('Zero-Knowledge (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_coresvp_q)
        tmp += '\n[+] Proof Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Common Random String (KB)', '|crs|', self.crs_size / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Incompressible elements (KB)', '|π_1|', self.incompressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Compressible elements (KB)', '|π_2|', self.compressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Challenge (KB)', '|π_3|', self.challenge_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Total proof bitsize (KB)', '|π|', self.proof_bitsize / 2 ** 13.)
        return tmp

class Show_ZKP_Parameters:
    """
    Main class containing all the parameters for the 
    zero-knowledge proof system (Show protocol). 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_1_start:int, n_attr:int, sig):
        """
        Computing all the zero-knowledge argument parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m_2             -- Commitment randomness dimension
            (int)   q_1_start       -- Starting search modulus (q_1 largest adequately splitted prime below q_1_start)
            (int)   n_attr          -- Number of disclosed attributes
            (SEP_Parameters) sig    -- Signature parameters
        - output:
            Show_ZKP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q_1     -- Modulus factor
                    (int)   q_2     -- Second modulus factor (signature modulus)
                    (int)   q_min   -- Smallest modulus factor
                    (int)   q       -- Modulus
                    (int)   l       -- Parameter for soundness amplification
                    (int)   n_attr  -- Number of disclosed attributes
                    (int)   m_1     -- Dimension of witness
                    (int)   m_2     -- Dimension of ABDLOP commitment randomness
                    (int)   xi_s2   -- Infinity norm of ABDLOP commitment randomness
                    (int)   rho     -- Infinity norm of challenges
                    (int)   eta     -- Manhattan-like norm of challenges
                    (int)   challenge_space_size -- Size of challenge space
                    For x = 1, 2, 3
                        (flt)   M_x         -- Rejection sampling repetition rate for z_x
                        (flt)   log_eps_x   -- Rejection sampling (log) loss for z_x
                        (flt)   alpha_x     -- Rejection sampling slack for z_x
                        (flt)   s_x         -- Gaussian width for y_x
                    (flt)   bound_witness -- Euclidean norm bound on the witness
                    (flt)   B_3     -- Euclidean norm bound for approximate range proofs
                [+] Zero-Knowledge Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical CoreSVP hardness
                    (int)   mlwe_coresvp_q  -- Quantum CoreSVP hardness
                [+] Soundness Security
                    (flt)   msis_beta       -- M-SIS Euclidean norm bound
                    (int)   msis_subdim     -- Optimal M-SIS subdimension
                    (flt)   msis_rhf        -- Root Hermite Factor
                    (int)   msis_coresvp_c  -- Classical CoreSVP hardness
                    (int)   msis_coresvp_q  -- Quantum CoreSVP hardness
                    (flt)   soundness_error -- Soundness error
                [+] Sizes
                    (int)   incompressible_bitsize  -- Bitsize of incompressible elements (in R_q)
                    (int)   compressible_bitsize    -- Bitsize of compressible elements (Gaussians)
                    (int)   challenge_bitsize       -- Bitsize of challenge (c)
                    (int)   proof_bitsize           -- Bitsize of proof (π)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R'' = Z[X]/<X^n' + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Number of disclosed attributes
        self.n_attr = n_attr

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = sig.kappa

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q_1 = prevprime(q_1_start + 1)
        while (q_1%(4*self.kappa) != 2*self.kappa + 1) or (q_1 < (2*sqrt(self.kappa)) ** self.kappa):
            q_1 = prevprime(q_1)
        self.q_1 = q_1
        self.q_2 = sig.q
        self.q = self.q_1 * sig.q
        self.q_min = min(self.q_1, sig.q)

        # Repetition for soundness amplification
        self.l = ceil(self.sec / log2(self.q_min))

        # Infinity norm bound on the challenges
        self.rho = ceil(1/2 * (2 ** (2*(self.sec + 1)/self.n) - 1))

        # Manhattan-like norm bound on the challenges (hardcoded)
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]

        # Size of challenge space
        self.challenge_space_size = (2 * self.rho + 1) ** (self.n // 2) / 2

        # Subring gap
        self.k = sig.n // self.n

        # Witness dimension
        self.m_1 = (2 * sig.d * self.k + 1) + (sig.k * sig.d * self.k + 1) + (sig.k * self.k + 1) + self.k + (sig.m - self.n_attr) * self.k 

        # Commitment randomness dimension and infinity norm bound (hardcoded)
        self.m_2 = m_2
        self.xi_s2 = 1

        # Bound on Euclidean norm of the witness
        self.bound_witness = sqrt(sig.B_1_prime_s + sig.B_2_s + sig.B_3_s + sig.w + sig.n * (sig.m - self.n_attr))

        # Rejection sampling parameters (hardcoded)
        self.M_1 = 2
        self.M_2 = 2
        self.M_3 = 2
        self.log_eps_1 = -130
        self.log_eps_2 = -130
        self.log_eps_3 = -130
        self.alpha_1 = sqrt(pi)/log(self.M_1) * (sqrt(-self.log_eps_1 * log(2) + log(self.M_1)) + sqrt(-self.log_eps_1 * log(2)))
        self.alpha_2 = sqrt(pi)/log(self.M_2) * (sqrt(-self.log_eps_2 * log(2) + log(self.M_2)) + sqrt(-self.log_eps_2 * log(2)))
        self.alpha_3 = sqrt(pi)/log(self.M_2) * (sqrt(-self.log_eps_3 * log(2) + log(self.M_3)) + sqrt(-self.log_eps_3 * log(2)))

        # Gaussian widths
        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * self.bound_witness

        # Checking approximate range proofs bounds
        self.B_256_s = floor(c_star(256, self.sec + 3) ** 2 * self.s_3 ** 2 * 256)
        self.B_256 = sqrt(self.B_256_s)
        assert self.q > self.bound_witness ** 2
        assert self.q > 2 * 41/sqrt(26) * self.n * self.m_1 * self.B_256
        assert self.q > 2 * self.B_256_s / 13 - self.B_256

        # Square Verification bounds
        self.B_1_s = floor(c_star(self.m_1 * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (self.m_1 * self.n))
        self.B_2_s = floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))

        ### Security

        # Computing hardness of M-LWE for zero-knowledge
        # Xs = ND.Uniform(-self.xi_s2,self.xi_s2)
        # Xe = ND.Uniform(-self.xi_s2,self.xi_s2)
        Xs = ND.CenteredBinomial(self.xi_s2)
        Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n = self.n * (self.m_2 - (self.d + floor(256/self.n) + self.l + 1)),
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.m_2,
                            cost_model = COST_MODEL)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for soundness
        self.msis_beta = 8 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (self.m_1 + self.m_2),
                        q = self.q,
                        beta = self.msis_beta,
                        cost_model = COST_MODEL)
        self.msis_subdim = res[0]
        self.msis_rhf = res[1]
        self.msis_bkz = res[2]
        self.msis_coresvp_c = res[3]
        self.msis_coresvp_q = res[4]

        self.soundness_error = self.q_min ** (-self.l) + self.q_min ** (-self.n/self.kappa) + 2 / self.challenge_space_size + \
                                    2 ** (-(self.msis_coresvp_q if QUANTUM else self.msis_coresvp_c))

        ### Efficiency

        # CRS size
        ajtai_crs_size = self.d * (self.m_1 + self.m_2) * (self.n * ceil(log2(self.q)))
        bdlop_crs_size = (256 / self.n + self.l + 1) * self.m_2 * (self.n * ceil(log2(self.q)))
        self.crs_size  = ajtai_crs_size + bdlop_crs_size

        # Size of incompressible elements (those uniform in R_q)
        self.incompressible_bitsize = self.n * (self.d + 256/self.n + 2 * self.l + 1) * ceil(log2(self.q))

        # Size of compressible elements (those Gaussians): 
        self.compressible_bitsize = ceil(self.n * self.m_1 * (1/2 + log2(self.s_1))) + \
                                    ceil(self.n * self.m_2 * (1/2 + log2(self.s_2))) + \
                                    ceil(256 * (1/2 + log2(self.s_3)))
        # Size of challenge
        self.challenge_bitsize = self.n * ceil(log2(2 * self.rho + 1))

        # Size of proof
        self.proof_bitsize = self.incompressible_bitsize + self.compressible_bitsize + self.challenge_bitsize
    
    def __repr__(self):
        """
        Printing a Show_ZKP_Parameters object
        """
        tmp = '\n[+] Zero-Knowledge Proof Parameters\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring gap factor (n_signature / n_proof)', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus factor', 'q_1', self.q_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Smallest modulus factor', 'q_min', self.q_min)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Parameter for soundness amplification', 'l', self.l)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Number of disclosed attributes', '|I|', self.n_attr)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of witness', 'm_1', self.m_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of ABDLOP commitment randomness', 'm_2', self.m_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of commitment randomness', 'ξ', self.xi_s2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of challenges', 'ρ', self.rho)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Manhattan-like norm of challenges', 'η', self.eta)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Size of challenge space', '|C|', floor(log2(self.challenge_space_size)))        
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 1', 'M_1', self.M_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 2', 'M_2', self.M_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3', 'M_3', self.M_3)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 1', 'ε_1', self.log_eps_1)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 2', 'ε_2', self.log_eps_2)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3', 'ε_3', self.log_eps_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 1 (security proof)', 'α_1', self.alpha_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 2 (security proof)', 'α_2', self.alpha_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (security proof)', 'α_3', self.alpha_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_1', 'σ_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_2', 'σ_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_3', 'σ_3', self.s_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness', 'B', self.bound_witness)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Norm bound for approximate range proof', 'B_256', self.B_256)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 1: B^2', 'Cond. 1', ceil(self.bound_witness ** 2))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 2: 16*n*m_1*B_256', 'Cond. 2', ceil(2 * 41/sqrt(26) * self.n * self.m_1 * self.B_256))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 3: 2*B_256^2 / 13 - B_256', 'Cond. 3', ceil(2 * self.B_256 ** 2 / 13 - self.B_256))
        tmp += '| {:60s} | {:^10s} | 2^{:<18.5f} |\n'.format('Soundness error', 'δ_s', log2(self.soundness_error))

        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Soundness (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β', self.msis_beta)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim',self.msis_subdim)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.msis_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.msis_bkz)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec', self.msis_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.msis_coresvp_q)
        tmp += '{:-^100s}'.format('Zero-Knowledge (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_coresvp_q)
        tmp += '\n[+] Proof Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Common Random String (KB)', '|crs|', self.crs_size / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Incompressible elements (KB)', '|π_1|', self.incompressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Compressible elements (KB)', '|π_2|', self.compressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Challenge (KB)', '|π_3|', self.challenge_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Total proof bitsize (KB)', '|π|', self.proof_bitsize / 2 ** 13.)
        return tmp

sig_pms = SEP_Parameters(target_bitsec=128, n=256, d=4, m=10, q_start=ceil(2 ** 18.7))
issue_zkp_pms = Issue_ZKP_Parameters(target_bitsec=128, n=64, d=20, m_2=58, q_1_start=ceil(2 ** 19), n_attr=0, sig=sig_pms)
show_zkp_pms = Show_ZKP_Parameters(target_bitsec=128, n=64, d=23, m_2=74, q_1_start=ceil(2 ** 39), n_attr=0, sig=sig_pms)
print(sig_pms)
print(issue_zkp_pms)
print(show_zkp_pms)
