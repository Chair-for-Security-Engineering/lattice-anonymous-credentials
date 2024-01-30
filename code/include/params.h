#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* Domain separators for XOF expansion
**************************************************/
// sep
#define DOMAIN_SEPARATOR_A 0
#define DOMAIN_SEPARATOR_R 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3
#define DOMAIN_SEPARATOR_D 4
// osig
#define DOMAIN_SEPARATOR_DS 5
#define DOMAIN_SEPARATOR_S 6
#define DOMAIN_SEPARATOR_RAND 7
#define DOMAIN_SEPARATOR_A1_ISS 8
#define DOMAIN_SEPARATOR_A2_ISS 9
#define DOMAIN_SEPARATOR_BYG_ISS 10
#define DOMAIN_SEPARATOR_B_ISS 11
#define DOMAIN_SEPARATOR_CHAL1_ISS 12
#define DOMAIN_SEPARATOR_CHAL2_ISS 13
#define DOMAIN_SEPARATOR_CHAL3_ISS 14
#define DOMAIN_SEPARATOR_CHAL4_ISS 15
#define DOMAIN_SEPARATOR_RAND_S2_ISS 16
#define DOMAIN_SEPARATOR_RAND_G_ISS 17
// show
#define DOMAIN_SEPARATOR_A1_SHOW 18
#define DOMAIN_SEPARATOR_A2_SHOW 19
#define DOMAIN_SEPARATOR_BYG_SHOW 20
#define DOMAIN_SEPARATOR_B_SHOW 21
#define DOMAIN_SEPARATOR_CHAL1_SHOW 22
#define DOMAIN_SEPARATOR_CHAL2_SHOW 23
#define DOMAIN_SEPARATOR_CHAL3_SHOW 24
#define DOMAIN_SEPARATOR_CHAL4_SHOW 25
#define DOMAIN_SEPARATOR_RAND_S2_SHOW 26
#define DOMAIN_SEPARATOR_RAND_G_SHOW 27

/*************************************************
* Signature parameters
**************************************************/
// Ring degree for the signature
#define PARAM_N 256
// Modulus for the signature
#define PARAM_Q 425801L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 20
// Module rank for the signature
#define PARAM_D 4
// Gadget dimension
#define PARAM_K 5
// Gadget base
#define PARAM_B 14
// Dimension of the message vector (without usk)
#define PARAM_M 10
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 5
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 7390.20585207465319399489
// Gaussian parameter s_2 for v_2 and v_3
#define PARAM_S2 68.17015305110869860528
// Squared Gaussian parameter s_1^2
#define PARAM_S1SQ 34270592.82034289091825485229
// Gaussian width for p_2
#define PARAM_SQRT_S2SQ_SGSQ 48.26510947391209782609
// Negated ratio -1/(1/s_G^2 - 1/s_2^2)
#define PARAM_SGINVSQ_S2INVSQ -4623.48663266090716206236
// Negated ratio -s_G^2/(s_2^2 - s_G^2)
#define PARAM_NEGSGSQ_DIV_S2SQ_SGSQ -0.99490375098435301915
// Squared verification bound on v_1 (hiding case)
#define PARAM_B1SQ 16568582601UL
// Squared verification bound on v_2
#define PARAM_B2SQ 4886925UL
// Squared verification bound on v_3
#define PARAM_B3SQ 1544265UL

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the public seed for CRS expansion
#define CRS_SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

/*************************************************
* [ISSUANCE] Zero-Knowledge proof parameters
**************************************************/
// Ring degree for the issuance proof
#define PARAM_N_ISS 64
// Ring degree gap between the issuance proof and the signature (subring embedding)
#define PARAM_K_ISS 4
// Modulus for the issuance proof
#define PARAM_Q_ISS 223205310001L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_ISS_BITLEN 38
// Modulus factor for the issuance proof
#define PARAM_Q1_ISS 524201L
// Modulus bit-length upper bound for uniform sampling mod q_1
#define PARAM_Q1_ISS_BITLEN 20
// Second modulus factor for the issuance proof
#define PARAM_Q2_ISS PARAM_Q
// Inverse of q_1 modulo q_2
#define PARAM_Q1_INVMOD_Q2_ISS 343579L
// Inverse of q_2 modulo q_1
#define PARAM_Q2_INVMOD_Q1_ISS 101223L
// Module rank for the issuance proof
#define PARAM_D_ISS 20
// Witness dimension
#define PARAM_M1_ISS 104
// Scaled witness dimension (m_1 / k_hat)
#define PARAM_M1_K_ISS 26
// ABDLOP commitment randomness dimension
#define PARAM_M2_ISS 58
// Soundness amplification dimension
#define PARAM_L_ISS 7
// Dimension for Approximate Range Proof
#define PARAM_ARP_ISS 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_ISS 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_ISS 11
// Gaussian mask width for cs_1
#define PARAM_S1_ISS 369050.89730269293067976832
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_ISS 136198564799.92280578613281250000
// Gaussian mask width for cs_2
#define PARAM_S2_ISS 275602.77920886297943070531
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_ISS 75956891907.64927673339843750000
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_ISS 72848.10643310110026504844
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_ISS 5306846610.88842582702636718750
// Rejection sampling rate for y_1
#define PARAM_REJ1_ISS 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_ISS 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_ISS 2
// Squared verification bound for z_1
#define PARAM_B1SQ_ISS 180657566352359UL
// Squared verification bound for z_2
#define PARAM_B2SQ_ISS 60411097594469UL
// Squared verification bound for z_3
#define PARAM_B3SQ_ISS 584702794673UL
// Infinity norm of challenges
#define PARAM_RHO_ISS 8
// Manhattan-like norm of challenges
#define PARAM_ETA_ISS 93

/*************************************************
* [SHOW] Zero-Knowledge proof parameters
**************************************************/
// Ring degree for the show proof
#define PARAM_N_SHOW 64
// Ring degree gap between the issuance proof and the signature (subring embedding)
#define PARAM_K_SHOW 4
// Modulus for the issuance proof
#define PARAM_Q_SHOW 234086575306343681L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_SHOW_BITLEN 58
// Modulus factor for the show proof
#define PARAM_Q1_SHOW 549755813881L
// Modulus bit-length upper bound for uniform sampling mod q_1
#define PARAM_Q1_SHOW_BITLEN 40
// Second modulus factor for the issuance proof
#define PARAM_Q2_SHOW PARAM_Q
// Inverse of q_1 modulo q_2
#define PARAM_Q1_INVMOD_Q2_SHOW 99299L
// Inverse of q_2 modulo q_1
#define PARAM_Q2_INVMOD_Q1_SHOW 421549908863L
// Module rank for the show proof
#define PARAM_D_SHOW 23
// Witness dimension
#define PARAM_M1_SHOW 211
// ABDLOP commitment randomness dimension
#define PARAM_M2_SHOW 74
// Soundness amplification dimension
#define PARAM_L_SHOW 7
// Dimension for Approximate Range Proof
#define PARAM_ARP_SHOW 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_SHOW 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_SHOW 11
// Gaussian mask width for cs_1
#define PARAM_S1_SHOW 582380223.29294335842132568359
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_SHOW 339166724482738560.00000000000000000000
// Gaussian mask width for cs_2
#define PARAM_S2_SHOW 311304.54102290823357179761
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_SHOW 96910517261.48355102539062500000
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_SHOW 114957846.73890274763107299805
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_SHOW 13215306526845054.00000000000000000000
// Rejection sampling rate for y_1
#define PARAM_REJ1_SHOW 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_SHOW 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_SHOW 2
// Squared verification bound for z_1 (high bits)
#define PARAM_B1SQ_SHOW_LOW64 6567424658286313472UL
// Squared verification bound for z_1 (low bits)
#define PARAM_B1SQ_SHOW_HIGH64 46UL
// Squared verification bound for z_2
#define PARAM_B2SQ_SHOW 95184984511325UL
// Squared verification bound for z_3
#define PARAM_B3SQ_SHOW 1456048615171063808UL
// Infinity norm of challenges
#define PARAM_RHO_SHOW 8
// Manhattan-like norm of challenges
#define PARAM_ETA_SHOW 93

/*************************************************
* Testing
**************************************************/

// The modulus factor
#define Q_1_MOD 47
// The modular inverse of Q_1_MOD mod Q_MIN_MOD (used for CRT reconstruction)
#define Q_1_MOD_INV 2
// The smallest modulus factor
#define Q_MIN_MOD 31
// The modular inverse of Q_MIN_MOD mod Q_1 (used for CRT reconstruction)
#define Q_MIN_MOD_INV 44
// The modulus for the proof system
#define Q_HAT_MOD (Q_1_MOD * Q_MIN_MOD)
// The degree of the ring for proofs
#define N_HAT_RING 4

#endif /* PARAMS_H */
