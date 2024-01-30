#include "arith.h"
#include "randombytes.h"
#include "poly_qiss_sampling.h"
#include "sep.h"
#include "osig.h"
#include "macros.h"
#include "fips202.h"

#include <math.h>

/*************************************************
* Name:        _reject_exp [static]
*
* Description: Rejection sampling given a probability threshold. Accepts
*              with probability threshold
*
* Arguments:   - double threshold: acceptance probability threshold
* 
* Returns 1 if the sample should be rejected (ie if udbl > threshold), 0 otherwise
**************************************************/
static int _reject_exp(double threshold) {
  uint64_t u = 0;
  double udbl = 0;
  randombytes((uint8_t*)&u, 6); // 48 random bits
  udbl = (double)u / (double)(1ul << 48);
  if (udbl > threshold){
    return 1;
  }
  return 0;
}

/*************************************************
* Name:        osig_user_prove_round1 [static]
*
* Description: Compute round 1 of zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
*              - poly_qiss_vec_k *chal_1: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL1_ISS_INPUT_BYTES bytes)
*              - poly_qiss_vec_m2 s2: polynomial vectors to host ABDLOP commitment randomness (initialized)
*              - poly_qiss_vec_k *y1: array of polynomial vectors to host mask for c.s1 (initialized)
*              - poly_qiss_vec_m2 y2: polynomial vectors to host mask for c.s2 (initialized)
*              - poly_qiss_vec_256_l y3_g: polynomial vectors to host mask for R.s1 and automorphism eq. (initialized)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_mat_d_k *A1: array of polynomial matrices for A1 (CRS)
*              - const poly_qiss_mat_d_m2 A2: polynomial matrix for A2 (CRS)
*              - const poly_qiss_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const uint8_t *randomness_seed: seed to extract randomness from (allocated SEED_BYTES bytes)
*              - const uint32_t kappa: XOF domain separator due to rejections 
**************************************************/
static void osig_user_prove_round1(
    osig_proof_t                *proof,
    poly_qiss_vec_k             chal_1[PARAM_ARP_ISS][PARAM_M1_K_ISS],
    uint8_t                     buf[CHAL1_ISS_INPUT_BYTES], 
    poly_qiss_vec_m2            s2, 
    poly_qiss_vec_k             y1[PARAM_M1_K_ISS], 
    poly_qiss_vec_m2            y2, 
    poly_qiss_vec_256_l         y3_g,
    const poly_qiss_vec_k       s1[PARAM_M1_K_ISS],
    const poly_qiss_mat_d_k     A1[PARAM_M1_K_ISS],
    const poly_qiss_mat_d_m2    A2,
    const poly_qiss_mat_256l_m2 Byg,
    const uint8_t               randomness_seed[SEED_BYTES],
    const uint32_t              kappa) {
  size_t i, j;
  uint32_t kpp = kappa;
  coeff_qiss tmp_coeff;
  poly_qiss_vec_d w, tmp_vec_d;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors
  poly_qiss_vec_d_init(w);
  poly_qiss_vec_d_init(tmp_vec_d);

  // sampling ABDLOP commitment randomness
  poly_qiss_vec_m2_binomial(s2, randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_S2_ISS);

  // sampling Gaussian masks for c.s1, c.s2, R.s1
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_sample_gaussian_s1(y1[i]);
  }
  poly_qiss_vec_m2_sample_gaussian_s2(y2);
  for (i = 0; i < PARAM_ARP_DIV_N_ISS; i++) {
    for (j = 0; j < PARAM_N_ISS; j++) {
      tmp_coeff = SampleZ(0, PARAM_S3_ISS);
      poly_qiss_set_coeff(y3_g->entries[i], j, tmp_coeff);
    }
  }

  // sampling uniform mask for equations with automorphisms
  for (i = PARAM_ARP_DIV_N_ISS; i < PARAM_ARP_DIV_N_L_ISS; i++) {
    poly_qiss_uniform_but_zero(y3_g->entries[i], randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_G_ISS);
  }

  // computing commitments tA = A1.s1 + A2.s2, w = A1.y1 + A2.y2, tB = Byg.s2 + y3_g
  poly_qiss_mat_d_k_mul_vec_k(proof->tA, A1[0], s1[0]);
  poly_qiss_mat_d_k_mul_vec_k(w, A1[0], y1[0]);
  for (i = 1; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_mul_vec_k(tmp_vec_d, A1[i], s1[i]); // A1.s1 (for tA)
    poly_qiss_vec_d_add(proof->tA, proof->tA, tmp_vec_d);
    poly_qiss_mat_d_k_mul_vec_k(tmp_vec_d, A1[i], y1[i]); // A1.y1 (for w)
    poly_qiss_vec_d_add(w, w, tmp_vec_d);
  }
  poly_qiss_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, s2); // A2.s2 (for tA)
  poly_qiss_vec_d_add(proof->tA, proof->tA, tmp_vec_d);
  poly_qiss_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, y2); // A2.y2 (for w)
  poly_qiss_vec_d_add(w, w, tmp_vec_d);
  poly_qiss_mat_256l_m2_mul_vec_m2(proof->tB, Byg, s2); // Byg.s2 (for tB)
  poly_qiss_vec_256_l_add(proof->tB, proof->tB, y3_g);

  // computing first challenge
  buf[0] = 1;
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES, proof->tA);
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES + POLYQISS_VECD_PACKEDBYTES, w);
  poly_qiss_vec_256_l_pack(buf + ISS_CHALLENGE_BASE_BYTES + 2*POLYQISS_VECD_PACKEDBYTES, proof->tB);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_ISS_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_binomial(chal_1[i][j], challenge_seed, DOMAIN_SEPARATOR_CHAL1_ISS, i, j, SEED_BYTES);
    }
  }

  // clean up vectors
  poly_qiss_vec_d_clear(w);
  poly_qiss_vec_d_clear(tmp_vec_d);
}

/*************************************************
* Name:        osig_user_prove_round2 [static]
*
* Description: Compute round 2 of zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
*              - coeff_qiss *chal_2: array of coeff_qiss to host second challenge (gamma_{i,j}) (allocated)
*              - uint8_t *buf: array for XOF input (allocated CHAL2_ISS_INPUT_BYTES bytes)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_vec_k *chal_1: array of polynomial vectors, first challenge R = R0 - R1
*              - const poly_qiss_vec_256_l y3_g: polynomial vector, mask for R.s1 (ARP)
* 
* Returns 1 if round 2 passes, and 0 if it rejects
**************************************************/
static int osig_user_prove_round2(
    osig_proof_t              *proof,
    coeff_qiss                chal_2[PARAM_L_ISS][PARAM_ARP_ISS + 1],
    uint8_t                   buf[CHAL2_ISS_INPUT_BYTES], 
    const poly_qiss_vec_k     s1[PARAM_M1_K_ISS],
    const poly_qiss_vec_k     chal_1[PARAM_ARP_ISS][PARAM_M1_K_ISS],
    const poly_qiss_vec_256_l y3_g) {
  size_t i,j,k,l;
  int64_t tmp;
  uint64_t sq_norm_y3 = 0, sq_norm_z3 = 0;
  coeff_qiss tmp_coeff;
  uint8_t challenge_seed[SEED_BYTES];

  // computing z3 (in Z not in ring) and square norms
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    tmp = 0;
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      for (k = 0; k < PARAM_K_ISS; k++) {
        for (l = 0; l < PARAM_N_ISS; l++) {
          tmp += poly_qiss_get_coeff_centered(chal_1[i][j]->entries[k], l) * poly_qiss_get_coeff_centered(s1[j]->entries[k], l);
        }
      }
    }
    tmp_coeff = poly_qiss_get_coeff_centered(y3_g->entries[i / PARAM_N_ISS], i % PARAM_N_ISS);
    tmp += tmp_coeff; 
    CHK_UI_OVF_ADDITION(sq_norm_y3, tmp_coeff * tmp_coeff);
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp * tmp); 
    proof->z3[i] = tmp;
  }

  // rejection sampling
  // sample u3 uniform in (0,1), goto reject if u3 > exp(pi * (sq_norm_y3 - sq_norm_z3) / PARAM_S3SQ_ISS) / PARAM_REJ3_ISS
  if (_reject_exp(exp(M_PI * (double)(sq_norm_y3 - sq_norm_z3) / (double)PARAM_S3SQ_ISS)/(double)PARAM_REJ3_ISS)) {
    return 0;
  }

  // computing second challenge
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    coeff_qiss_pack(buf + CHAL1_ISS_INPUT_BYTES + i*COEFFQISS_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_ISS_INPUT_BYTES);
  for (i = 0; i < PARAM_L_ISS; i++) {
    vec_qiss_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_ISS, i, SEED_BYTES); // writes PARAM_ARP_ISS + 1 uniform numbers to the first argument
  }
  return 1;
}

/*************************************************
* Name:        osig_user_prove_round3 [static]
*
* Description: Compute round 3 of zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
*              - poly_qiss_vec_l chal_3_l: polynomial vector to host third challenge (mu_i)_{i < l} (initialized)
*              - poly_qiss_vec_k *chal_3_2dk: array of polynomial vectors to host third challenge (mu_i)_{l <= i < 2dk+l} (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL3_ISS_INPUT_BYTES bytes)
*              - poly_qiss_vec_256 *sum_gamma_e_star: array polynomial vectors to host Sum_j gamma_{ij}e_j* (initialized)
*              - poly_qiss_vec_k *sum_gamma_r_star: array polynomial vectors to host Sum_j gamma_{ij}r_j* (initialized)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_vec_k *chal_1: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const coeff_qiss *chal_2: array of coeff_qiss, second challenge (gamma_{i,j})
*              - const poly_qiss_vec_256_l y3_g: polynomial vectors, mask for R.s1 and automorphism eq.
*              - const poly_qiss s1_star_s1_one: polynomial, precomputation of <s1*, s1 - one>
**************************************************/
static void osig_user_prove_round3(
    osig_proof_t              *proof,
    poly_qiss_vec_l           chal_3_l,
    poly_qiss_vec_k           chal_3_2dk[2*PARAM_D],
    uint8_t                   buf[CHAL3_ISS_INPUT_BYTES], 
    poly_qiss_vec_256         sum_gamma_e_star[PARAM_L_ISS],
    poly_qiss_vec_k           sum_gamma_r_star[PARAM_L_ISS][PARAM_M1_K_ISS],
    const poly_qiss_vec_k     s1[PARAM_M1_K_ISS],
    const poly_qiss_vec_k     chal_1[PARAM_ARP_ISS][PARAM_M1_K_ISS],
    const coeff_qiss          chal_2[PARAM_L_ISS][PARAM_ARP_ISS + 1],
    const poly_qiss_vec_256_l y3_g,
    const poly_qiss           s1_star_s1_one) {
  size_t i,j,k;
  coeff_qiss tmp_coeff;
  poly_qiss tmp_poly;
  poly_qiss_vec_k tmp_vec_k;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors and polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_vec_k_init(tmp_vec_k);

  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_set(tmp_poly, y3_g->entries[PARAM_ARP_DIV_N_ISS + i]); // g_i
    tmp_coeff = 0; // g_i constant coefficient is 0

    // sum of -gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_ISS; j++) {
      tmp_coeff -= chal_2[i][j] * proof->z3[j]; // possible overflow for larger parameters
    }
    poly_qiss_set_coeff(tmp_poly, 0, tmp_coeff);

    // sum of gamma_{ij}r_j*
    for (k = 0; k < PARAM_M1_K_ISS; k++) {
      poly_qiss_vec_k_conjugate(sum_gamma_r_star[i][k], chal_1[0][k]);
      poly_qiss_vec_k_mul_scalar(sum_gamma_r_star[i][k], sum_gamma_r_star[i][k], chal_2[i][0]);
      for (j = 1; j < PARAM_ARP_ISS; j++) {
        poly_qiss_vec_k_conjugate(tmp_vec_k, chal_1[j][k]);
        poly_qiss_vec_k_mul_scalar(tmp_vec_k, tmp_vec_k, chal_2[i][j]);
        poly_qiss_vec_k_add(sum_gamma_r_star[i][k], sum_gamma_r_star[i][k], tmp_vec_k);
      }
    }

    // sum of gamma_{ij}e_j* = conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}]))
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      poly_qiss_set_coeff(sum_gamma_e_star[i]->entries[j], 0, chal_2[i][j * PARAM_N_ISS]);
      for (k = 1; k < PARAM_N_ISS; k++) {
        poly_qiss_set_coeff(sum_gamma_e_star[i]->entries[j], k, - chal_2[i][(j + 1) * PARAM_N_ISS - k]); // set conjugate directly
      }
    }

    /**********************************************
    * Computing
    *   h_i = g_i - sum_j gamma_{ij}.z3_j 
    *         + sum_j gamma_{ij}.e_j*.y3 
    *         + sum_j gamma_{ij}.r_j*.s1 
    *         + gamma_{i,256}.<s1*, s1-one>
    **********************************************/
    poly_qiss_set(proof->h->entries[i], tmp_poly); // g_i - sum_j gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) { // <sum_e_gamma_star,y3>
      poly_qiss_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], y3_g->entries[j]);
      poly_qiss_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }
    for (j = 0; j < PARAM_M1_K_ISS; j++) { // <sum_r_gamma_star,s1>
      poly_qiss_vec_k_mul_inner(tmp_poly, sum_gamma_r_star[i][j], s1[j]);
      poly_qiss_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }
    poly_qiss_mul_scalar(tmp_poly, s1_star_s1_one, chal_2[i][PARAM_ARP_ISS]); // gamma_{i,256}.<s1*, s1-one>
    poly_qiss_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
  }

  // computing third challenge
  buf[0] = 3;
  poly_qiss_vec_l_pack(buf + CHAL2_ISS_INPUT_BYTES, proof->h);
  poly_qiss_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_ISS, CHAL3_ISS_INPUT_BYTES);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_ISS_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_uniform(chal_3_2dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, i+PARAM_L_ISS, SEED_BYTES);
  }

  // clean up vectors and polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_vec_k_clear(tmp_vec_k);
}

/*************************************************
* Name:        osig_user_prove_round4 [static]
*
* Description: Compute round 4 of zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
*              - uint8_t *buf: array for XOF input (allocated CHAL4_ISS_INPUT_BYTES bytes)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qiss_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const poly_qiss_vec_m2 b: polynomial vector for b (CRS)
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qiss_vec_256 *sum_gamma_e_star: array polynomial vectors, Sum_j gamma_{ij}e_j*
*              - const poly_qiss_vec_k *sum_gamma_r_star: array polynomial vectors, Sum_j gamma_{ij}r_j*
*              - const poly_qiss_vec_k *y1: array of polynomial vectors, mask for c.s1
*              - const poly_qiss_vec_m2 y2: polynomial vector, mask for c.s2
*              - const coeff_qiss *chal_2: array of coeff_qiss, second challenge (gamma_{i,j})
*              - const poly_qiss_vec_l chal_3_l: polynomial vector, third challenge (mu_i)_{i < l} 
*              - const poly_qiss_vec_k *chal_3_2dk: array of polynomial vectors, third challenge (mu_i)_{l <= i < 2dk+l}
*              - const poly_qiss_vec_k one: polynomial vector with all ones
**************************************************/
static void osig_user_prove_round4(
    osig_proof_t                *proof,
    uint8_t                     buf[CHAL4_ISS_INPUT_BYTES],
    const poly_qiss_vec_k       s1[PARAM_M1_K_ISS],
    const poly_qiss_vec_m2      s2,
    const poly_qiss_mat_256l_m2 Byg, 
    const poly_qiss_vec_m2      b,
    const poly_qiss_mat_k_k     A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k     Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qiss_mat_k_k     D_embed[PARAM_D][PARAM_M], 
    const poly_qiss_vec_256     sum_gamma_e_star[PARAM_L_ISS],
    const poly_qiss_vec_k       sum_gamma_r_star[PARAM_L_ISS][PARAM_M1_K_ISS],
    const poly_qiss_vec_k       y1[PARAM_M1_K_ISS],
    const poly_qiss_vec_m2      y2,
    const coeff_qiss            chal_2[PARAM_L_ISS][PARAM_ARP_ISS + 1],
    const poly_qiss_vec_l       chal_3_l,
    const poly_qiss_vec_k       chal_3_2dk[2*PARAM_D],
    const poly_qiss_vec_k       one) {
  size_t i,j;
  uint32_t kappa_c;
  poly_qiss sum_mu_gamma, tmp_poly, e0, e1, y1_star_one, t0;
  poly_qiss_vec_k tmp_vec_k, tmp_C_y1;
  poly_qiss_vec_256_l tmp_vec_256_l;
  uint8_t challenge_seed[SEED_BYTES];
  
  // init vectors and polynomials
  poly_qiss_vec_k_init(tmp_vec_k);
  poly_qiss_vec_k_init(tmp_C_y1);
  poly_qiss_vec_256_l_init(tmp_vec_256_l);
  poly_qiss_init(sum_mu_gamma);
  poly_qiss_init(tmp_poly);
  poly_qiss_init(e0);
  poly_qiss_init(e1);
  poly_qiss_init(y1_star_one);
  poly_qiss_init(t0);

  // Computing sum_i mu_i gamma_{i,257}
  poly_qiss_mul_scalar(sum_mu_gamma, chal_3_l->entries[0], chal_2[0][PARAM_ARP_ISS]);
  for (i = 1; i < PARAM_L_ISS; i++) {
    poly_qiss_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_ISS]);
    poly_qiss_add(sum_mu_gamma, sum_mu_gamma, tmp_poly);
  }

  /**********************************************
  * Computing garbage terms e0 and e1
  *   e0 = sum_mu_gamma.<y1*, y1> 
  *   
  *   e1 = sum_mu_gamma.(<y1*, s1> + <y1*,s1>* - <y1*, one>)
  *        + sum_{i < l} mu_i.(<sum_gamma_r_star[i], y1> 
  *        - [Byg.y2]_{256/n + i} - <sum_gamma_e_star[i], [Byg.y2]_{:256/n})
  *        + sum_{i < dk} mu_{l+i}.([q1.I | A_embed]y1_{:2dk} + D_embed.y1_{4dk:})_i
  *        + sum_{i < dk} mu_{l+dk+i}.(Ds_embed.y1_{2dk:4dk})_i
  **********************************************/
  poly_qiss_zero(e0);
  poly_qiss_zero(e1);
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_conjugate(tmp_vec_k, y1[i]); // y1*
    poly_qiss_vec_k_mul_inner(tmp_poly, tmp_vec_k, y1[i]);
    poly_qiss_add(e0, e0, tmp_poly); // <y1*, y1>
    poly_qiss_vec_k_mul_inner(tmp_poly, tmp_vec_k, s1[i]);
    poly_qiss_add(e1, e1, tmp_poly); // <y1*, s1>
    poly_qiss_vec_k_mul_inner(tmp_poly, tmp_vec_k, one);
    poly_qiss_add(y1_star_one, y1_star_one, tmp_poly); // <y1*, one>
  }
  poly_qiss_conjugate(tmp_poly, e1); // contains <y1, s1*>
  poly_qiss_add(e1, e1, tmp_poly); // e1 = <y1*, s1> + <y1, s1*>
  poly_qiss_sub(e1, e1, y1_star_one); // e1 = <y1*, s1> + <y1, s1*> - <y1*, one>

  poly_qiss_mul(e0, e0, sum_mu_gamma);
  poly_qiss_mul(e1, e1, sum_mu_gamma);

  poly_qiss_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, y2); // B_{y,g}.y2
  for (i = 0; i < PARAM_L_ISS; i++) {
    // y1_star_one can be used as a temp variable now
    poly_qiss_vec_k_mul_inner(y1_star_one, sum_gamma_r_star[i][0], y1[0]);
    for (j = 1; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_mul_inner(tmp_poly, sum_gamma_r_star[i][j], y1[j]);
      poly_qiss_add(y1_star_one, y1_star_one, tmp_poly);
    }
    poly_qiss_sub(y1_star_one, y1_star_one, tmp_vec_256_l->entries[PARAM_ARP_DIV_N_ISS + i]);
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      poly_qiss_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], tmp_vec_256_l->entries[j]);
      poly_qiss_sub(y1_star_one, y1_star_one, tmp_poly);
    }
    poly_qiss_mul(y1_star_one, y1_star_one, chal_3_l->entries[i]);
    poly_qiss_add(e1, e1, y1_star_one);
  }

  // computing C.y1 and adding to e1
  for (i = 0; i < PARAM_D; i++) {
    // Top part of C.y1
    poly_qiss_vec_k_mul_scalar(tmp_C_y1, y1[i], PARAM_Q1_ISS); // q_1*I_{dk} x y1[:dk]
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, A_embed[i][j], y1[PARAM_D + j]); // part of y1 corresponding to r_{12}
      poly_qiss_vec_k_add(tmp_C_y1, tmp_C_y1, tmp_vec_k);
    }
    for (j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, D_embed[i][j], y1[4*PARAM_D + j]); // part of y1 corresponding to m
      poly_qiss_vec_k_add(tmp_C_y1, tmp_C_y1, tmp_vec_k);
    }

    // adding mu_{l+i} * [[q_1.I | A_embed | 0 | D_embed]y_1]_i to e1
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(tmp_poly, tmp_C_y1->entries[j], chal_3_2dk[i]->entries[j]);
      poly_qiss_add(e1, e1, tmp_poly);
    }

    // Bottom part of C.y1
    poly_qiss_mat_k_k_mul_vec_k(tmp_C_y1, Ds_embed[i][0], y1[2*PARAM_D + 0]);
    for (j = 1; j < 2*PARAM_D; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, Ds_embed[i][j], y1[2*PARAM_D + j]); // part of y1 corresponding to usk
      poly_qiss_vec_k_add(tmp_C_y1, tmp_C_y1, tmp_vec_k);
    }
    // adding mu_{l+dk+j} * [[0 | Ds_embed | 0]y_1]_j to e1
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(tmp_poly, tmp_C_y1->entries[j], chal_3_2dk[PARAM_D + i]->entries[j]);
      poly_qiss_add(e1, e1, tmp_poly);
    }
  }

  // committing to garbage terms
  poly_qiss_vec_m2_mul_inner(t0, b, y2);
  poly_qiss_add(t0, t0, e0);
  poly_qiss_vec_m2_mul_inner(tmp_poly, b, s2);
  poly_qiss_add(proof->t1, tmp_poly, e1);

  // computing fourth challenge
  kappa_c = 0;
  buf[0] = 4;
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES, t0);
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES + POLYQISS_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_ISS_INPUT_BYTES);
  do {
    poly_qiss_sample_challenge(proof->c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_ISS, kappa_c++, SEED_BYTES);
  } while (challenge_size_iss(proof->c) > PARAM_ETA_ISS);
  proof->ctr_c = kappa_c - 1;

  // clean up vectors and polynomials
  poly_qiss_vec_k_clear(tmp_vec_k);
  poly_qiss_vec_k_clear(tmp_C_y1);
  poly_qiss_vec_256_l_clear(tmp_vec_256_l);
  poly_qiss_clear(sum_mu_gamma);
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(e0);
  poly_qiss_clear(e1);
  poly_qiss_clear(y1_star_one);
  poly_qiss_clear(t0);
}

/*************************************************
* Name:        osig_user_prove_round5 [static]
*
* Description: Compute round 5 of zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qiss_vec_k *y1: array of polynomial vectors, mask for c.s1
*              - const poly_qiss_vec_m2 y2: polynomial vector, mask for c.s2
* 
* Returns 1 if round 5 passes, and 0 if it rejects
**************************************************/
static int osig_user_prove_round5(
    osig_proof_t           *proof,
    const poly_qiss_vec_k  s1[PARAM_M1_K_ISS],
    const poly_qiss_vec_m2 s2,
    const poly_qiss_vec_k  y1[PARAM_M1_K_ISS],
    const poly_qiss_vec_m2 y2) {
  size_t i;
  uint64_t sq_norm_y1 = 0, sq_norm_z1 = 0, sq_norm_y2, sq_norm_z2;

  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 and relevant square l2 norms
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_mul_poly_qiss(proof->z1[i], s1[i], proof->c); 
    poly_qiss_vec_k_add(proof->z1[i], proof->z1[i], y1[i]);
    CHK_UI_OVF_ADDITION(sq_norm_y1, poly_qiss_vec_k_norm2(y1[i]));
    CHK_UI_OVF_ADDITION(sq_norm_z1, poly_qiss_vec_k_norm2(proof->z1[i]));
  }
  poly_qiss_vec_m2_mul_poly_qiss(proof->z2, s2, proof->c); 
  poly_qiss_vec_m2_add(proof->z2, proof->z2, y2);
  sq_norm_y2 = poly_qiss_vec_m2_norm2(y2);
  sq_norm_z2 = poly_qiss_vec_m2_norm2(proof->z2);

  // rejection sampling
  // sample u1 uniform in (0,1), goto reject if u1 > exp(pi * (sq_norm_y1 - sq_norm_z1) / PARAM_S1SQ_ISS) / PARAM_REJ1_ISS
  // sample u2 uniform in (0,1), goto reject if u2 > exp(pi * (sq_norm_y2 - sq_norm_z2) / PARAM_S2SQ_ISS) / PARAM_REJ2_ISS
  if (_reject_exp(
    exp(M_PI * (double)(sq_norm_y1 - sq_norm_z1) / (double)PARAM_S1SQ_ISS) / (double)PARAM_REJ1_ISS
  ) || _reject_exp(
    exp(M_PI * (double)(sq_norm_y2 - sq_norm_z2) / (double)PARAM_S2SQ_ISS) / (double)PARAM_REJ2_ISS
  )) {
    return 0;
  }
  return 1;
}

/*************************************************
* Name:        osig_user_prove
*
* Description: Compute zero-knowledge proof of commitment opening
*              and user registration by user (anonymous credentials issuance)
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure (initialized)
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qiss_vec_k *u: array of polynomial vectors, subring embedding of q1.(cmt-upk | upk)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, subring embedding of (r|usk|msg), witness
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void osig_user_prove(
    osig_proof_t            *proof, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qiss_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qiss_vec_k   u[2*PARAM_D], 
    const poly_qiss_vec_k   s1[PARAM_M1_K_ISS], 
    const uint8_t           crs_seed[CRS_SEED_BYTES],
    const uint8_t           seed[SEED_BYTES]) {
  size_t i,j,k;
  uint8_t randomness_seed[SEED_BYTES];
  uint8_t buf[CHAL4_ISS_INPUT_BYTES] = {0};
  uint32_t kappa;
  poly_qiss tmp_poly, s1_star_s1_one;
  poly_qiss_mat_d_k A1[PARAM_M1_K_ISS];
  poly_qiss_mat_d_m2 A2;
  poly_qiss_mat_256l_m2 Byg;
  poly_qiss_vec_k y1[PARAM_M1_K_ISS], one, s1_star_i;
  poly_qiss_vec_k sum_gamma_r_star[PARAM_L_ISS][PARAM_M1_K_ISS];
  poly_qiss_vec_256 sum_gamma_e_star[PARAM_L_ISS];
  poly_qiss_vec_m2 s2, y2, b;
  poly_qiss_vec_256_l y3_g;

  // challenges
  poly_qiss_vec_k chal_1[PARAM_ARP_ISS][PARAM_M1_K_ISS];
  coeff_qiss chal_2[PARAM_L_ISS][PARAM_ARP_ISS + 1];
  poly_qiss_vec_l chal_3_l;
  poly_qiss_vec_k chal_3_2dk[2*PARAM_D];

  // init
  // init polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_init(s1_star_s1_one);
  // init vectors and matrices
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_init(A1[i]);
    poly_qiss_vec_k_init(y1[i]);
  }
  for (i = 0; i < PARAM_L_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_init(sum_gamma_r_star[i][j]);
    }
    poly_qiss_vec_256_init(sum_gamma_e_star[i]);    
  }
  poly_qiss_vec_k_init(one);
  poly_qiss_vec_k_init(s1_star_i);
  poly_qiss_mat_d_m2_init(A2);
  poly_qiss_mat_256l_m2_init(Byg);
  poly_qiss_vec_m2_init(s2);
  poly_qiss_vec_m2_init(y2);
  poly_qiss_vec_m2_init(b);
  poly_qiss_vec_256_l_init(y3_g);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_init(chal_1[i][j]);
    }
  }
  poly_qiss_vec_l_init(chal_3_l);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(chal_3_2dk[i]);
  }

  // generate random secret seed
  randombytes(randomness_seed, SEED_BYTES);

  // expanding CRS
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_uniform(A1[i], crs_seed, DOMAIN_SEPARATOR_A1_ISS, i*PARAM_K_ISS);
  }
  poly_qiss_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_ISS);
  poly_qiss_mat_256l_m2_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_ISS);
  poly_qiss_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_ISS);

  // byte-packing the statement
  // the matrices A_embed, Ds_embed, D_embed are derived from seed
  // the statement is thus (seed || byte-packing(u)).
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_pack(buf + 1 + CRS_SEED_BYTES + SEED_BYTES + i*POLYQISS_VECK_PACKEDBYTES, u[i]);
  }

  // precomputations
  poly_qiss_zero(s1_star_s1_one);
  for (j = 0; j < PARAM_K_ISS; j++) {
    for (k = 0; k < PARAM_N_ISS; k++) {
      poly_qiss_set_coeff(one->entries[j], k, 1); // poly with all one
    }
  }
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_conjugate(s1_star_i, s1[i]); // s1*

    poly_qiss_vec_k_sub(y1[0], s1[i], one); // using y1[0] as temp variable 
    poly_qiss_vec_k_mul_inner(tmp_poly, s1_star_i, y1[0]); // <s1*, s1 - one>
    poly_qiss_add(s1_star_s1_one, s1_star_s1_one, tmp_poly);
  }
  poly_qiss_vec_k_zero(y1[0]); // resetting y1[0]


  kappa = 0;
reject:
  /****** first round ******/
  osig_user_prove_round1(proof, chal_1, buf, s2, y1, y2, y3_g, s1, A1, A2, Byg, randomness_seed, kappa);
  kappa += PARAM_ARP_DIV_N_L_ISS - PARAM_ARP_DIV_N_ISS + 1;

  /****** second round ******/
  if (!osig_user_prove_round2(proof, chal_2, buf, s1, chal_1, y3_g)) {
    goto reject;
  }

  /****** third round ******/
  osig_user_prove_round3(proof, chal_3_l, chal_3_2dk, buf, sum_gamma_e_star, sum_gamma_r_star, s1, chal_1, chal_2, y3_g, s1_star_s1_one);
  
  /****** fourth round ******/
  osig_user_prove_round4(proof, buf, s1, s2, Byg, b, A_embed, Ds_embed, D_embed, sum_gamma_e_star, sum_gamma_r_star, y1, y2, chal_2, chal_3_l, chal_3_2dk, one);

  /****** fifth round ******/
  if (!osig_user_prove_round5(proof, s1, s2, y1, y2)) {
    goto reject;
  }

  // clean up
  // clean up polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(s1_star_s1_one);
  // clean up vectors and matrices
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_clear(A1[i]);
    poly_qiss_vec_k_clear(y1[i]);
  }
  for (i = 0; i < PARAM_L_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_clear(sum_gamma_r_star[i][j]);
    }
    poly_qiss_vec_256_clear(sum_gamma_e_star[i]);    
  }
  poly_qiss_vec_k_clear(one);
  poly_qiss_vec_k_clear(s1_star_i);
  poly_qiss_mat_d_m2_clear(A2);
  poly_qiss_mat_256l_m2_clear(Byg);
  poly_qiss_vec_m2_clear(s2);
  poly_qiss_vec_m2_clear(y2);
  poly_qiss_vec_m2_clear(b);
  poly_qiss_vec_256_l_clear(y3_g);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_clear(chal_1[i][j]);
    }
  }
  poly_qiss_vec_l_clear(chal_3_l);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(chal_3_2dk[i]);
  }
}
