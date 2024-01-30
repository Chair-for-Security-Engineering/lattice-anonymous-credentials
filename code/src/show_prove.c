#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "sep.h"
#include "show.h"
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
  if (udbl > threshold) {
    return 1;
  }
  return 0;
}

/*************************************************
* Name:        show_user_prove_round1 [static]
*
* Description: Compute round 1 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - poly_qshow_vec_m1 *chal_1: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL1_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_m2 s2: polynomial vector to host ABDLOP commitment randomness (initialized)
*              - poly_qshow_vec_m1 y1: polynomial vector to host mask for c.s1 (initialized)
*              - poly_qshow_vec_m2 y2: polynomial vector to host mask for c.s2 (initialized)
*              - poly_qshow_vec_256_l y3_g: polynomial vectors to host mask for R.s1 and automorphism eq. (initialized)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_mat_d_m1 A1: polynomial matrix for A1 (CRS)
*              - const poly_qshow_mat_d_m2 A2: polynomial matrix for A2 (CRS)
*              - const poly_qshow_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const uint8_t *randomness_seed: seed to extract randomness from (allocated SEED_BYTES bytes)
*              - const uint32_t kappa: XOF domain separator due to rejections 
**************************************************/
static void show_user_prove_round1(
    show_proof_t                 *proof,
    poly_qshow_vec_m1            chal_1[PARAM_ARP_SHOW],
    uint8_t                      buf[CHAL1_SHOW_INPUT_BYTES], 
    poly_qshow_vec_m2            s2, 
    poly_qshow_vec_m1            y1, 
    poly_qshow_vec_m2            y2, 
    poly_qshow_vec_256_l         y3_g,
    const poly_qshow_vec_m1      s1,
    const poly_qshow_mat_d_m1    A1,
    const poly_qshow_mat_d_m2    A2,
    const poly_qshow_mat_256l_m2 Byg,
    const uint8_t                randomness_seed[SEED_BYTES],
    const uint32_t               kappa) {
  size_t i, j;
  uint32_t kpp = kappa;
  coeff_qshow tmp_coeff;
  poly_qshow_vec_d w, tmp_vec_d;
  uint8_t binomial_seed[SEED_BYTES];

  // init vectors
  poly_qshow_vec_d_init(w);
  poly_qshow_vec_d_init(tmp_vec_d);

  // sampling ABDLOP commitment randomness
  poly_qshow_vec_m2_binomial(s2, randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_S2_SHOW);

  // sampling Gaussian masks for c.s1, c.s2, R.s1
  poly_qshow_vec_m1_sample_gaussian_s1(y1);
  poly_qshow_vec_m2_sample_gaussian_s2(y2);
  for (i = 0; i < PARAM_ARP_DIV_N_SHOW; i++) {
    for (j = 0; j < PARAM_N_SHOW; j++) {
      tmp_coeff = SampleZ(0, PARAM_S3_SHOW);
      poly_qshow_set_coeff(y3_g->entries[i], j, tmp_coeff);
    }
  }

  // sampling uniform mask for equations with automorphisms
  for (i = PARAM_ARP_DIV_N_SHOW; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
    poly_qshow_uniform_but_zero(y3_g->entries[i], randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_G_SHOW);
  }

  // computing commitments tA = A1.s1 + A2.s2, w = A1.y1 + A2.y2, tB = Byg.s2 + y3_g
  // tA
  poly_qshow_mat_d_m1_mul_vec_m1(proof->tA, A1, s1);
  poly_qshow_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, s2); 
  poly_qshow_vec_d_add(proof->tA, proof->tA, tmp_vec_d);
  // w
  poly_qshow_mat_d_m1_mul_vec_m1(w, A1, y1);
  poly_qshow_mat_d_m2_mul_vec_m2(tmp_vec_d, A2, y2); 
  poly_qshow_vec_d_add(w, w, tmp_vec_d);
  // tB
  poly_qshow_mat_256l_m2_mul_vec_m2(proof->tB, Byg, s2);
  poly_qshow_vec_256_l_add(proof->tB, proof->tB, y3_g);

  // computing first challenge
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof->tA);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof->tB);
  shake256(binomial_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], binomial_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }

  // clean up vectors
  poly_qshow_vec_d_clear(w);
  poly_qshow_vec_d_clear(tmp_vec_d);
}

/*************************************************
* Name:        show_user_prove_round2 [static]
*
* Description: Compute round 2 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - coeff_qshow *chal_2: array of coeff_qshow to host second challenge (gamma_{i,j}) (allocated)
*              - uint8_t *buf: array for XOF input (allocated CHAL2_SHOW_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge R = R0 - R1
*              - const poly_qshow_vec_256_l y3_g: polynomial vector, mask for R.s1 (ARP)
* 
* Returns 1 if round 2 passes, and 0 if it rejects
**************************************************/
static int show_user_prove_round2(
    show_proof_t               *proof,
    coeff_qshow                chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],
    uint8_t                    buf[CHAL2_SHOW_INPUT_BYTES], 
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const poly_qshow_vec_256_l y3_g) {
  size_t i,j,k;
  int64_t tmp;
  uint64_t sq_norm_y3 = 0, sq_norm_z3 = 0;
  uint8_t challenge_seed[SEED_BYTES];
  coeff_qshow tmp_coeff;

  // computing z3 (in Z not in ring) and square norms
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    tmp = 0;
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      for (k = 0; k < PARAM_N_SHOW; k++) {
          tmp += poly_qshow_get_coeff_centered(chal_1[i]->entries[j], k) * poly_qshow_get_coeff_centered(s1->entries[j], k);
      }
    }
    tmp_coeff = poly_qshow_get_coeff_centered(y3_g->entries[i / PARAM_N_SHOW], i % PARAM_N_SHOW);
    tmp += tmp_coeff; 
    CHK_UI_OVF_ADDITION(sq_norm_y3, tmp_coeff * tmp_coeff);
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp * tmp); 
    proof->z3[i] = tmp;
  }

  // rejection sampling
  // sample u3 uniform in (0,1), goto reject if u3 > exp(pi * (sq_norm_y3 - sq_norm_z3) / PARAM_S3SQ_SHOW) / PARAM_REJ3_SHOW
  if (_reject_exp(exp(M_PI * (double)(sq_norm_y3 - sq_norm_z3) / (double)PARAM_S3SQ_SHOW)/(double)PARAM_REJ3_SHOW)) {
    return 0;
  }

  // computing second challenge
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 6 uniform numbers to the first argument
  }
  return 1;
}

/*************************************************
* Name:        show_user_prove_round3 [static]
*
* Description: Compute round 3 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - poly_qshow_vec_l chal_3_l: polynomial vector to host third challenge (mu_i)_{i < l} (initialized)
*              - poly_qshow_vec_k *chal_3_dk: array of polynomial vectors to host third challenge (mu_i)_{l <= i < dk+l} (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL3_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors to host Sum_j gamma_{ij}e_j* (initialized)
*              - poly_qshow_vec_k *sum_gamma_r_star: array polynomial vectors to host Sum_j gamma_{ij}r_j* (initialized)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_256_l y3_g: polynomial vectors, mask for R.s1 and automorphism eq.
*              - const poly_qshow *quadratic_precomp: array of polynomials containing precomputations
**************************************************/
static void show_user_prove_round3(
    show_proof_t               *proof,
    poly_qshow_vec_l           chal_3_l,
    poly_qshow_vec_k           chal_3_dk[PARAM_D],
    uint8_t                    buf[CHAL3_SHOW_INPUT_BYTES], 
    poly_qshow_vec_256         sum_gamma_e_star[PARAM_L_SHOW],
    poly_qshow_vec_m1          sum_gamma_r_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const coeff_qshow          chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],
    const poly_qshow_vec_256_l y3_g,
    const poly_qshow           quadratic_precomp[6]) {
  size_t i,j,k;
  poly_qshow tmp_poly;
  poly_qshow_vec_m1 tmp_vec_m1, rj_star;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_vec_m1_init(tmp_vec_m1);
  poly_qshow_vec_m1_init(rj_star);

  // sum of gamma_{ij}r_j*
  for (j = 0; j < PARAM_ARP_SHOW; j++) {
    poly_qshow_vec_m1_conjugate(rj_star, chal_1[j]);
    for (i = 0; i < PARAM_L_SHOW; i++) {
      if (j == 0) {
        poly_qshow_vec_m1_mul_scalar(sum_gamma_r_star[i], rj_star, chal_2[i][j]);
      } else {
        poly_qshow_vec_m1_mul_scalar(tmp_vec_m1, rj_star, chal_2[i][j]);
        poly_qshow_vec_m1_add(sum_gamma_r_star[i], sum_gamma_r_star[i], tmp_vec_m1);
      }
    }    
  }

  // sum of gamma_{ij}e_j* = conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}]))
  for (i = 0; i < PARAM_L_SHOW; i++) {
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_set_coeff(sum_gamma_e_star[i]->entries[j], 0, chal_2[i][j * PARAM_N_SHOW]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_e_star[i]->entries[j], k, - chal_2[i][(j + 1) * PARAM_N_SHOW - k]); // set conjugate directly
      }
    }
  }

  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_set(tmp_poly, y3_g->entries[PARAM_ARP_DIV_N_SHOW + i]); // g_i
    // sum of -gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(tmp_poly, chal_2[i][j], -proof->z3[j]); // adds -(chal_2[i][j] * proof->z3[j]) to the constant coefficient without overflow
    }
    
    /**********************************************
    * Computing
    *   h_i = g_i - sum_j gamma_{ij}.z3_j 
    *         + sum_j gamma_{ij}.e_j*.y3 
    *         + sum_j gamma_{ij}.r_j*.s1 
    *         + gamma_{i,256}.(<v1*,v1> - B1^2)
    *         + gamma_{i,257}.(<v2*,v2> - B2^2)
    *         + gamma_{i,258}.(<v3*,v3> - B3^2)
    *         + gamma_{i,259}.(<t*,t> - w)
    *         + gamma_{i,260}.<t*,t - one>
    *         + gamma_{i,261}.<(usk|msg)*,(usk|msg) - one>
    **********************************************/
    poly_qshow_set(proof->h->entries[i], tmp_poly); // g_i - sum_j gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) { // <sum_e_gamma_star,y3>
      poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], y3_g->entries[j]);
      poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }
    poly_qshow_vec_m1_mul_inner(tmp_poly, sum_gamma_r_star[i], s1); // <sum_r_gamma_star,s1>
    poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    for (j = 0; j < 6; j++) { // precomputed quadratic terms
      poly_qshow_mul_scalar(tmp_poly, quadratic_precomp[j], chal_2[i][PARAM_ARP_SHOW + j]);
      poly_qshow_add(proof->h->entries[i], proof->h->entries[i], tmp_poly);
    }
  }

  // computing third challenge
  buf[0] = 3;
  poly_qshow_vec_l_pack(buf + CHAL2_SHOW_INPUT_BYTES, proof->h);
  poly_qshow_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_SHOW, CHAL3_SHOW_INPUT_BYTES);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_SHOW, i+PARAM_L_SHOW, SEED_BYTES);
  }

  // clean up vectors and polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_vec_m1_clear(tmp_vec_m1);
  poly_qshow_vec_m1_clear(rj_star);
}

/*************************************************
* Name:        show_user_prove_round4 [static]
*
* Description: Compute round 4 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
*              - uint8_t *buf: array for XOF input (allocated CHAL4_ISS_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const poly_qshow_vec_m2 b: polynomial vector for b (CRS)
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors, Sum_j gamma_{ij}e_j*
*              - const poly_qshow_vec_k *sum_gamma_r_star: array polynomial vectors, Sum_j gamma_{ij}r_j*
*              - const poly_qshow_vec_m1 y1: polynomial vector, mask for c.s1
*              - const poly_qshow_vec_m2 y2: polynomial vector, mask for c.s2
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_l chal_3_l: polynomial vector, third challenge (mu_i)_{i < l} 
*              - const poly_qshow_vec_k *chal_3_dk: array of polynomial vectors, third challenge (mu_i)_{l <= i < dk+l}
*              - const poly_qshow one: polynomial with all ones
**************************************************/
static void show_user_prove_round4(
    show_proof_t                *proof,
    uint8_t                      buf[CHAL4_SHOW_INPUT_BYTES],
    const poly_qshow_vec_m1      s1,
    const poly_qshow_vec_m2      s2, 
    const poly_qshow_mat_256l_m2 Byg,
    const poly_qshow_vec_m2      b,
    const poly_qshow_mat_k_k     A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k     B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k     A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k     Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k     D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_vec_256     sum_gamma_e_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1      sum_gamma_r_star[PARAM_L_SHOW],
    const poly_qshow_vec_m1      y1,
    const poly_qshow_vec_m2      y2,
    const coeff_qshow            chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6],
    const poly_qshow_vec_l       chal_3_l,
    const poly_qshow_vec_k       chal_3_dk[PARAM_D],
    const poly_qshow             one) {
  size_t i,j,k,i_k_quot,i_k_rem;
  uint32_t kappa_c;
  int64_t bexpi;
  uint8_t challenge_seed[SEED_BYTES];
  poly_qshow tmp_poly, e0, e1, y1i_star, y1s_y1, y1s_s1, y1s_one, t0, sum_mu_gamma[6];
  poly_qshow_vec_256_l tmp_vec_256_l;
  poly_qshow_vec_k tmp_vec_k, Gy1_v2[PARAM_D], Gs1_v2[PARAM_D];
  poly_qshow_mat_k_k chal_3_quad_matrix[PARAM_D];
  
  // init matrices, vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(e0);
  poly_qshow_init(e1);
  poly_qshow_init(y1i_star);
  poly_qshow_init(y1s_y1);
  poly_qshow_init(y1s_s1);
  poly_qshow_init(y1s_one);
  poly_qshow_init(t0);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_k_init(tmp_vec_k);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_init(chal_3_quad_matrix[i]);
    poly_qshow_vec_k_init(Gy1_v2[i]);
    poly_qshow_vec_k_zero(Gy1_v2[i]);
    poly_qshow_vec_k_init(Gs1_v2[i]);
    poly_qshow_vec_k_zero(Gs1_v2[i]);
  }
  for (i = 0; i < 6; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }

  // compute gadget quadratic matrix necessary to compute sum_i mu_{l + i} G_i"
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_chal_3_embed(chal_3_quad_matrix[i], chal_3_dk[i]);
  }

  // compute sum_i mu_i gamma_{i,256 + j} for j in [6]
  for (i = 0; i < PARAM_L_SHOW; i++) {
    for (j = 0; j < 6; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_SHOW + j]);
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
  }
  poly_qshow_add(sum_mu_gamma[3], sum_mu_gamma[3], sum_mu_gamma[4]); // sum_mu_gamma[4] will be used later

  // quadratic terms
  // <y1_v1"*,y1_v1">, <y1_v1"*,s1_v1">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(e0, sum_mu_gamma[0], y1s_y1); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y1*,s1> + <y1*,s1>*
  poly_qshow_mul(e1, sum_mu_gamma[0], y1s_s1);


  // <y1_v2"*,y1_v2">, <y1_v2"*,s1_v2">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);

  // <y1_v3"*,y1_v3">, <y1_v3"*,s1_v3">
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);

  // <y1_t'*,y1_t'>, <y1_t'*,s1_t'>, <y1_t'*,one>
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  poly_qshow_zero(y1s_one);
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, one);
    poly_qshow_add(y1s_one, y1s_one, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);
  poly_qshow_mul(tmp_poly, sum_mu_gamma[4], y1s_one); 
  poly_qshow_sub(e1, e1, tmp_poly); // substract linear part containing sum_mu_gamma[4].<y1*, one>

  // <y1_m'*,y1_m'>, <y1_m'*,s1_m'>, <y1_m'*,one> 
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  poly_qshow_zero(y1s_one);
  for (i = IDX_USK_SHOW; i < PARAM_M1_SHOW; i++) {// same for usk and ALL message attributes --> no selective disclosure (!)
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, one);
    poly_qshow_add(y1s_one, y1s_one, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[5], y1s_y1); 
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(tmp_poly, tmp_poly, y1s_s1);
  poly_qshow_sub(tmp_poly, tmp_poly, y1s_one); // contains <y1_m*, s1_m> + <y1_m*, s1_m>* - <y1_m*, one>
  poly_qshow_mul(tmp_poly, sum_mu_gamma[5], tmp_poly); 
  poly_qshow_add(e1, e1, tmp_poly);
  
  // Linear part depending on sum_gamma_r_star and sum_gamma_e_star
  poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, y2);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    // y1s_y1 can be used as a temp variable
    poly_qshow_vec_m1_mul_inner(y1s_y1, sum_gamma_r_star[i], y1);
    poly_qshow_sub(y1s_y1, y1s_y1, tmp_vec_256_l->entries[PARAM_ARP_DIV_N_SHOW + i]);
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], tmp_vec_256_l->entries[j]);
      poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly);
    }
    poly_qshow_mul(y1s_y1, y1s_y1, chal_3_l->entries[i]);
    poly_qshow_add(e1, e1, y1s_y1);
  }
  
  // quadratic part depending on chal_3_quad_matrix
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gy1_v2[i]->entries[j]);
      poly_qshow_zero(Gs1_v2[i]->entries[j]);
      bexpi = PARAM_Q1_SHOW;
      for (k = 0; k < PARAM_K; k++) {
        poly_qshow_mul_scalar(tmp_poly, y1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gy1_v2[i]->entries[j], Gy1_v2[i]->entries[j], tmp_poly);
        poly_qshow_mul_scalar(tmp_poly, s1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gs1_v2[i]->entries[j], Gs1_v2[i]->entries[j], tmp_poly);
        bexpi *= PARAM_B;
      }
    }
  }
  // computing (sum_i mu_{l+i}G_i").y1_{v2} within the first vec_k of Gy1_v2, and (sum_i mu_{l+i}G_i").s1_{v2} within the first vec_k of Gs1_v2
  poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[0], Gy1_v2[0]);
  poly_qshow_vec_k_set(Gy1_v2[0], tmp_vec_k);
  poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[0], Gs1_v2[0]);
  poly_qshow_vec_k_set(Gs1_v2[0], tmp_vec_k);
  for (i = 1; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[i], Gy1_v2[i]);
    poly_qshow_vec_k_add(Gy1_v2[0], Gy1_v2[0], tmp_vec_k);
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[i], Gs1_v2[i]);
    poly_qshow_vec_k_add(Gs1_v2[0], Gs1_v2[0], tmp_vec_k);
  }
  // <y1_t, (sum_i mu_{l+i}G_i").y1_{v2}>, <s1_t, (sum_i mu_{l+i}G_i").y1_{v2}>, and <y1_t, (sum_i mu_{l+i}G_i").s1_{v2}> and adding to e0/e1
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul(tmp_poly, y1->entries[IDX_TAG_SHOW + i], Gy1_v2[0]->entries[i]);
    poly_qshow_add(e0, e0, tmp_poly);
    poly_qshow_mul(tmp_poly, y1->entries[IDX_TAG_SHOW + i], Gs1_v2[0]->entries[i]);
    poly_qshow_add(e1, e1, tmp_poly);
    poly_qshow_mul(tmp_poly, s1->entries[IDX_TAG_SHOW + i], Gy1_v2[0]->entries[i]);
    poly_qshow_add(e1, e1, tmp_poly);
  }

  // linear part depending on embedded matrices A', B, A3, Ds, D
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    i_k_quot = i / PARAM_K_SHOW;
    i_k_rem = i % PARAM_K_SHOW;
    // y1s_y1 used as temp variable to host ([A|-B|A3|-Ds|-D].[y1_v1|y1_v2|y1_v3|y1_s|y1_m])_i
    poly_qshow_mul_scalar(y1s_y1, y1->entries[IDX_V1_SHOW + i], PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_V12_SHOW + j]);
      poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    }
    for (j = 0; j < PARAM_D*PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, B_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_V2_SHOW + j]);
      poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A3_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_V3_SHOW + j]);
      poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    }
    for (j = 0; j < 2*PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, Ds_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_USK_SHOW + j]);
      poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_M*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, D_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], y1->entries[IDX_M_SHOW + j]);
      poly_qshow_sub(y1s_y1, y1s_y1, tmp_poly); // substraction
    }
    poly_qshow_mul(y1s_y1, y1s_y1, chal_3_dk[i_k_quot]->entries[i_k_rem]);
    poly_qshow_add(e1, e1, y1s_y1);
  }

  // committing to garbage terms
  poly_qshow_vec_m2_mul_inner(t0, b, y2);
  poly_qshow_add(t0, t0, e0);
  poly_qshow_vec_m2_mul_inner(tmp_poly, b, s2);
  poly_qshow_add(proof->t1, tmp_poly, e1);

  // computing fourth challenge
  kappa_c = 0;
  buf[0] = 4;
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  do {
    poly_qshow_sample_challenge(proof->c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, kappa_c++, SEED_BYTES);
  } while (challenge_size_show(proof->c) > PARAM_ETA_SHOW);
  proof->ctr_c = kappa_c - 1;

  // clean up matrices, vectors and polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(e0);
  poly_qshow_clear(e1);
  poly_qshow_clear(y1i_star);
  poly_qshow_clear(y1s_y1);
  poly_qshow_clear(y1s_s1);
  poly_qshow_clear(y1s_one);
  poly_qshow_clear(t0);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_vec_k_clear(tmp_vec_k);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_clear(chal_3_quad_matrix[i]);
    poly_qshow_vec_k_clear(Gy1_v2[i]);
    poly_qshow_vec_k_clear(Gs1_v2[i]);
  }
  for (i = 0; i < 6; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
  }
}

/*************************************************
* Name:        show_user_prove_round5 [static]
*
* Description: Compute round 5 of zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof_t *proof: pointer to showuance proof structure
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m2 s2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_vec_m1 y1: polynomial vector, mask for c.s1
*              - const poly_qshow_vec_m2 y2: polynomial vector, mask for c.s2
* 
* Returns 1 if round 5 passes, and 0 if it rejects
**************************************************/
static int show_user_prove_round5(
    show_proof_t            *proof,
    const poly_qshow_vec_m1 s1,
    const poly_qshow_vec_m2 s2,
    const poly_qshow_vec_m1 y1,
    const poly_qshow_vec_m2 y2) {
  uint64_t sq_norm_y1, sq_norm_z1, sq_norm_y2, sq_norm_z2;

  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 (and square norms)
  poly_qshow_vec_m1_mul_poly_qshow(proof->z1, s1, proof->c); 
  poly_qshow_vec_m1_add(proof->z1, proof->z1, y1);
  sq_norm_y1 = poly_qshow_vec_m1_norm2(y1);
  sq_norm_z1 = poly_qshow_vec_m1_norm2(proof->z1);

  poly_qshow_vec_m2_mul_poly_qshow(proof->z2, s2, proof->c); 
  poly_qshow_vec_m2_add(proof->z2, proof->z2, y2);
  sq_norm_y2 = poly_qshow_vec_m2_norm2(y2);
  sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);

  // rejection sampling
  // sample u1 uniform in (0,1), goto reject if u1 > exp(pi * (sq_norm_y1 - sq_norm_z1) / PARAM_S1SQ_SHOW) / PARAM_REJ1_SHOW
  // sample u2 uniform in (0,1), goto reject if u2 > exp(pi * (sq_norm_y2 - sq_norm_z2) / PARAM_S2SQ_SHOW) / PARAM_REJ2_SHOW
  if (_reject_exp(
    exp(M_PI * (double)(sq_norm_y1 - sq_norm_z1) / (double)PARAM_S1SQ_SHOW) / (double)PARAM_REJ1_SHOW
  ) || _reject_exp(
    exp(M_PI * (double)(sq_norm_y2 - sq_norm_z2) / (double)PARAM_S2SQ_SHOW) / (double)PARAM_REJ2_SHOW
  )) {
    return 0;
  }
  return 1;
}

/*************************************************
* Name:        show_user_prove
*
* Description: Compute zero-knowledge proof of valid credential
*              (anonymous credentials show)
*
* Arguments:   - show_proof *proof: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*                   s1 = [theta(v1)|a1 | theta(v2)|a2 | theta(v3)|a3 | theta(tag) | theta(usk) | theta(msg)]
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void show_user_prove(
    show_proof_t             *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i;
  uint8_t randomness_seed[SEED_BYTES];
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0};
  uint32_t kappa;
  poly_qshow s1i_star, tmp_poly, one, quadratic_precomp[6];
  poly_qshow_mat_d_m1 A1;
  poly_qshow_mat_d_m2 A2;
  poly_qshow_mat_256l_m2 Byg;
  poly_qshow_vec_m1 y1;
  poly_qshow_vec_m1 sum_gamma_r_star[PARAM_L_SHOW];
  poly_qshow_vec_256 sum_gamma_e_star[PARAM_L_SHOW];
  poly_qshow_vec_m2 s2, y2, b;
  poly_qshow_vec_256_l y3_g;

  // challenges
  poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
  coeff_qshow chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6];
  poly_qshow_vec_l chal_3_l;
  poly_qshow_vec_k chal_3_dk[PARAM_D];

  // init
  // init polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(s1i_star);
  poly_qshow_init(one);
  for (i = 0; i < 6; i++) {
    poly_qshow_init(quadratic_precomp[i]);
    poly_qshow_zero(quadratic_precomp[i]);
  }
  // init vectors and matrices
  poly_qshow_mat_d_m1_init(A1);
  poly_qshow_vec_m1_init(y1);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_init(sum_gamma_r_star[i]);
    poly_qshow_vec_256_init(sum_gamma_e_star[i]);    
  }
  poly_qshow_mat_d_m2_init(A2);
  poly_qshow_mat_256l_m2_init(Byg);
  poly_qshow_vec_m2_init(s2);
  poly_qshow_vec_m2_init(y2);
  poly_qshow_vec_m2_init(b);
  poly_qshow_vec_256_l_init(y3_g);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_init(chal_1[i]);
  }
  poly_qshow_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(chal_3_dk[i]);
  }

  // generate random secret seed
  randombytes(randomness_seed, SEED_BYTES);

  // expanding CRS
  poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_SHOW);
  poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_SHOW);
  poly_qshow_mat_256l_m2_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_SHOW);
  poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_SHOW);

  // byte-packing the statement
  // public parameters are derived from seed
  // the signer public key does not have to be hash (single-issuer setting)
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }

  // precomputations for quadratic equations
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(one, i, 1);
  }
  // <v1"*,v1"> - B_1'^2
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[0], quadratic_precomp[0], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[0], -PARAM_B1SQ, 1); 

  // <v2"*,v2"> - B_2^2
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[1], quadratic_precomp[1], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[1], -PARAM_B2SQ, 1); 

  // <v3"*,v3"> - B_3^2
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[2], quadratic_precomp[2], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[2], -PARAM_B3SQ, 1); 

  // <tag'*,tag'> - w and <tag'*,tag' - one>
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_mul(tmp_poly, s1i_star, s1->entries[i]);
    poly_qshow_add(quadratic_precomp[3], quadratic_precomp[3], tmp_poly); // <tag'*,tag'>
    poly_qshow_sub(tmp_poly, s1->entries[i], one);
    poly_qshow_mul(tmp_poly, s1i_star, tmp_poly);
    poly_qshow_add(quadratic_precomp[4], quadratic_precomp[4], tmp_poly); // <tag'*,tag' - one>
  }
  poly_qshow_muladd_constant(quadratic_precomp[3], -PARAM_W, 1); 

  // <m'*,m' - one> (includes usk and message)
  for (i = IDX_USK_SHOW; i < PARAM_M1_SHOW; i++) {// same for usk and ALL message attributes --> no selective disclosure (!)
    poly_qshow_conjugate(s1i_star, s1->entries[i]);
    poly_qshow_sub(tmp_poly, s1->entries[i], one);
    poly_qshow_mul(tmp_poly, s1i_star, tmp_poly);
    poly_qshow_add(quadratic_precomp[5], quadratic_precomp[5], tmp_poly);
  }

  kappa = 0;
reject:
  /****** first round ******/
  show_user_prove_round1(proof, chal_1, buf, s2, y1, y2, y3_g, s1, A1, A2, Byg, randomness_seed, kappa);
  kappa += PARAM_ARP_DIV_N_L_SHOW - PARAM_ARP_DIV_N_SHOW + 1;

  /****** second round ******/
  if (!show_user_prove_round2(proof, chal_2, buf, s1, chal_1, y3_g)) {
    goto reject;
  }

  /****** third round ******/
  show_user_prove_round3(proof, chal_3_l, chal_3_dk, buf, sum_gamma_e_star, sum_gamma_r_star, s1, chal_1, chal_2, y3_g, quadratic_precomp);
  
  /****** fourth round ******/
  show_user_prove_round4(proof, buf, s1, s2, Byg, b, A_embed, B_embed, A3_embed, Ds_embed, D_embed, sum_gamma_e_star, sum_gamma_r_star, y1, y2, chal_2, chal_3_l, chal_3_dk, one);

  /****** fifth round ******/
  if (!show_user_prove_round5(proof, s1, s2, y1, y2)) {
    goto reject;
  }
  
  // clean up
  // clean up polynomials
  poly_qshow_clear(s1i_star);
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(one);
  for (i = 0; i < 6; i++) {
    poly_qshow_clear(quadratic_precomp[i]);
  }
  // clean up vectors and matrices
  poly_qshow_mat_d_m1_clear(A1);
  poly_qshow_vec_m1_clear(y1);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_clear(sum_gamma_r_star[i]);
    poly_qshow_vec_256_clear(sum_gamma_e_star[i]);    
  }
  poly_qshow_mat_d_m2_clear(A2);
  poly_qshow_mat_256l_m2_clear(Byg);
  poly_qshow_vec_m2_clear(s2);
  poly_qshow_vec_m2_clear(y2);
  poly_qshow_vec_m2_clear(b);
  poly_qshow_vec_256_l_clear(y3_g);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_clear(chal_1[i]);
  }
  poly_qshow_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(chal_3_dk[i]);
  }
}
