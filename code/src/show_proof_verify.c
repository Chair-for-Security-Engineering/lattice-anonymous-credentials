#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "sep.h"
#include "show.h"
#include "macros.h"
#include "fips202.h"

/*************************************************
* Name:        show_verify
*
* Description: Verification of the zero-knowledge proof of valid 
*              credentials (anonymous credentials show)
*
* Arguments:   - const osig_proof *proof: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices, subring embedding of q1.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices, subring embedding of q1.A3
*              - const poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qshow_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qshow_vec_k *u_embed: array of polynomial vectors, subring embedding of q1.u
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
* 
* Returns 1 if proof could be verified correctly and 0 otherwise
**************************************************/
int show_verify(
    const show_proof_t       *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i,j,k,i_k_quot,i_k_rem;
  int is_valid = 1;
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0}, challenge_seed[SEED_BYTES];
  uint64_t sq_norm_z3, sq_norm_z2;
  uint128 sq_norm_z1;
  int64_t bexpi;
  coeff_qshow chal_2[PARAM_L_SHOW][PARAM_ARP_SHOW + 6];
  coeff_qshow tmp_coeff;
  poly_qshow tmp_poly, z1s_z1, z1s_c_one, sum_mu_gamma[7], t0, c, c_one;
  poly_qshow_mat_d_m1 A1;
  poly_qshow_mat_d_m2 A2;
  poly_qshow_mat_256l_m2 Byg;
  poly_qshow_vec_d tmp_vec_d, w;
  poly_qshow_vec_k tmp_vec_k; 
  poly_qshow_vec_m1 sum_gamma_r_star_i, tmp_vec_m1;
  poly_qshow_vec_k Gz1_v2[PARAM_D], chal_3_dk[PARAM_D];
  poly_qshow sum_gamma_e_star_ij;
  poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
  poly_qshow_vec_m2 b;
  poly_qshow_vec_256_l tmp_vec_256_l, Z;
  poly_qshow_vec_l chal_3_l;
  poly_qshow_mat_k_k chal_3_quad_matrix[PARAM_D];

  // init
  // init polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(z1s_z1);
  poly_qshow_init(z1s_c_one);
  poly_qshow_init(c);
  poly_qshow_init(c_one);
  poly_qshow_init(t0);
  // init vectors and matrices
  poly_qshow_mat_d_m1_init(A1);
  poly_qshow_mat_d_m2_init(A2);
  poly_qshow_mat_256l_m2_init(Byg);
  poly_qshow_vec_k_init(tmp_vec_k);
  poly_qshow_vec_m1_init(tmp_vec_m1);
  poly_qshow_vec_m1_init(sum_gamma_r_star_i);
  poly_qshow_vec_d_init(tmp_vec_d);
  poly_qshow_vec_d_init(w);
  poly_qshow_init(sum_gamma_e_star_ij);    
  poly_qshow_vec_m2_init(b);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_256_l_init(Z);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_init(chal_1[i]);
  }
  poly_qshow_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_init(chal_3_quad_matrix[i]);
    poly_qshow_vec_k_init(chal_3_dk[i]);
    poly_qshow_vec_k_init(Gz1_v2[i]);
  }
  for (i = 0; i < 7; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }
  ////////////////////////////////////// NO INITIALIZATIONS BELOW THIS LINE

  // expanding CRS
  poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_SHOW);
  poly_qshow_mat_d_m2_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_SHOW);
  poly_qshow_mat_256l_m2_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_SHOW);
  poly_qshow_vec_m2_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_SHOW);

  // byte-packing the statement (only seeds)
  // the signer public key does not have to be hash (single-issuer setting)
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }

  // precomputation
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(c_one, i, 1); // poly with all one
  }
  poly_qshow_mul(c_one, c_one, proof->c); // holds c.one

  /****** reconstructing first message ******/
  // recomputing w = A1.z1 + A2.z2 - c.tA
  poly_qshow_vec_d_mul_poly_qshow(tmp_vec_d, proof->tA, proof->c);
  poly_qshow_mat_d_m2_mul_vec_m2(w, A2, proof->z2);
  poly_qshow_vec_d_sub(w, w, tmp_vec_d);
  poly_qshow_mat_d_m1_mul_vec_m1(tmp_vec_d, A1, proof->z1);
  poly_qshow_vec_d_add(w, w, tmp_vec_d);
  // recomputing first challenge
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof->tA);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof->tB);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], challenge_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }

  /****** reconstructing second message ******/
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 6 uniform numbers to the first argument
  }

  /****** reconstructing third message ******/
  buf[0] = 3;
  poly_qshow_vec_l_pack(buf + CHAL2_SHOW_INPUT_BYTES, proof->h);
  poly_qshow_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_SHOW, CHAL3_SHOW_INPUT_BYTES);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_SHOW, i+PARAM_L_SHOW, SEED_BYTES);
  }

  /****** reconstructing fourth message ******/
  /**********************************************
  * Recomputing commitment t0
  *   Z = c.tB - Byg.z2
  *   
  *   t0 = zFz + c.<f,z> + c^2.f + <b,z2> - c.t1
  *      = <b, z2> - c.t1 + sum_mu_gamma[0].(<z1_v1*, z1_v1> - c^2.B1^2)
  *         + sum_mu_gamma[1].(<z1_v2*, z1_v2> - c^2.B2^2)
  *         + sum_mu_gamma[2].(<z1_v3*, z1_v3> - c^2.B3^2) 
  *         + sum_mu_gamma[6].<z1_t*, z1_t>
  *         - sum_mu_gamma[3].c^2.w - sum_mu_gamma[4].<z1_t*, c.one>
  *         + sum_mu_gamma[5].<z1_m*, z1_m - c.one>
  *         + c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n + i}
  *         + <sum_gamma_e_star[i], Z_{:256/n}> - c.(sum_gamma_z3 + hi))
  *         + <z1_t, sum_{i < dk} mu_{l+i} G_i".z1_v2> 
  *         + c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed].z1_v1 - B_embed.z1_v2 + A3_embed.z1_v3 - Ds_embed.z1_usk - D_embed.z1_m - c.u_embed)_i
  **********************************************/
  // computing <b, z2> - c.t1
  poly_qshow_vec_m2_mul_inner(t0, b, proof->z2);
  poly_qshow_mul(tmp_poly, proof->c, proof->t1);
  poly_qshow_sub(t0, t0, tmp_poly);

  // c holds proof->c^2 for now
  poly_qshow_mul(c, proof->c, proof->c);

  // compute sum_i mu_i gamma_{i,256 + j} for j in [6]
  for (i = 0; i < PARAM_L_SHOW; i++) {
    for (j = 0; j < 6; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_SHOW + j]);
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
  }
  poly_qshow_add(sum_mu_gamma[6], sum_mu_gamma[3], sum_mu_gamma[4]); // both sum_mu_gamma[3] and sum_mu_gamma[4] will be used later

  // computing Z = c.tB - Byg.z2
  poly_qshow_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, proof->z2);
  poly_qshow_vec_256_l_mul_poly_qshow(Z, proof->tB, proof->c);
  poly_qshow_vec_256_l_sub(Z, Z, tmp_vec_256_l);

  // quadratic terms
  // <z1_v1"*,z1_v1"> - c^2 B1'^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V1_SHOW; i < IDX_V2_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B1SQ); // c^2 B1'^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[0]); 
  poly_qshow_add(t0, t0, z1s_z1); 

  // <z1_v2"*,z1_v2"> - c^2 B2^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V2_SHOW; i < IDX_V3_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B2SQ); // c^2 B2^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[1]); 
  poly_qshow_add(t0, t0, z1s_z1); 

  // <z1_v3"*,z1_v3"> - c^2 B3^2
  poly_qshow_zero(z1s_z1);
  for (i = IDX_V3_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(tmp_poly, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_B3SQ); // c^2 B3^2
  poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[2]); 
  poly_qshow_add(t0, t0, z1s_z1); 

  // -c^2.w
  poly_qshow_mul_scalar(tmp_poly, c, PARAM_W); // c^2 w
  poly_qshow_mul(tmp_poly, tmp_poly, sum_mu_gamma[3]);
  poly_qshow_sub(t0, t0, tmp_poly); // substract

  // no more need for c^2, c can be used as temp variable
  // <z1_t'*,z1_t'> and <z1_t'*,c.one>
  poly_qshow_zero(z1s_z1);
  poly_qshow_zero(z1s_c_one);
  for (i = IDX_TAG_SHOW; i < IDX_USK_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_mul(c, tmp_poly, proof->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, c);
    poly_qshow_mul(c, tmp_poly, c_one);
    poly_qshow_add(z1s_c_one, z1s_c_one, c);
  }
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[6]);
  poly_qshow_add(t0, t0, z1s_z1);
  poly_qshow_mul(z1s_c_one, z1s_c_one, sum_mu_gamma[4]);
  poly_qshow_sub(t0, t0, z1s_c_one); // substract

  // z1s_c_one can be used as temp variable
  // <z1_m'*,z1_m' - c.one>
  poly_qshow_zero(z1s_z1);
  for (i = IDX_USK_SHOW; i < PARAM_M1_SHOW; i++) { // same for usk and ALL message attributes --> no selective disclosure (!)
    poly_qshow_conjugate(tmp_poly, proof->z1->entries[i]);
    poly_qshow_sub(z1s_c_one, proof->z1->entries[i], c_one);
    poly_qshow_mul(tmp_poly, tmp_poly, z1s_c_one);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul(z1s_z1, z1s_z1, sum_mu_gamma[5]);
  poly_qshow_add(t0, t0, z1s_z1);

  // quadratic part depending on chal_3_quad_matrix
  // compute gadget quadratic matrix necessary to compute sum_i mu_{l + i} G_i"
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_chal_3_embed(chal_3_quad_matrix[i], chal_3_dk[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gz1_v2[i]->entries[j]);
      bexpi = PARAM_Q1_SHOW;
      for (k = 0; k < PARAM_K; k++) {
        poly_qshow_mul_scalar(tmp_poly, proof->z1->entries[IDX_V2_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gz1_v2[i]->entries[j], Gz1_v2[i]->entries[j], tmp_poly);
        bexpi *= PARAM_B;
      }
    }
  }
  // computing (sum_i mu_{l+i}G_i").z1_{v2} within the first vec_k of Gz1_v2
  poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[0], Gz1_v2[0]);
  poly_qshow_vec_k_set(Gz1_v2[0], tmp_vec_k);
  for (i = 1; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix[i], Gz1_v2[i]);
    poly_qshow_vec_k_add(Gz1_v2[0], Gz1_v2[0], tmp_vec_k);
  }
  // <z1_t, (sum_i mu_{l+i}G_i").z1_{v2}>
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul(tmp_poly, proof->z1->entries[IDX_TAG_SHOW + i], Gz1_v2[0]->entries[i]);
    poly_qshow_add(t0, t0, tmp_poly);
  }

  // Linear part depending on sum_gamma_r_star and sum_gamma_e_star
  // computing Y = c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n + i} + <sum_gamma_e_star[i], Z_{:256/n}> - c.(sum_gamma_z3 + hi))
  // conjugating all chal_1
  for (j = 0; j < PARAM_ARP_SHOW; j++) {
    poly_qshow_vec_m1_conjugate(tmp_vec_m1, chal_1[j]);
    poly_qshow_vec_m1_set(chal_1[j], tmp_vec_m1);
  }
  
  // z1s_z1, z1s_c_one can be used as temp variables
  poly_qshow_zero(z1s_c_one); // will hold Y/c
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_set(tmp_poly, proof->h->entries[i]); // h_i
    is_valid = is_valid && (poly_qshow_get_coeff(tmp_poly, 0) == 0);
    if (!is_valid) {
      goto show_verify_cleanup; // h_i does not have a zero coefficient, is_valid = 0
    }
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(tmp_poly, chal_2[i][j], proof->z3[j]); // adds (chal_2[i][j] * proof->z3[j]) to the constant coefficient without overflow
    }
    poly_qshow_mul(tmp_poly, tmp_poly, proof->c); 
    poly_qshow_sub(tmp_poly, Z->entries[PARAM_ARP_DIV_N_SHOW + i], tmp_poly); // Z_{256/n + i} - c.(sum_gamma_z3 + hi) 

    // sum of gamma_{ij}r_j*z1
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      if (j == 0) {
        poly_qshow_vec_m1_mul_scalar(sum_gamma_r_star_i, chal_1[j], chal_2[i][j]);
      } else {
        poly_qshow_vec_m1_mul_scalar(tmp_vec_m1, chal_1[j], chal_2[i][j]);
        poly_qshow_vec_m1_add(sum_gamma_r_star_i, sum_gamma_r_star_i, tmp_vec_m1);
      }
    }
    poly_qshow_vec_m1_mul_inner(z1s_z1, sum_gamma_r_star_i, proof->z1);
    poly_qshow_add(tmp_poly, tmp_poly, z1s_z1);

    // sum of gamma_{ij}e_j*Z = <conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}])), Z_{:256/n}>
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_set_coeff(sum_gamma_e_star_ij, 0, chal_2[i][j * PARAM_N_SHOW + 0]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_e_star_ij, k, - chal_2[i][j * PARAM_N_SHOW + (PARAM_N_SHOW - k)]); // set conjugate directly
      }
      poly_qshow_mul(sum_gamma_e_star_ij, sum_gamma_e_star_ij, Z->entries[j]);
      poly_qshow_add(tmp_poly, tmp_poly, sum_gamma_e_star_ij);
    }
    poly_qshow_mul(tmp_poly, tmp_poly, chal_3_l->entries[i]);
    poly_qshow_add(z1s_c_one, z1s_c_one, tmp_poly);
  }
  poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
  poly_qshow_add(t0, t0, z1s_c_one);

  // computing c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed]z1_v1 - B_embed.z1_v2 + A3_embed.z1_v3 - Ds_embed.z1_usk - D_embed.z1_m - c.u_embed)_i
  // z1s_z1, z1s_c_one can be used as temp variables
  poly_qshow_zero(z1s_c_one); 
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    i_k_quot = i / PARAM_K_SHOW;
    i_k_rem = i % PARAM_K_SHOW;
    // y1s_y1 used as temp variable to host ([A|-B|A3|-Ds|-D].[y1_v1|y1_v2|y1_v3|y1_s|y1_m])_i
    poly_qshow_mul_scalar(z1s_z1, proof->z1->entries[IDX_V1_SHOW + i], PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_V12_SHOW + j]);
      poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
    }
    for (j = 0; j < PARAM_D*PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, B_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_V2_SHOW + j]);
      poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A3_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_V3_SHOW + j]);
      poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
    }
    for (j = 0; j < 2*PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, Ds_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_USK_SHOW + j]);
      poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_M*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, D_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof->z1->entries[IDX_M_SHOW + j]);
      poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly); // substraction
    }
    poly_qshow_mul(tmp_poly, proof->c, u_embed[i_k_quot]->entries[i_k_rem]);
    poly_qshow_sub(z1s_z1, z1s_z1, tmp_poly);
    poly_qshow_mul(z1s_z1, z1s_z1, chal_3_dk[i_k_quot]->entries[i_k_rem]);
    poly_qshow_add(z1s_c_one, z1s_c_one, z1s_z1);
  }
  poly_qshow_mul(z1s_c_one, z1s_c_one, proof->c);
  poly_qshow_add(t0, t0, z1s_c_one);

  // computing fourth challenge
  poly_qshow_zero(c);
  buf[0] = 4;
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  poly_qshow_sample_challenge(c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, proof->ctr_c, SEED_BYTES);

  /****** checks ******/
  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 and square l2 norms
  sq_norm_z1 = poly_qshow_vec_m1_norm2(proof->z1);
  sq_norm_z2 = poly_qshow_vec_m2_norm2(proof->z2);
  sq_norm_z3 = 0;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    tmp_coeff = proof->z3[i];
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp_coeff * tmp_coeff);
  }
  is_valid = is_valid && (sq_norm_z1 <= ((uint128)PARAM_B1SQ_SHOW_LOW64 + (((uint128)PARAM_B1SQ_SHOW_HIGH64) << 64)));
  is_valid = is_valid && (sq_norm_z2 <= PARAM_B2SQ_SHOW) && (sq_norm_z3 <= PARAM_B3SQ_SHOW) && poly_qshow_equal(proof->c, c);

  // clean up
show_verify_cleanup:
  // clean up polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(z1s_z1);
  poly_qshow_clear(z1s_c_one);
  poly_qshow_clear(c);
  poly_qshow_clear(c_one);
  poly_qshow_clear(t0);
  // clean up vectors and matrices
  poly_qshow_mat_d_m1_clear(A1);
  poly_qshow_mat_d_m2_clear(A2);
  poly_qshow_mat_256l_m2_clear(Byg);
  poly_qshow_vec_k_clear(tmp_vec_k);
  poly_qshow_vec_m1_clear(tmp_vec_m1);
  poly_qshow_vec_m1_clear(sum_gamma_r_star_i);
  poly_qshow_vec_d_clear(tmp_vec_d);
  poly_qshow_vec_d_clear(w);
  poly_qshow_clear(sum_gamma_e_star_ij);    
  poly_qshow_vec_m2_clear(b);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_vec_256_l_clear(Z);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_clear(chal_1[i]);
  }
  poly_qshow_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_clear(chal_3_quad_matrix[i]);
    poly_qshow_vec_k_clear(chal_3_dk[i]);
    poly_qshow_vec_k_clear(Gz1_v2[i]);
  }
  for (i = 0; i < 7; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
  }

  return is_valid;
}
