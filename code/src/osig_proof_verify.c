#include "arith.h"
#include "randombytes.h"
#include "poly_qiss_sampling.h"
#include "sep.h"
#include "osig.h"
#include "macros.h"
#include "fips202.h"

/*************************************************
* Name:        osig_signer_verify
*
* Description: Signer verification of the zero-knowledge proof of commitment 
*              opening and user registration (anonymous credentials issuance)
*
* Arguments:   - const osig_proof *proof: pointer to issuance proof structure
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *Ds_embed: array of polynomial matrices, subring embedding of q1.Ds
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices, subring embedding of q1.D
*              - const poly_qiss_vec_k *u: array of polynomial vectors, subring embedding of q1.(cmt-upk | upk)
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
* 
* Returns 1 if proof could be verified correctly and 0 otherwise
**************************************************/
int osig_signer_verify(
    const osig_proof_t      *proof, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qiss_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qiss_vec_k   u[2*PARAM_D], 
    const uint8_t           crs_seed[CRS_SEED_BYTES], 
    const uint8_t           seed[SEED_BYTES]) {
  size_t i,j,k;
  int is_valid = 1;
  uint8_t buf[CHAL4_ISS_INPUT_BYTES] = {0}, challenge_seed[SEED_BYTES];
  uint64_t sq_norm_z3, sq_norm_z2, sq_norm_z1;
  coeff_qiss chal_2[PARAM_L_ISS][PARAM_ARP_ISS + 1];
  coeff_qiss tmp_coeff;
  poly_qiss tmp_poly, z1_star_z1, z1_star_one, sum_mu_gamma, t0, c;
  poly_qiss_mat_d_k A1[PARAM_M1_K_ISS];
  poly_qiss_mat_d_m2 A2;
  poly_qiss_mat_256l_m2 Byg;
  poly_qiss_vec_d tmp_vec_d, w;
  poly_qiss_vec_k tmp_vec_k, tmp_C_z1, one, sum_gamma_r_star_ik;
  poly_qiss_vec_k chal_3_2dk[2*PARAM_D];
  poly_qiss sum_gamma_e_star_ij;
  poly_qiss_vec_k chal_1[PARAM_ARP_ISS][PARAM_M1_K_ISS];
  poly_qiss_vec_m2 b;
  poly_qiss_vec_256_l tmp_vec_256_l, Z;
  poly_qiss_vec_l chal_3_l;

  // init
  // init polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_init(z1_star_z1);
  poly_qiss_init(z1_star_one);
  poly_qiss_init(c);
  poly_qiss_init(sum_mu_gamma);
  poly_qiss_init(t0);
  // init vectors and matrices
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_init(A1[i]);
  }
  poly_qiss_mat_d_m2_init(A2);
  poly_qiss_mat_256l_m2_init(Byg);
  poly_qiss_vec_k_init(tmp_vec_k);
  poly_qiss_vec_k_init(tmp_C_z1);
  poly_qiss_vec_k_init(one);
  poly_qiss_vec_k_init(sum_gamma_r_star_ik);
  poly_qiss_vec_d_init(tmp_vec_d);
  poly_qiss_vec_d_init(w);
  poly_qiss_init(sum_gamma_e_star_ij);    
  poly_qiss_vec_m2_init(b);
  poly_qiss_vec_256_l_init(tmp_vec_256_l);
  poly_qiss_vec_256_l_init(Z);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_init(chal_1[i][j]);
    }
  }
  poly_qiss_vec_l_init(chal_3_l);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(chal_3_2dk[i]);
  }

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

  // precomputation
  for (j = 0; j < PARAM_K_ISS; j++) {
    for (k = 0; k < PARAM_N_ISS; k++) {
      poly_qiss_set_coeff(one->entries[j], k, 1); // poly with all one
    }
  }
  ////////////////////////////////////// NO INITIALIZATIONS BELOW THIS LINE

  /****** reconstructing first message ******/
  // recomputing w = A1.z1 + A2.z2 - c.tA
  poly_qiss_vec_d_mul_poly_qiss(tmp_vec_d, proof->tA, proof->c);
  poly_qiss_mat_d_m2_mul_vec_m2(w, A2, proof->z2);
  poly_qiss_vec_d_sub(w, w, tmp_vec_d);
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_mul_vec_k(tmp_vec_d, A1[i], proof->z1[i]); // A1.z1 
    poly_qiss_vec_d_add(w, w, tmp_vec_d);
  }
  // recomputing first challenge
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

  /****** reconstructing second message ******/
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    coeff_qiss_pack(buf + CHAL1_ISS_INPUT_BYTES + i*COEFFQISS_PACKEDBYTES, proof->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_ISS_INPUT_BYTES);
  // recomputing second challenge
  for (i = 0; i < PARAM_L_ISS; i++) {
    vec_qiss_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_ISS, i, SEED_BYTES); // writes PARAM_ARP_ISS + 1 uniform numbers to the first argument
  }

  /****** reconstructing third message ******/
  buf[0] = 3;
  poly_qiss_vec_l_pack(buf + CHAL2_ISS_INPUT_BYTES, proof->h);
  // recomputing third challenge
  poly_qiss_vec_l_uniform(chal_3_l, buf, DOMAIN_SEPARATOR_CHAL3_ISS, CHAL3_ISS_INPUT_BYTES);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_ISS_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_uniform(chal_3_2dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, i+PARAM_L_ISS, SEED_BYTES);
  }
  
  /****** reconstructing fourth message ******/
  /**********************************************
  * Recomputing commitment t0
  *   Z = c.tB - Byg.z2
  *   
  *   t0 = zFz + c.<f,z> + c^2.f + <b,z2> - c.t1
  *      = <b, z2> - c.t1 + sum_mu_gamma.(<z1*, z1> - c.<z1*, one>)
  *        + c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n:}
  *        + <sum_gamma_e_star[i], Z_{:256/n} - c.(sum_gamma_z3 + hi))
  *        + c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed]z1_{:2dk} + D_embed.z1_{4dk:} - c.u_{:dk})_i
  *        + c.sum_{i < dk} mu_{l+dk+i}.(Ds_embed.y1_{2dk:4dk} - c.u_{dk:})_i
  **********************************************/
  // computing <b, z2> - c.t1
  poly_qiss_vec_m2_mul_inner(t0, b, proof->z2);
  poly_qiss_mul(tmp_poly, proof->c, proof->t1);
  poly_qiss_sub(t0, t0, tmp_poly);

  // computing sum_i mu_i gamma_{i,257}
  poly_qiss_mul_scalar(sum_mu_gamma, chal_3_l->entries[0], chal_2[0][PARAM_ARP_ISS]);
  for (i = 1; i < PARAM_L_ISS; i++) {
    poly_qiss_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[i][PARAM_ARP_ISS]);
    poly_qiss_add(sum_mu_gamma, sum_mu_gamma, tmp_poly);
  }

  // computing sum_mu_gamma.(<z1*, z1> - c.<z1*, one>)
  poly_qiss_zero(z1_star_z1);
  poly_qiss_zero(z1_star_one);
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_conjugate(tmp_vec_k, proof->z1[i]); // z1*
    poly_qiss_vec_k_mul_inner(tmp_poly, tmp_vec_k, proof->z1[i]);
    poly_qiss_add(z1_star_z1, z1_star_z1, tmp_poly); // <z1*, z1>
    poly_qiss_vec_k_mul_inner(tmp_poly, tmp_vec_k, one);
    poly_qiss_add(z1_star_one, z1_star_one, tmp_poly); // <z1*, one>
  }
  poly_qiss_mul(z1_star_one, z1_star_one, proof->c);
  poly_qiss_sub(tmp_poly, z1_star_z1, z1_star_one);
  poly_qiss_mul(tmp_poly, tmp_poly, sum_mu_gamma);
  poly_qiss_add(t0, t0, tmp_poly);

  // computing Z = c.tB - Byg.z2
  poly_qiss_mat_256l_m2_mul_vec_m2(tmp_vec_256_l, Byg, proof->z2);
  poly_qiss_vec_256_l_mul_poly_qiss(Z, proof->tB, proof->c);
  poly_qiss_vec_256_l_sub(Z, Z, tmp_vec_256_l);

  // computing c.sum_{i < l} mu_i.(<sum_gamma_r_star[i], z1> + Z_{256/n + i} + <sum_gamma_e_star[i], Z_{:256/n}> - c.(sum_gamma_z3 + hi))
  for (i = 0; i < PARAM_L_ISS; i++) {
    // z1_star_one can be used as a temporary variable
    poly_qiss_set(tmp_poly, proof->h->entries[i]); // h_i
    is_valid = is_valid && (poly_qiss_get_coeff(tmp_poly, 0) == 0);
    if (!is_valid) {
      goto osig_verify_cleanup; // h_i does not have a zero coefficient, is_valid = 0
    }
    tmp_coeff = 0; // h_i constant coefficient is 0
    // sum_j gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_ISS; j++) {
      tmp_coeff += chal_2[i][j] * proof->z3[j];
    }
    poly_qiss_set_coeff(tmp_poly, 0, tmp_coeff);
    poly_qiss_mul(tmp_poly, tmp_poly, proof->c); 
    poly_qiss_sub(tmp_poly, Z->entries[PARAM_ARP_DIV_N_ISS + i], tmp_poly); // Z_{256/n + i} - c.(sum_gamma_z3 + hi) 

    // sum of gamma_{ij}r_j*z1
    for (k = 0; k < PARAM_M1_K_ISS; k++) {
      poly_qiss_vec_k_conjugate(sum_gamma_r_star_ik, chal_1[0][k]);
      poly_qiss_vec_k_mul_scalar(sum_gamma_r_star_ik, sum_gamma_r_star_ik, chal_2[i][0]);
      for (j = 1; j < PARAM_ARP_ISS; j++) {
        poly_qiss_vec_k_conjugate(tmp_vec_k, chal_1[j][k]);
        poly_qiss_vec_k_mul_scalar(tmp_vec_k, tmp_vec_k, chal_2[i][j]);
        poly_qiss_vec_k_add(sum_gamma_r_star_ik, sum_gamma_r_star_ik, tmp_vec_k);
      }
      poly_qiss_vec_k_mul_inner(z1_star_one, sum_gamma_r_star_ik, proof->z1[k]);
      poly_qiss_add(tmp_poly, tmp_poly, z1_star_one);
    }

    // sum of gamma_{ij}e_j*Z = <conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}])), Z_{:256/n}>
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      poly_qiss_set_coeff(sum_gamma_e_star_ij, 0, chal_2[i][j * PARAM_N_ISS + 0]);
      for (k = 1; k < PARAM_N_ISS; k++) {
        poly_qiss_set_coeff(sum_gamma_e_star_ij, k, - chal_2[i][j * PARAM_N_ISS + (PARAM_N_ISS - k)]); // set conjugate directly
      }
      poly_qiss_mul(sum_gamma_e_star_ij, sum_gamma_e_star_ij, Z->entries[j]);
      poly_qiss_add(tmp_poly, tmp_poly, sum_gamma_e_star_ij);
    }
    poly_qiss_mul(tmp_poly, tmp_poly, chal_3_l->entries[i]);
    poly_qiss_mul(tmp_poly, tmp_poly, proof->c);
    poly_qiss_add(t0, t0, tmp_poly);
  }

  // computing c.sum_{i < dk} mu_{l+i}.([q1.I | A_embed]z1_{:2dk} + D_embed.z1_{4dk:} - c.u_{:dk})_i
  // + c.sum_{i < dk} mu_{l+dk+i}.(Ds_embed.z1_{2dk:4dk} - c.u_{dk:})_i
  // z1_star_one can be used as a temporary variable
  poly_qiss_zero(z1_star_one);
  for (i = 0; i < PARAM_D; i++) {
    // Top part of C.z1
    poly_qiss_vec_k_mul_scalar(tmp_C_z1, proof->z1[i], PARAM_Q1_ISS); // q_1*I_{dk} x z1[:dk]
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, A_embed[i][j], proof->z1[PARAM_D + j]); // part of z1 corresponding to r_{12}
      poly_qiss_vec_k_add(tmp_C_z1, tmp_C_z1, tmp_vec_k);
    }
    for (j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, D_embed[i][j], proof->z1[4*PARAM_D + j]); // part of z1 corresponding to m
      poly_qiss_vec_k_add(tmp_C_z1, tmp_C_z1, tmp_vec_k);
    }
    // adding mu_{l+i} * ( [[q_1*I | A_embed | 0 | D_embed]z_1]_i - c.u_i)
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(tmp_poly, u[i]->entries[j], proof->c);
      poly_qiss_sub(tmp_poly, tmp_C_z1->entries[j], tmp_poly);
      poly_qiss_mul(tmp_poly, tmp_poly, chal_3_2dk[i]->entries[j]);
      poly_qiss_add(z1_star_one, z1_star_one, tmp_poly);
    }

    // Bottom part of C.z1
    poly_qiss_mat_k_k_mul_vec_k(tmp_C_z1, Ds_embed[i][0], proof->z1[2*PARAM_D + 0]);
    for (j = 1; j < 2*PARAM_D; j++) {
      poly_qiss_mat_k_k_mul_vec_k(tmp_vec_k, Ds_embed[i][j], proof->z1[2*PARAM_D + j]); // part of z1 corresponding to usk
      poly_qiss_vec_k_add(tmp_C_z1, tmp_C_z1, tmp_vec_k);
    }
    // adding mu_{l+dk+i} * ([[0 | Ds_embed | 0]z_1]_i - c.u_i)
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(tmp_poly, u[PARAM_D + i]->entries[j], proof->c);
      poly_qiss_sub(tmp_poly, tmp_C_z1->entries[j], tmp_poly);
      poly_qiss_mul(tmp_poly, tmp_poly, chal_3_2dk[PARAM_D + i]->entries[j]);
      poly_qiss_add(z1_star_one, z1_star_one, tmp_poly);
    }
  }
  poly_qiss_mul(z1_star_one, z1_star_one, proof->c);
  poly_qiss_add(t0, t0, z1_star_one);

  // computing fourth challenge
  buf[0] = 4;
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES, t0);
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES + POLYQISS_PACKEDBYTES, proof->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_ISS_INPUT_BYTES);
  poly_qiss_sample_challenge(c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_ISS, proof->ctr_c, SEED_BYTES);

  /****** checks ******/
  // computing z1 = y1 + c.s1, z2 = y2 + c.s2 and square l2 norms
  sq_norm_z1 = 0;
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    CHK_UI_OVF_ADDITION(sq_norm_z1, poly_qiss_vec_k_norm2(proof->z1[i]));
  }
  sq_norm_z2 = poly_qiss_vec_m2_norm2(proof->z2);
  sq_norm_z3 = 0;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    tmp_coeff = proof->z3[i];
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp_coeff * tmp_coeff);
  }

  is_valid = is_valid && (sq_norm_z1 <= PARAM_B1SQ_ISS) && (sq_norm_z2 <= PARAM_B2SQ_ISS) && (sq_norm_z3 <= PARAM_B3SQ_ISS) && poly_qiss_equal(proof->c, c);

osig_verify_cleanup:
  // clean up
  // clear polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(z1_star_z1);
  poly_qiss_clear(z1_star_one);
  poly_qiss_clear(c);
  poly_qiss_clear(sum_mu_gamma);
  poly_qiss_clear(t0);
  // clear vectors and matrices
  for (i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_mat_d_k_clear(A1[i]);
  }
  poly_qiss_mat_d_m2_clear(A2);
  poly_qiss_mat_256l_m2_clear(Byg);
  poly_qiss_vec_k_clear(tmp_vec_k);
  poly_qiss_vec_k_clear(tmp_C_z1);
  poly_qiss_vec_k_clear(one);
  poly_qiss_vec_k_clear(sum_gamma_r_star_ik);
  poly_qiss_vec_d_clear(tmp_vec_d);
  poly_qiss_vec_d_clear(w);
  poly_qiss_clear(sum_gamma_e_star_ij);
  poly_qiss_vec_m2_clear(b);
  poly_qiss_vec_256_l_clear(tmp_vec_256_l);
  poly_qiss_vec_256_l_clear(Z);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    for (j = 0; j < PARAM_M1_K_ISS; j++) {
      poly_qiss_vec_k_clear(chal_1[i][j]);
    }
  }
  poly_qiss_vec_l_clear(chal_3_l);
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(chal_3_2dk[i]);
  }

  return is_valid;
}
