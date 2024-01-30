#ifndef OSIG_H
#define OSIG_H

#include <stdint.h>
#include "params.h"

// Byte-length of language packing 
#define ISS_CHALLENGE_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES + 2*PARAM_D*POLYQISS_VECK_PACKEDBYTES)
// Byte-length of first challenge input
#define CHAL1_ISS_INPUT_BYTES (ISS_CHALLENGE_BASE_BYTES + 2*POLYQISS_VECD_PACKEDBYTES + POLYQISS_VEC256L_PACKEDBYTES)
// Byte-length of second challenge input
#define CHAL2_ISS_INPUT_BYTES (CHAL1_ISS_INPUT_BYTES + PARAM_ARP_ISS*COEFFQISS_PACKEDBYTES)
// Byte-length of third challenge input 
#define CHAL3_ISS_INPUT_BYTES (CHAL2_ISS_INPUT_BYTES + POLYQISS_VECL_PACKEDBYTES)
// Byte-length of fourth challenge input 
#define CHAL4_ISS_INPUT_BYTES (CHAL3_ISS_INPUT_BYTES + 2*POLYQISS_PACKEDBYTES)

typedef struct {
  poly_q_vec_d s[2];
} user_sk_t;

typedef struct {
  poly_q_vec_d t;
  uint8_t seed[SEED_BYTES];
} user_pk_t;

typedef struct {
  poly_qiss_vec_d tA;
  poly_qiss_vec_256_l tB;
  coeff_qiss z3[PARAM_ARP_ISS];
  poly_qiss_vec_l h;
  poly_qiss t1;
  poly_qiss c;
  uint32_t ctr_c;
  poly_qiss_vec_k z1[PARAM_M1_K_ISS];
  poly_qiss_vec_m2 z2;
} osig_proof_t;

void user_keys_init(user_pk_t *upk, user_sk_t *usk);
void user_keys_clear(user_pk_t *upk, user_sk_t *usk);
void osig_proof_init(osig_proof_t *proof);
void osig_proof_clear(osig_proof_t *proof);

void osig_user_keygen(user_pk_t *upk, user_sk_t *usk, const uint8_t seed[SEED_BYTES]);
void osig_user_commit(poly_q_vec_d r[2], poly_q_vec_d cmt, const uint8_t msg[PARAM_M*PARAM_N/8], const user_pk_t *upk);
void osig_user_embed(
    poly_qiss_mat_k_k  A_embed[PARAM_D][PARAM_D], 
    poly_qiss_mat_k_k  Ds_embed[PARAM_D][2*PARAM_D], 
    poly_qiss_mat_k_k  D_embed[PARAM_D][PARAM_M], 
    poly_qiss_vec_k    u[2*PARAM_D], 
    poly_qiss_vec_k    s1[PARAM_M1_K_ISS], 
    const user_pk_t    *upk, 
    const user_sk_t    *usk, 
    const poly_q_vec_d cmt, 
    const poly_q_vec_d r[2], 
    const uint8_t      msg[PARAM_M*PARAM_N/8]);
void osig_user_prove(
    osig_proof_t            *proof, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qiss_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qiss_vec_k   u[2*PARAM_D], 
    const poly_qiss_vec_k   s1[PARAM_M1_K_ISS], 
    const uint8_t           crs_seed[CRS_SEED_BYTES], 
    const uint8_t           seed[SEED_BYTES]);
int osig_signer_verify(
    const osig_proof_t      *proof, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qiss_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qiss_vec_k   u[2*PARAM_D], 
    const uint8_t           crs_seed[CRS_SEED_BYTES], 
    const uint8_t           seed[SEED_BYTES]);
void osig_signer_sign_commitment(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const poly_q_vec_d cmt);
void osig_user_sig_complete(sep_sig_t *sig, const poly_q_vec_d r[2]);
int osig_user_verify(const sep_sig_t *sig, const sep_pk_t *pk, const user_pk_t *upk, const uint8_t msg[PARAM_M*PARAM_N/8]);

#endif
