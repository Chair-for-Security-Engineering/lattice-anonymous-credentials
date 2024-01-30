#ifndef SHOW_H
#define SHOW_H

#include <stdint.h>
#include "osig.h"
#include "params.h"

// Byte-length of language packing 
#define SHOW_CHALLENGE_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES) // no need to pack signer public key
// Byte-length of first challenge input
#define CHAL1_SHOW_INPUT_BYTES (SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES + POLYQSHOW_VEC256L_PACKEDBYTES)
// Byte-length of second challenge input 
#define CHAL2_SHOW_INPUT_BYTES (CHAL1_SHOW_INPUT_BYTES + PARAM_ARP_SHOW*COEFFQSHOW_PACKEDBYTES)
// Byte-length of third challenge input
#define CHAL3_SHOW_INPUT_BYTES (CHAL2_SHOW_INPUT_BYTES + POLYQSHOW_VECL_PACKEDBYTES)
// Byte-length of fourth challenge input
#define CHAL4_SHOW_INPUT_BYTES (CHAL3_SHOW_INPUT_BYTES + 2*POLYQSHOW_PACKEDBYTES)

// Starting index for v1 in s1
#define IDX_V1_SHOW 0
// Starting index for v12 in s1
#define IDX_V12_SHOW (IDX_V1_SHOW + PARAM_D*PARAM_K_SHOW)
// Starting index for v2 in s1
#define IDX_V2_SHOW (IDX_V12_SHOW + PARAM_D*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for v3 in s1 
#define IDX_V3_SHOW (IDX_V2_SHOW + PARAM_D*PARAM_K*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for tag in s1 
#define IDX_TAG_SHOW (IDX_V3_SHOW + PARAM_K*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for usk in s1 
#define IDX_USK_SHOW (IDX_TAG_SHOW + PARAM_K_SHOW)
// Starting index for m in s1 
#define IDX_M_SHOW (IDX_USK_SHOW + 2*PARAM_D*PARAM_K_SHOW)

typedef struct {
  poly_qshow_vec_d tA;
  poly_qshow_vec_256_l tB;
  coeff_qshow z3[PARAM_ARP_SHOW];
  poly_qshow_vec_l h;
  poly_qshow t1;
  poly_qshow c;
  uint32_t ctr_c;
  poly_qshow_vec_m1 z1;
  poly_qshow_vec_m2 z2;
} show_proof_t;

void show_proof_init(show_proof_t *proof);
void show_proof_clear(show_proof_t *proof);

void show_user_embed(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    poly_qshow_vec_k   u_embed[PARAM_D], 
    poly_qshow_vec_m1  s1, 
    const user_pk_t    *upk, 
    const user_sk_t    *usk, 
    const sep_pk_t     *pk, 
    const sep_sig_t    *sig,
    const uint8_t      msg[PARAM_M*PARAM_N/8]);
void show_user_prove(
    show_proof_t             *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);
int show_verify(
    const show_proof_t       *proof, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    const poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);

#endif
