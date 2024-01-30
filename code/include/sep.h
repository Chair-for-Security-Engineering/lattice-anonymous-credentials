#ifndef SEP_H
#define SEP_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

typedef struct {
  poly_q_mat_d_d R[2][PARAM_K];
  poly_real_mat_2d_2d S;
} sep_sk_t;

typedef struct {
  poly_q_mat_d_d B[PARAM_K];
  uint8_t seed[SEED_BYTES];
} sep_pk_t;

typedef struct {
  poly_q tag;
  poly_q_vec_d v12;
  poly_q_vec_d v2[PARAM_K];
  poly_q_vec_k v3;
} sep_sig_t;

void sep_keys_init(sep_pk_t *pk, sep_sk_t *sk);
void sep_keys_clear(sep_pk_t *pk, sep_sk_t *sk);
void sep_sig_init(sep_sig_t *sig);
void sep_sig_clear(sep_sig_t *sig);

void sep_keygen(sep_pk_t *pk, sep_sk_t *sk);
void _sep_sign_commitment(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const poly_q_vec_d cmt);
void sep_sign(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const uint8_t msg[PARAM_M*PARAM_N/8]);
int _sep_verify_from_commitment(const sep_sig_t *sig, const poly_q_vec_d cmt, const sep_pk_t *pk);
int sep_verify(const sep_sig_t *sig, const uint8_t msg[PARAM_M*PARAM_N/8], const sep_pk_t *pk);

#endif
