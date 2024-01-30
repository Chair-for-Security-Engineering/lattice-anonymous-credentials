#ifndef SAMPLING__H
#define SAMPLING__H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void poly_q_mat_d_d_uniform(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset);
void poly_q_mat_d_k_uniform(poly_q_mat_d_k mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_q_mat_d_m_uniform(poly_q_mat_d_m mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_q_vec_d_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);

void poly_q_mat_d_d_binomial(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator);

void poly_q_vec_d_bin_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset);

void poly_q_vec_2d_dk_sample_pre(poly_q_vec_d v1[2], poly_q_vec_d v2[PARAM_K], const poly_q_mat_d_d R[2][PARAM_K], const poly_q_mat_d_d A,
                                const poly_q_mat_d_d B[PARAM_K], const poly_q_vec_d u, const poly_q tag, const poly_real_mat_2d_2d S);

void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]);

#endif
