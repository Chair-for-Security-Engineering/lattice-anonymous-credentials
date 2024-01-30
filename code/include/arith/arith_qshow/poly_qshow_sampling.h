#ifndef POLY_QSHOW_SAMPLING_H
#define POLY_QSHOW_SAMPLING_H

#include <stdint.h>
#include "arith.h"

void vec_qshow_uniform(coeff_qshow out[PARAM_ARP_SHOW + 6], const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);
void poly_qshow_uniform_but_zero(poly_qshow out, const uint8_t seed[SEED_BYTES], uint32_t kappa, uint32_t domain_separator);
void poly_qshow_mat_d_m1_uniform(poly_qshow_mat_d_m1 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qshow_mat_d_m2_uniform(poly_qshow_mat_d_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qshow_mat_256l_m2_uniform(poly_qshow_mat_256l_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qshow_vec_m2_uniform(poly_qshow_vec_m2 vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qshow_vec_l_uniform(poly_qshow_vec_l vec, const uint8_t *buf, uint32_t domain_separator, size_t buflen);
void poly_qshow_vec_k_uniform(poly_qshow_vec_k vec, const uint8_t *buf, uint32_t domain_separator, uint32_t cnt, size_t buflen);
// TODO streamline the binomial sampling function signatures
void poly_qshow_vec_m1_binomial(poly_qshow_vec_m1 res, const uint8_t *buf, uint32_t domain_separator, uint32_t i, size_t inlen);
void poly_qshow_vec_m2_binomial(poly_qshow_vec_m2 res, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator);

void poly_qshow_vec_m1_sample_gaussian_s1(poly_qshow_vec_m1 res);
void poly_qshow_vec_m2_sample_gaussian_s2(poly_qshow_vec_m2 res);
void poly_qshow_sample_challenge(poly_qshow out, const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);

#endif
