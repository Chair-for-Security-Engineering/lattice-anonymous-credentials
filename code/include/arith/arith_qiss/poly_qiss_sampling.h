#ifndef POLY_QISS_SAMPLING_H
#define POLY_QISS_SAMPLING_H

#include <stdint.h>
#include "arith.h"

void vec_qiss_uniform(coeff_qiss out[PARAM_ARP_ISS + 1], const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);
void poly_qiss_uniform_but_zero(poly_qiss out, const uint8_t seed[SEED_BYTES], uint32_t kappa, uint32_t domain_separator);
void poly_qiss_mat_d_k_uniform(poly_qiss_mat_d_k mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset);
void poly_qiss_mat_d_m2_uniform(poly_qiss_mat_d_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_mat_256l_m2_uniform(poly_qiss_mat_256l_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_vec_m2_uniform(poly_qiss_vec_m2 vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_vec_l_uniform(poly_qiss_vec_l vec, const uint8_t *buf, uint32_t domain_separator, size_t buflen);
void poly_qiss_vec_k_uniform(poly_qiss_vec_k vec, const uint8_t *buf, uint32_t domain_separator, uint32_t cnt, size_t buflen);
// TODO streamline the binomial sampling function signatures
void poly_qiss_vec_k_binomial(poly_qiss_vec_k res, const uint8_t *buf, uint32_t domain_separator, uint32_t i, uint32_t j, size_t inlen);
void poly_qiss_vec_m2_binomial(poly_qiss_vec_m2 res, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator);

void poly_qiss_vec_k_sample_gaussian_s1(poly_qiss_vec_k res);
void poly_qiss_vec_m2_sample_gaussian_s2(poly_qiss_vec_m2 res);
void poly_qiss_sample_challenge(poly_qiss out, const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);

#endif
