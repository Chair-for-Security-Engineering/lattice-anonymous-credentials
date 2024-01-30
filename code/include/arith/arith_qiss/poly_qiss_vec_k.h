#ifndef POLY_QISS_VEC_K_H
#define POLY_QISS_VEC_K_H

#include "poly_qiss.h"

#define POLYQISS_VECK_PACKEDBYTES (PARAM_K_ISS * POLYQISS_PACKEDBYTES)

typedef struct {
	poly_qiss entries[PARAM_K_ISS];
} __poly_qiss_vec_k;

typedef __poly_qiss_vec_k poly_qiss_vec_k[1];

void poly_qiss_vec_k_setup(void);
void poly_qiss_vec_k_teardown(void);

void poly_qiss_vec_k_init(poly_qiss_vec_k arg);
void poly_qiss_vec_k_clear(poly_qiss_vec_k arg);
void poly_qiss_vec_k_zero(poly_qiss_vec_k arg);
void poly_qiss_vec_k_set(poly_qiss_vec_k res, const poly_qiss_vec_k arg);
void poly_qiss_vec_k_get_poly(poly_qiss res, const poly_qiss_vec_k arg, size_t pos);
void poly_qiss_vec_k_set_poly(poly_qiss_vec_k res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_k_neg(poly_qiss_vec_k res, const poly_qiss_vec_k arg);
void poly_qiss_vec_k_add(poly_qiss_vec_k res, const poly_qiss_vec_k lhs, const poly_qiss_vec_k rhs);
void poly_qiss_vec_k_sub(poly_qiss_vec_k res, const poly_qiss_vec_k lhs, const poly_qiss_vec_k rhs);
void poly_qiss_vec_k_mul_scalar(poly_qiss_vec_k res, const poly_qiss_vec_k arg, const coeff_qiss fac);
void poly_qiss_vec_k_mul_poly_qiss(poly_qiss_vec_k res, const poly_qiss_vec_k lhs, const poly_qiss rhs);
void poly_qiss_vec_k_mul_inner(poly_qiss res, const poly_qiss_vec_k lhs, const poly_qiss_vec_k rhs);
void poly_qiss_vec_k_conjugate(poly_qiss_vec_k res, const poly_qiss_vec_k arg);
uint64_t poly_qiss_vec_k_norm2(const poly_qiss_vec_k arg);
int poly_qiss_vec_k_equal(const poly_qiss_vec_k lhs, const poly_qiss_vec_k rhs);
void poly_qiss_vec_k_dump(const poly_qiss_vec_k arg);
void poly_qiss_vec_k_pack(uint8_t buf[POLYQISS_VECK_PACKEDBYTES], const poly_qiss_vec_k arg);

#endif /* POLY_QISS_VEC_K_H */
