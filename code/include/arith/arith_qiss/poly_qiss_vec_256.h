#ifndef POLY_QISS_VEC_256_H
#define POLY_QISS_VEC_256_H

#include "poly_qiss.h"

typedef struct {
	poly_qiss entries[PARAM_ARP_DIV_N_ISS];
} __poly_qiss_vec_256;

typedef __poly_qiss_vec_256 poly_qiss_vec_256[1];

void poly_qiss_vec_256_setup(void);
void poly_qiss_vec_256_teardown(void);

void poly_qiss_vec_256_init(poly_qiss_vec_256 arg);
void poly_qiss_vec_256_clear(poly_qiss_vec_256 arg);
void poly_qiss_vec_256_zero(poly_qiss_vec_256 arg);
void poly_qiss_vec_256_set(poly_qiss_vec_256 res, const poly_qiss_vec_256 arg);
void poly_qiss_vec_256_get_poly(poly_qiss res, const poly_qiss_vec_256 arg, size_t pos);
void poly_qiss_vec_256_set_poly(poly_qiss_vec_256 res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_256_neg(poly_qiss_vec_256 res, const poly_qiss_vec_256 arg);
void poly_qiss_vec_256_add(poly_qiss_vec_256 res, const poly_qiss_vec_256 lhs, const poly_qiss_vec_256 rhs);
void poly_qiss_vec_256_sub(poly_qiss_vec_256 res, const poly_qiss_vec_256 lhs, const poly_qiss_vec_256 rhs);
void poly_qiss_vec_256_mul_scalar(poly_qiss_vec_256 res, const poly_qiss_vec_256 arg, const coeff_qiss fac);
void poly_qiss_vec_256_mul_inner(poly_qiss res, const poly_qiss_vec_256 lhs, const poly_qiss_vec_256 rhs);
uint64_t poly_qiss_vec_256_norm2(const poly_qiss_vec_256 arg);
int poly_qiss_vec_256_equal(const poly_qiss_vec_256 lhs, const poly_qiss_vec_256 rhs);
void poly_qiss_vec_256_dump(const poly_qiss_vec_256 arg);

#endif /* POLY_QISS_VEC_256_H */
