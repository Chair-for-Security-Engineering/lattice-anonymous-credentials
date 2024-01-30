#ifndef POLY_QISS_VEC_M2_H
#define POLY_QISS_VEC_M2_H

#include "poly_qiss.h"

typedef struct {
	poly_qiss entries[PARAM_M2_ISS];
} __poly_qiss_vec_m2;

typedef __poly_qiss_vec_m2 poly_qiss_vec_m2[1];

void poly_qiss_vec_m2_setup(void);
void poly_qiss_vec_m2_teardown(void);

void poly_qiss_vec_m2_init(poly_qiss_vec_m2 arg);
void poly_qiss_vec_m2_clear(poly_qiss_vec_m2 arg);
void poly_qiss_vec_m2_zero(poly_qiss_vec_m2 arg);
void poly_qiss_vec_m2_set(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 arg);
void poly_qiss_vec_m2_get_poly(poly_qiss res, const poly_qiss_vec_m2 arg, size_t pos);
void poly_qiss_vec_m2_set_poly(poly_qiss_vec_m2 res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_m2_neg(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 arg);
void poly_qiss_vec_m2_add(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 lhs, const poly_qiss_vec_m2 rhs);
void poly_qiss_vec_m2_sub(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 lhs, const poly_qiss_vec_m2 rhs);
void poly_qiss_vec_m2_mul_scalar(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 arg, const coeff_qiss fac);
void poly_qiss_vec_m2_mul_poly_qiss(poly_qiss_vec_m2 res, const poly_qiss_vec_m2 lhs, const poly_qiss rhs);
void poly_qiss_vec_m2_mul_inner(poly_qiss res, const poly_qiss_vec_m2 lhs, const poly_qiss_vec_m2 rhs);
uint64_t poly_qiss_vec_m2_norm2(const poly_qiss_vec_m2 arg);
int poly_qiss_vec_m2_equal(const poly_qiss_vec_m2 lhs, const poly_qiss_vec_m2 rhs);
void poly_qiss_vec_m2_dump(const poly_qiss_vec_m2 arg);

#endif /* POLY_QISS_VEC_M2_H */
