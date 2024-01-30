#ifndef POLY_QSHOW_VEC_M1_H
#define POLY_QSHOW_VEC_M1_H

#include "poly_qshow.h"

typedef struct {
	poly_qshow entries[PARAM_M1_SHOW];
} __poly_qshow_vec_m1;

typedef __poly_qshow_vec_m1 poly_qshow_vec_m1[1];

void poly_qshow_vec_m1_setup(void);
void poly_qshow_vec_m1_teardown(void);

void poly_qshow_vec_m1_init(poly_qshow_vec_m1 arg);
void poly_qshow_vec_m1_clear(poly_qshow_vec_m1 arg);
void poly_qshow_vec_m1_zero(poly_qshow_vec_m1 arg);
void poly_qshow_vec_m1_set(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 arg);
void poly_qshow_vec_m1_get_poly(poly_qshow res, const poly_qshow_vec_m1 arg, size_t pos);
void poly_qshow_vec_m1_set_poly(poly_qshow_vec_m1 res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_m1_neg(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 arg);
void poly_qshow_vec_m1_add(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 lhs, const poly_qshow_vec_m1 rhs);
void poly_qshow_vec_m1_sub(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 lhs, const poly_qshow_vec_m1 rhs);
void poly_qshow_vec_m1_mul_scalar(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 arg, const coeff_qshow fac);
void poly_qshow_vec_m1_mul_poly_qshow(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 lhs, const poly_qshow rhs);
void poly_qshow_vec_m1_mul_inner(poly_qshow res, const poly_qshow_vec_m1 lhs, const poly_qshow_vec_m1 rhs);
void poly_qshow_vec_m1_conjugate(poly_qshow_vec_m1 res, const poly_qshow_vec_m1 arg);
uint128 poly_qshow_vec_m1_norm2(const poly_qshow_vec_m1 arg);
int poly_qshow_vec_m1_equal(const poly_qshow_vec_m1 lhs, const poly_qshow_vec_m1 rhs);
void poly_qshow_vec_m1_dump(const poly_qshow_vec_m1 arg);

#endif /* POLY_QSHOW_VEC_M1_H */
