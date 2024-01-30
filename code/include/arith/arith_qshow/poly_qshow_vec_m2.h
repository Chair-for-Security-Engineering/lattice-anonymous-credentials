#ifndef POLY_QSHOW_VEC_M2_H
#define POLY_QSHOW_VEC_M2_H

#include "poly_qshow.h"

typedef struct {
	poly_qshow entries[PARAM_M2_SHOW];
} __poly_qshow_vec_m2;

typedef __poly_qshow_vec_m2 poly_qshow_vec_m2[1];

void poly_qshow_vec_m2_setup(void);
void poly_qshow_vec_m2_teardown(void);

void poly_qshow_vec_m2_init(poly_qshow_vec_m2 arg);
void poly_qshow_vec_m2_clear(poly_qshow_vec_m2 arg);
void poly_qshow_vec_m2_zero(poly_qshow_vec_m2 arg);
void poly_qshow_vec_m2_set(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 arg);
void poly_qshow_vec_m2_get_poly(poly_qshow res, const poly_qshow_vec_m2 arg, size_t pos);
void poly_qshow_vec_m2_set_poly(poly_qshow_vec_m2 res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_m2_neg(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 arg);
void poly_qshow_vec_m2_add(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 lhs, const poly_qshow_vec_m2 rhs);
void poly_qshow_vec_m2_sub(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 lhs, const poly_qshow_vec_m2 rhs);
void poly_qshow_vec_m2_mul_scalar(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 arg, const coeff_qshow fac);
void poly_qshow_vec_m2_mul_poly_qshow(poly_qshow_vec_m2 res, const poly_qshow_vec_m2 lhs, const poly_qshow rhs);
void poly_qshow_vec_m2_mul_inner(poly_qshow res, const poly_qshow_vec_m2 lhs, const poly_qshow_vec_m2 rhs);
uint64_t poly_qshow_vec_m2_norm2(const poly_qshow_vec_m2 arg);
int poly_qshow_vec_m2_equal(const poly_qshow_vec_m2 lhs, const poly_qshow_vec_m2 rhs);
void poly_qshow_vec_m2_dump(const poly_qshow_vec_m2 arg);

#endif /* POLY_QSHOW_VEC_M2_H */
