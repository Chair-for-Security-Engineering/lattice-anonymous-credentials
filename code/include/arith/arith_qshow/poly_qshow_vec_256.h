#ifndef POLY_QSHOW_VEC_256_H
#define POLY_QSHOW_VEC_256_H

#include "poly_qshow.h"

typedef struct {
	poly_qshow entries[PARAM_ARP_DIV_N_SHOW];
} __poly_qshow_vec_256;

typedef __poly_qshow_vec_256 poly_qshow_vec_256[1];

void poly_qshow_vec_256_setup(void);
void poly_qshow_vec_256_teardown(void);

void poly_qshow_vec_256_init(poly_qshow_vec_256 arg);
void poly_qshow_vec_256_clear(poly_qshow_vec_256 arg);
void poly_qshow_vec_256_zero(poly_qshow_vec_256 arg);
void poly_qshow_vec_256_set(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg);
void poly_qshow_vec_256_get_poly(poly_qshow res, const poly_qshow_vec_256 arg, size_t pos);
void poly_qshow_vec_256_set_poly(poly_qshow_vec_256 res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_256_neg(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg);
void poly_qshow_vec_256_add(poly_qshow_vec_256 res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs);
void poly_qshow_vec_256_sub(poly_qshow_vec_256 res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs);
void poly_qshow_vec_256_mul_scalar(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg, const coeff_qshow fac);
void poly_qshow_vec_256_mul_inner(poly_qshow res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs);
uint64_t poly_qshow_vec_256_norm2(const poly_qshow_vec_256 arg);
int poly_qshow_vec_256_equal(const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs);
void poly_qshow_vec_256_dump(const poly_qshow_vec_256 arg);

#endif /* POLY_QSHOW_VEC_256_H */
