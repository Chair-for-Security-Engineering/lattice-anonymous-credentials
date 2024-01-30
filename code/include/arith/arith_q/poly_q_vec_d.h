#ifndef POLY_Q_VEC_D_H
#define POLY_Q_VEC_D_H

#include "poly_q.h"

typedef struct {
	poly_q entries[PARAM_D];
} __poly_q_vec_d;

typedef __poly_q_vec_d poly_q_vec_d[1];

void poly_q_vec_d_setup(void);
void poly_q_vec_d_teardown(void);

void poly_q_vec_d_init(poly_q_vec_d arg);
void poly_q_vec_d_clear(poly_q_vec_d arg);
void poly_q_vec_d_zero(poly_q_vec_d arg);
void poly_q_vec_d_set(poly_q_vec_d res, const poly_q_vec_d arg);
void poly_q_vec_d_get_poly(poly_q res, const poly_q_vec_d arg, size_t pos);
void poly_q_vec_d_set_poly(poly_q_vec_d res, const poly_q arg, size_t pos);
void poly_q_vec_d_neg(poly_q_vec_d res, const poly_q_vec_d arg);
void poly_q_vec_d_add(poly_q_vec_d res, const poly_q_vec_d lhs, const poly_q_vec_d rhs);
void poly_q_vec_d_sub(poly_q_vec_d res, const poly_q_vec_d lhs, const poly_q_vec_d rhs);
void poly_q_vec_d_mul_scalar(poly_q_vec_d res, const poly_q_vec_d arg, coeff_q fac);
void poly_q_vec_d_mul_poly(poly_q_vec_d res, const poly_q_vec_d arg, const poly_q arg2);
void poly_q_vec_d_mul_inner(poly_q res, const poly_q_vec_d lhs, const poly_q_vec_d rhs);
uint64_t poly_q_vec_d_norm2(const poly_q_vec_d arg);
void poly_q_vec_d_gaussian_sqrt_s2sq_sGsq(poly_q_vec_d res);
int poly_q_vec_d_equal(const poly_q_vec_d lhs, const poly_q_vec_d rhs);
void poly_q_vec_d_dump(const poly_q_vec_d arg);

#endif /* POLY_Q_VEC_D_H */
