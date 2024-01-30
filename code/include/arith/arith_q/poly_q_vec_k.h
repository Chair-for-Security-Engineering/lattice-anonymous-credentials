#ifndef POLY_Q_VEC_K_H
#define POLY_Q_VEC_K_H

#include "poly_q.h"

typedef struct {
	poly_q entries[PARAM_K];
} __poly_q_vec_k;

typedef __poly_q_vec_k poly_q_vec_k[1];

void poly_q_vec_k_setup(void);
void poly_q_vec_k_teardown(void);

void poly_q_vec_k_init(poly_q_vec_k arg);
void poly_q_vec_k_clear(poly_q_vec_k arg);
void poly_q_vec_k_zero(poly_q_vec_k arg);
void poly_q_vec_k_set(poly_q_vec_k res, const poly_q_vec_k arg);
void poly_q_vec_k_get_poly(poly_q res, const poly_q_vec_k arg, size_t pos);
void poly_q_vec_k_set_poly(poly_q_vec_k res, const poly_q arg, size_t pos);
void poly_q_vec_k_neg(poly_q_vec_k res, const poly_q_vec_k arg);
void poly_q_vec_k_add(poly_q_vec_k res, const poly_q_vec_k lhs, const poly_q_vec_k rhs);
void poly_q_vec_k_sub(poly_q_vec_k res, const poly_q_vec_k lhs, const poly_q_vec_k rhs);
void poly_q_vec_k_mul_scalar(poly_q_vec_k res, const poly_q_vec_k arg, const poly_q fac);
void poly_q_vec_k_mul_inner(poly_q res, const poly_q_vec_k lhs, const poly_q_vec_k rhs);
uint64_t poly_q_vec_k_norm2(const poly_q_vec_k arg);
void poly_q_vec_k_sample_gaussian_s2(poly_q_vec_k res);
int poly_q_vec_k_equal(const poly_q_vec_k lhs, const poly_q_vec_k rhs);
void poly_q_vec_k_dump(const poly_q_vec_k arg);

#endif /* POLY_Q_VEC_K_H */
