#ifndef POLY_Q_MAT_D_K_H
#define POLY_Q_MAT_D_K_H

#include "arith_q.h"

typedef struct {
    poly_q_vec_k rows[PARAM_D];
} __poly_q_mat_d_k;

typedef __poly_q_mat_d_k poly_q_mat_d_k[1];

void poly_q_mat_d_k_setup(void);
void poly_q_mat_d_k_teardown(void);

void poly_q_mat_d_k_init(poly_q_mat_d_k res);
void poly_q_mat_d_k_clear(poly_q_mat_d_k res);
void poly_q_mat_d_k_zero(poly_q_mat_d_k res);
void poly_q_mat_d_k_set(poly_q_mat_d_k res, const poly_q_mat_d_k arg);
void poly_q_mat_d_k_neg(poly_q_mat_d_k res, const poly_q_mat_d_k arg);
void poly_q_mat_d_k_add(poly_q_mat_d_k res, const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs);
void poly_q_mat_d_k_sub(poly_q_mat_d_k res, const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs);
void poly_q_mat_d_k_mul_vec_k(poly_q_vec_d res, const poly_q_mat_d_k lhs, const poly_q_vec_k rhs);
int poly_q_mat_d_k_equal(const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs);
void poly_q_mat_d_k_dump(const poly_q_mat_d_k arg);

#endif /* POLY_Q_MAT_D_K_H */
