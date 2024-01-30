#ifndef POLY_QSHOW_MAT_D_K_H
#define POLY_QSHOW_MAT_D_K_H

#include "arith_qshow.h"

typedef struct {
    poly_qshow_vec_k rows[PARAM_D_SHOW];
} __poly_qshow_mat_d_k;

typedef __poly_qshow_mat_d_k poly_qshow_mat_d_k[1];

void poly_qshow_mat_d_k_setup(void);
void poly_qshow_mat_d_k_teardown(void);

void poly_qshow_mat_d_k_init(poly_qshow_mat_d_k res);
void poly_qshow_mat_d_k_clear(poly_qshow_mat_d_k res);
void poly_qshow_mat_d_k_zero(poly_qshow_mat_d_k res);
void poly_qshow_mat_d_k_set(poly_qshow_mat_d_k res, const poly_qshow_mat_d_k arg);
void poly_qshow_mat_d_k_neg(poly_qshow_mat_d_k res, const poly_qshow_mat_d_k arg);
void poly_qshow_mat_d_k_add(poly_qshow_mat_d_k res, const poly_qshow_mat_d_k lhs, const poly_qshow_mat_d_k rhs);
void poly_qshow_mat_d_k_sub(poly_qshow_mat_d_k res, const poly_qshow_mat_d_k lhs, const poly_qshow_mat_d_k rhs);
void poly_qshow_mat_d_k_mul_vec_k(poly_qshow_vec_d res, const poly_qshow_mat_d_k lhs, const poly_qshow_vec_k rhs);
int poly_qshow_mat_d_k_equal(const poly_qshow_mat_d_k lhs, const poly_qshow_mat_d_k rhs);
void poly_qshow_mat_d_k_dump(const poly_qshow_mat_d_k arg);

#endif /* POLY_QSHOW_MAT_D_K_H */
