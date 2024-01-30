#ifndef POLY_QISS_MAT_K_K_H
#define POLY_QISS_MAT_K_K_H

#include "arith_qiss.h"

typedef struct {
    poly_qiss_vec_k rows[PARAM_K_ISS];
} __poly_qiss_mat_k_k;

typedef __poly_qiss_mat_k_k poly_qiss_mat_k_k[1];

void poly_qiss_mat_k_k_setup(void);
void poly_qiss_mat_k_k_teardown(void);

void poly_qiss_mat_k_k_init(poly_qiss_mat_k_k res);
void poly_qiss_mat_k_k_clear(poly_qiss_mat_k_k res);
void poly_qiss_mat_k_k_zero(poly_qiss_mat_k_k res);
void poly_qiss_mat_k_k_set(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k arg);
void poly_qiss_mat_k_k_neg(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k arg);
void poly_qiss_mat_k_k_add(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs);
void poly_qiss_mat_k_k_sub(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs);
void poly_qiss_mat_k_k_mul_vec_k(poly_qiss_vec_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_vec_k rhs);
int poly_qiss_mat_k_k_equal(const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs);
void poly_qiss_mat_k_k_dump(const poly_qiss_mat_k_k arg);

#endif /* POLY_QISS_MAT_K_K_H */
