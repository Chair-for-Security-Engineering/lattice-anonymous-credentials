#ifndef POLY_Q_MAT_D_D_H
#define POLY_Q_MAT_D_D_H

#include "poly_q.h"
#include "poly_q_vec_d.h"

typedef struct {
    poly_q_vec_d rows[PARAM_D];
} __poly_q_mat_d_d;

typedef __poly_q_mat_d_d poly_q_mat_d_d[1];

void poly_q_mat_d_d_setup(void);
void poly_q_mat_d_d_teardown(void);

void poly_q_mat_d_d_init(poly_q_mat_d_d res);
void poly_q_mat_d_d_clear(poly_q_mat_d_d res);
void poly_q_mat_d_d_zero(poly_q_mat_d_d res);
void poly_q_mat_d_d_set(poly_q_mat_d_d res, const poly_q_mat_d_d arg);
void poly_q_mat_d_d_neg(poly_q_mat_d_d res, const poly_q_mat_d_d arg);
void poly_q_mat_d_d_add(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs);
void poly_q_mat_d_d_sub(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs);
void poly_q_mat_d_d_mul_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs);
void poly_q_mat_d_d_muladd_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs);
void poly_q_mat_d_d_mulsub_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs);
void poly_q_mat_d_d_mul_mat_d_d(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs);
void poly_q_mat_d_d_conjugate(poly_q_mat_d_d res, const poly_q_mat_d_d arg);
int polq_q_mat_d_d_equal(const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs);
void poly_q_mat_d_d_dump(const poly_q_mat_d_d arg);

#endif /* POLY_Q_MAT_D_D_H */
