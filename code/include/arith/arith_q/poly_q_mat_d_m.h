#ifndef POLY_Q_MAT_D_M
#define POLY_Q_MAT_D_M

#include "arith_q.h"

typedef struct {
    poly_q_vec_m rows[PARAM_D];
} __poly_q_mat_d_m;

typedef __poly_q_mat_d_m poly_q_mat_d_m[1];

void poly_q_mat_d_m_setup(void);
void poly_q_mat_d_m_teardown(void);

void poly_q_mat_d_m_init(poly_q_mat_d_m res);
void poly_q_mat_d_m_clear(poly_q_mat_d_m res);
void poly_q_mat_d_m_zero(poly_q_mat_d_m res);
void poly_q_mat_d_m_set(poly_q_mat_d_m res, const poly_q_mat_d_m arg);
void poly_q_mat_d_m_neg(poly_q_mat_d_m res, const poly_q_mat_d_m arg);
void poly_q_mat_d_m_add(poly_q_mat_d_m res, const poly_q_mat_d_m lhs, const poly_q_mat_d_m rhs);
void poly_q_mat_d_m_sub(poly_q_mat_d_m res, const poly_q_mat_d_m lhs, const poly_q_mat_d_m rhs);
void poly_q_mat_d_m_mul_vec_m(poly_q_vec_d res, const poly_q_mat_d_m lhs, const poly_q_vec_m rhs);
int poly_q_mat_d_m_equal(const poly_q_mat_d_m lhs, const poly_q_mat_d_m rhs);
void poly_q_mat_d_m_dump(const poly_q_mat_d_m arg);

#endif /* POLY_Q_MAT_D_M */
