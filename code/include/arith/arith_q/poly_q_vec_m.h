#ifndef POLY_Q_VEC_M_H
#define POLY_Q_VEC_M_H

#include "poly_q.h"

typedef struct {
	poly_q entries[PARAM_M];
} __poly_q_vec_m;

typedef __poly_q_vec_m poly_q_vec_m[1];

void poly_q_vec_m_setup(void);
void poly_q_vec_m_teardown(void);

void poly_q_vec_m_init(poly_q_vec_m arg);
void poly_q_vec_m_clear(poly_q_vec_m arg);
void poly_q_vec_m_zero(poly_q_vec_m arg);
void poly_q_vec_m_set(poly_q_vec_m res, const poly_q_vec_m arg);
void poly_q_vec_m_get_poly(poly_q res, const poly_q_vec_m arg, size_t pos);
void poly_q_vec_m_set_poly(poly_q_vec_m res, const poly_q arg, size_t pos);
void poly_q_vec_m_neg(poly_q_vec_m res, const poly_q_vec_m arg);
void poly_q_vec_m_add(poly_q_vec_m res, const poly_q_vec_m lhs, const poly_q_vec_m rhs);
void poly_q_vec_m_sub(poly_q_vec_m res, const poly_q_vec_m lhs, const poly_q_vec_m rhs);
void poly_q_vec_m_mul_scalar(poly_q_vec_m res, const poly_q_vec_m arg, const poly_q fac);
void poly_q_vec_m_mul_inner(poly_q res, const poly_q_vec_m lhs, const poly_q_vec_m rhs);
int poly_q_vec_m_equal(const poly_q_vec_m lhs, const poly_q_vec_m rhs);
void poly_q_vec_m_dump(const poly_q_vec_m arg);

#endif /* POLY_Q_VEC_D_H */
