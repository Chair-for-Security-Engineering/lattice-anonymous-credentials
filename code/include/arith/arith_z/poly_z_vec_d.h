#ifndef POLY_Z_VEC_D_H
#define POLY_Z_VEC_D_H

#include "poly_z.h"

typedef struct {
	poly_z entries[PARAM_D];
} __poly_z_vec_d;

typedef __poly_z_vec_d poly_z_vec_d[1];

void poly_z_vec_d_setup(void);
void poly_z_vec_d_teardown(void);

void poly_z_vec_d_init(poly_z_vec_d arg);
void poly_z_vec_d_clear(poly_z_vec_d arg);
void poly_z_vec_d_zero(poly_z_vec_d arg);
void poly_z_vec_d_set(poly_z_vec_d res, const poly_z_vec_d arg);
void poly_z_vec_d_get_poly(poly_z res, const poly_z_vec_d arg, size_t pos);
void poly_z_vec_d_set_poly(poly_z_vec_d res, const poly_z arg, size_t pos);
void poly_z_vec_d_neg(poly_z_vec_d res, const poly_z_vec_d arg);
void poly_z_vec_d_add(poly_z_vec_d res, const poly_z_vec_d lhs, const poly_z_vec_d rhs);
void poly_z_vec_d_sub(poly_z_vec_d res, const poly_z_vec_d lhs, const poly_z_vec_d rhs);
void poly_z_vec_d_mul_scalar(poly_z_vec_d res, const poly_z_vec_d arg, const poly_z fac);
void poly_z_vec_d_mul_inner(poly_z res, const poly_z_vec_d lhs, const poly_z_vec_d rhs);
int poly_z_vec_d_equal(const poly_z_vec_d lhs, const poly_z_vec_d rhs);
void poly_z_vec_d_dump(const poly_z_vec_d arg);

#endif /* POLY_Q_VEC_D_H */
