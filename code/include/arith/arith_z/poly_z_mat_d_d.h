#ifndef POLY_Z_MAT_D_D_H
#define POLY_Z_MAT_D_D_H

#include "poly_z.h"
#include "poly_z_vec_d.h"

typedef struct {
    poly_z_vec_d rows[PARAM_D];
} __poly_z_mat_d_d;

typedef __poly_z_mat_d_d poly_z_mat_d_d[1];

void poly_z_mat_d_d_setup(void);
void poly_z_mat_d_d_teardown(void);

void poly_z_mat_d_d_init(poly_z_mat_d_d res);
void poly_z_mat_d_d_clear(poly_z_mat_d_d res);
void poly_z_mat_d_d_zero(poly_z_mat_d_d res);
void poly_z_mat_d_d_set(poly_z_mat_d_d res, const poly_z_mat_d_d arg);
void poly_z_mat_d_d_neg(poly_z_mat_d_d res, const poly_z_mat_d_d arg);
void poly_z_mat_d_d_add(poly_z_mat_d_d res, const poly_z_mat_d_d lhs, const poly_z_mat_d_d rhs);
void poly_z_mat_d_d_sub(poly_z_mat_d_d res, const poly_z_mat_d_d lhs, const poly_z_mat_d_d rhs);
void poly_z_mat_d_d_mul_vec_d(poly_z_vec_d res, const poly_z_mat_d_d lhs, const poly_z_vec_d rhs);
void poly_z_mat_d_d_mul_mat_d_d(poly_z_mat_d_d res, const poly_z_mat_d_d lhs, const poly_z_mat_d_d rhs);
int poly_z_mat_d_d_equal(const poly_z_mat_d_d lhs, const poly_z_mat_d_d rhs);
void poly_z_mat_d_d_dump(const poly_z_mat_d_d arg);

#endif /* POLY_Q_MAT_D_D_H */
