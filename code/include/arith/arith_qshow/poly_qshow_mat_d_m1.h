#ifndef POL_QSHOW_MAT_D_M1_H
#define POL_QSHOW_MAT_D_M1_H

#include "arith_qshow.h"

typedef struct {
    poly_qshow_vec_m1 rows[PARAM_D_SHOW];
} __poly_qshow_mat_d_m1;

typedef __poly_qshow_mat_d_m1 poly_qshow_mat_d_m1[1];

void poly_qshow_mat_d_m1_setup(void);
void poly_qshow_mat_d_m1_teardown(void);

void poly_qshow_mat_d_m1_init(poly_qshow_mat_d_m1 res);
void poly_qshow_mat_d_m1_clear(poly_qshow_mat_d_m1 res);
void poly_qshow_mat_d_m1_zero(poly_qshow_mat_d_m1 res);
void poly_qshow_mat_d_m1_set(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 arg);
void poly_qshow_mat_d_m1_neg(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 arg);
void poly_qshow_mat_d_m1_add(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs);
void poly_qshow_mat_d_m1_sub(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs);
void poly_qshow_mat_d_m1_mul_vec_m1(poly_qshow_vec_d res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_vec_m1 rhs);
int poly_qshow_mat_d_m1_equal(const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs);
void poly_qshow_mat_d_m1_dump(const poly_qshow_mat_d_m1 arg);

#endif /* POL_QSHOW_MAT_D_M1_H */
