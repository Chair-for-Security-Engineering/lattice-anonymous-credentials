#ifndef POL_QSHOW_MAT_D_M2_H
#define POL_QSHOW_MAT_D_M2_H

#include "arith_qshow.h"

typedef struct {
    poly_qshow_vec_m2 rows[PARAM_D_SHOW];
} __poly_qshow_mat_d_m2;

typedef __poly_qshow_mat_d_m2 poly_qshow_mat_d_m2[1];

void poly_qshow_mat_d_m2_setup(void);
void poly_qshow_mat_d_m2_teardown(void);

void poly_qshow_mat_d_m2_init(poly_qshow_mat_d_m2 res);
void poly_qshow_mat_d_m2_clear(poly_qshow_mat_d_m2 res);
void poly_qshow_mat_d_m2_zero(poly_qshow_mat_d_m2 res);
void poly_qshow_mat_d_m2_set(poly_qshow_mat_d_m2 res, const poly_qshow_mat_d_m2 arg);
void poly_qshow_mat_d_m2_neg(poly_qshow_mat_d_m2 res, const poly_qshow_mat_d_m2 arg);
void poly_qshow_mat_d_m2_add(poly_qshow_mat_d_m2 res, const poly_qshow_mat_d_m2 lhs, const poly_qshow_mat_d_m2 rhs);
void poly_qshow_mat_d_m2_sub(poly_qshow_mat_d_m2 res, const poly_qshow_mat_d_m2 lhs, const poly_qshow_mat_d_m2 rhs);
void poly_qshow_mat_d_m2_mul_vec_m2(poly_qshow_vec_d res, const poly_qshow_mat_d_m2 lhs, const poly_qshow_vec_m2 rhs);
int poly_qshow_mat_d_m2_equal(const poly_qshow_mat_d_m2 lhs, const poly_qshow_mat_d_m2 rhs);
void poly_qshow_mat_d_m2_dump(const poly_qshow_mat_d_m2 arg);

#endif /* POL_QSHOW_MAT_D_M2_H */
