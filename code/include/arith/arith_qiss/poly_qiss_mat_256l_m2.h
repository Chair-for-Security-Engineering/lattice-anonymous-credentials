#ifndef POL_QISS_MAT_256L_M2_H
#define POL_QISS_MAT_256L_M2_H

#include "arith_qiss.h"

typedef struct {
    poly_qiss_vec_m2 rows[PARAM_ARP_DIV_N_L_ISS];
} __poly_qiss_mat_256l_m2;

typedef __poly_qiss_mat_256l_m2 poly_qiss_mat_256l_m2[1];

void poly_qiss_mat_256l_m2_setup(void);
void poly_qiss_mat_256l_m2_teardown(void);

void poly_qiss_mat_256l_m2_init(poly_qiss_mat_256l_m2 res);
void poly_qiss_mat_256l_m2_clear(poly_qiss_mat_256l_m2 res);
void poly_qiss_mat_256l_m2_zero(poly_qiss_mat_256l_m2 res);
void poly_qiss_mat_256l_m2_set(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 arg);
void poly_qiss_mat_256l_m2_neg(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 arg);
void poly_qiss_mat_256l_m2_add(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs);
void poly_qiss_mat_256l_m2_sub(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs);
void poly_qiss_mat_256l_m2_mul_vec_m2(poly_qiss_vec_256_l res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_vec_m2 rhs);
int poly_qiss_mat_256l_m2_equal(const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs);
void poly_qiss_mat_256l_m2_dump(const poly_qiss_mat_256l_m2 arg);

#endif /* POL_QISS_MAT_256L_M2_H */
