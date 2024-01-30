#ifndef POLY_REAL_MAT_2D_2D_H
#define POLY_REAL_MAT_2D_2D_H

#include "poly_real.h"
#include "poly_real_vec_2d.h"

typedef struct {
    poly_real_vec_2d rows[2*PARAM_D];
} __poly_real_mat_2d_2d;

typedef __poly_real_mat_2d_2d poly_real_mat_2d_2d[1];

void poly_real_mat_2d_2d_init(poly_real_mat_2d_2d res);
void poly_real_mat_2d_2d_clear(poly_real_mat_2d_2d res);

#endif /* POLY_REAL_MAT_2D_2D_H */
