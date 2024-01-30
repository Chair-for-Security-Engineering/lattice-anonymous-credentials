#ifndef POLY_QSHOW_VEC_D_H
#define POLY_QSHOW_VEC_D_H

#include "poly_qshow.h"

#define POLYQSHOW_VECD_PACKEDBYTES (PARAM_D_SHOW * POLYQSHOW_PACKEDBYTES)

typedef struct {
	poly_qshow entries[PARAM_D_SHOW];
} __poly_qshow_vec_d;

typedef __poly_qshow_vec_d poly_qshow_vec_d[1];

void poly_qshow_vec_d_setup(void);
void poly_qshow_vec_d_teardown(void);

void poly_qshow_vec_d_init(poly_qshow_vec_d arg);
void poly_qshow_vec_d_clear(poly_qshow_vec_d arg);
void poly_qshow_vec_d_zero(poly_qshow_vec_d arg);
void poly_qshow_vec_d_set(poly_qshow_vec_d res, const poly_qshow_vec_d arg);
void poly_qshow_vec_d_get_poly(poly_qshow res, const poly_qshow_vec_d arg, size_t pos);
void poly_qshow_vec_d_set_poly(poly_qshow_vec_d res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_d_neg(poly_qshow_vec_d res, const poly_qshow_vec_d arg);
void poly_qshow_vec_d_add(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs);
void poly_qshow_vec_d_sub(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs);
void poly_qshow_vec_d_mul_scalar(poly_qshow_vec_d res, const poly_qshow_vec_d arg, const coeff_qshow fac);
void poly_qshow_vec_d_mul_poly_qshow(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow rhs);
void poly_qshow_vec_d_mul_inner(poly_qshow res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs);
uint64_t poly_qshow_vec_d_norm2(const poly_qshow_vec_d arg);
int poly_qshow_vec_d_equal(const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs);
void poly_qshow_vec_d_dump(const poly_qshow_vec_d arg);
void poly_qshow_vec_d_pack(uint8_t buf[POLYQSHOW_VECD_PACKEDBYTES], const poly_qshow_vec_d arg);

#endif /* POLY_QSHOW_VEC_D_H */
