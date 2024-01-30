#ifndef POLY_QSHOW_VEC_L_H
#define POLY_QSHOW_VEC_L_H

#include "poly_qshow.h"

#define POLYQSHOW_VECL_PACKEDBYTES (PARAM_L_SHOW * POLYQSHOW_PACKEDBYTES)

typedef struct {
	poly_qshow entries[PARAM_L_SHOW];
} __poly_qshow_vec_l;

typedef __poly_qshow_vec_l poly_qshow_vec_l[1];

void poly_qshow_vec_l_setup(void);
void poly_qshow_vec_l_teardown(void);

void poly_qshow_vec_l_init(poly_qshow_vec_l arg);
void poly_qshow_vec_l_clear(poly_qshow_vec_l arg);
void poly_qshow_vec_l_zero(poly_qshow_vec_l arg);
void poly_qshow_vec_l_set(poly_qshow_vec_l res, const poly_qshow_vec_l arg);
void poly_qshow_vec_l_get_poly(poly_qshow res, const poly_qshow_vec_l arg, size_t pos);
void poly_qshow_vec_l_set_poly(poly_qshow_vec_l res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_l_neg(poly_qshow_vec_l res, const poly_qshow_vec_l arg);
void poly_qshow_vec_l_add(poly_qshow_vec_l res, const poly_qshow_vec_l lhs, const poly_qshow_vec_l rhs);
void poly_qshow_vec_l_sub(poly_qshow_vec_l res, const poly_qshow_vec_l lhs, const poly_qshow_vec_l rhs);
void poly_qshow_vec_l_mul_scalar(poly_qshow_vec_l res, const poly_qshow_vec_l arg, const coeff_qshow fac);
void poly_qshow_vec_l_mul_inner(poly_qshow res, const poly_qshow_vec_l lhs, const poly_qshow_vec_l rhs);
uint64_t poly_qshow_vec_l_norm2(const poly_qshow_vec_l arg);
int poly_qshow_vec_l_equal(const poly_qshow_vec_l lhs, const poly_qshow_vec_l rhs);
void poly_qshow_vec_l_dump(const poly_qshow_vec_l arg);
void poly_qshow_vec_l_pack(uint8_t buf[POLYQSHOW_VECL_PACKEDBYTES], const poly_qshow_vec_l arg);

#endif /* POLY_QSHOW_VEC_L_H */
