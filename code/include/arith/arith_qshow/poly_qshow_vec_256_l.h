#ifndef POLY_QSHOW_VEC_256_L_H
#define POLY_QSHOW_VEC_256_L_H

#include "poly_qshow.h"

#define POLYQSHOW_VEC256L_PACKEDBYTES (PARAM_ARP_DIV_N_L_SHOW * POLYQSHOW_PACKEDBYTES)

typedef struct {
	poly_qshow entries[PARAM_ARP_DIV_N_L_SHOW];
} __poly_qshow_vec_256_l;

typedef __poly_qshow_vec_256_l poly_qshow_vec_256_l[1];

void poly_qshow_vec_256_l_setup(void);
void poly_qshow_vec_256_l_teardown(void);

void poly_qshow_vec_256_l_init(poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_clear(poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_zero(poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_set(poly_qshow_vec_256_l res, const poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_get_poly(poly_qshow res, const poly_qshow_vec_256_l arg, size_t pos);
void poly_qshow_vec_256_l_set_poly(poly_qshow_vec_256_l res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_256_l_neg(poly_qshow_vec_256_l res, const poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_add(poly_qshow_vec_256_l res, const poly_qshow_vec_256_l lhs, const poly_qshow_vec_256_l rhs);
void poly_qshow_vec_256_l_sub(poly_qshow_vec_256_l res, const poly_qshow_vec_256_l lhs, const poly_qshow_vec_256_l rhs);
void poly_qshow_vec_256_l_mul_scalar(poly_qshow_vec_256_l res, const poly_qshow_vec_256_l arg, const coeff_qshow fac);
void poly_qshow_vec_256_l_mul_poly_qshow(poly_qshow_vec_256_l out, const poly_qshow_vec_256_l lhs, const poly_qshow rhs);
void poly_qshow_vec_256_l_mul_inner(poly_qshow res, const poly_qshow_vec_256_l lhs, const poly_qshow_vec_256_l rhs);
uint64_t poly_qshow_vec_256_l_norm2(const poly_qshow_vec_256_l arg);
int poly_qshow_vec_256_l_equal(const poly_qshow_vec_256_l lhs, const poly_qshow_vec_256_l rhs);
void poly_qshow_vec_256_l_dump(const poly_qshow_vec_256_l arg);
void poly_qshow_vec_256_l_pack(uint8_t buf[POLYQSHOW_VEC256L_PACKEDBYTES], const poly_qshow_vec_256_l arg);

#endif /* POLY_QSHOW_VEC_256_L_H */
