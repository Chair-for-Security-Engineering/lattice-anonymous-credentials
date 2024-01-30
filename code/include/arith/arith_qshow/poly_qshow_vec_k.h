#ifndef POLY_QSHOW_VEC_K_H
#define POLY_QSHOW_VEC_K_H

#include "poly_qshow.h"

#define POLYQSHOW_VECK_PACKEDBYTES (POLYQSHOW_PACKEDBYTES * PARAM_K_SHOW)

typedef struct {
	poly_qshow entries[PARAM_K_SHOW];
} __poly_qshow_vec_k;

typedef __poly_qshow_vec_k poly_qshow_vec_k[1];

void poly_qshow_vec_k_setup(void);
void poly_qshow_vec_k_teardown(void);

void poly_qshow_vec_k_init(poly_qshow_vec_k arg);
void poly_qshow_vec_k_clear(poly_qshow_vec_k arg);
void poly_qshow_vec_k_zero(poly_qshow_vec_k arg);
void poly_qshow_vec_k_set(poly_qshow_vec_k res, const poly_qshow_vec_k arg);
void poly_qshow_vec_k_get_poly(poly_qshow res, const poly_qshow_vec_k arg, size_t pos);
void poly_qshow_vec_k_set_poly(poly_qshow_vec_k res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_k_neg(poly_qshow_vec_k res, const poly_qshow_vec_k arg);
void poly_qshow_vec_k_add(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs);
void poly_qshow_vec_k_sub(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs);
void poly_qshow_vec_k_mul_scalar(poly_qshow_vec_k res, const poly_qshow_vec_k arg, const coeff_qshow fac);
void poly_qshow_vec_k_mul_poly_qshow(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow rhs);
void poly_qshow_vec_k_mul_inner(poly_qshow res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs);
void poly_qshow_vec_k_conjugate(poly_qshow_vec_k res, const poly_qshow_vec_k arg);
uint64_t poly_qshow_vec_k_norm2(const poly_qshow_vec_k arg);
int poly_qshow_vec_k_equal(const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs);
void poly_qshow_vec_k_dump(const poly_qshow_vec_k arg);
void poly_qshow_vec_k_pack(uint8_t buf[POLYQSHOW_VECK_PACKEDBYTES], const poly_qshow_vec_k arg);

#endif /* POLY_QSHOW_VEC_K_H */
