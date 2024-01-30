#ifndef POLY_QISS_VEC_L_H
#define POLY_QISS_VEC_L_H

#include "poly_qiss.h"

#define POLYQISS_VECL_PACKEDBYTES (PARAM_L_ISS * POLYQISS_PACKEDBYTES)

typedef struct {
	poly_qiss entries[PARAM_L_ISS];
} __poly_qiss_vec_l;

typedef __poly_qiss_vec_l poly_qiss_vec_l[1];

void poly_qiss_vec_l_setup(void);
void poly_qiss_vec_l_teardown(void);

void poly_qiss_vec_l_init(poly_qiss_vec_l arg);
void poly_qiss_vec_l_clear(poly_qiss_vec_l arg);
void poly_qiss_vec_l_zero(poly_qiss_vec_l arg);
void poly_qiss_vec_l_set(poly_qiss_vec_l res, const poly_qiss_vec_l arg);
void poly_qiss_vec_l_get_poly(poly_qiss res, const poly_qiss_vec_l arg, size_t pos);
void poly_qiss_vec_l_set_poly(poly_qiss_vec_l res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_l_neg(poly_qiss_vec_l res, const poly_qiss_vec_l arg);
void poly_qiss_vec_l_add(poly_qiss_vec_l res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs);
void poly_qiss_vec_l_sub(poly_qiss_vec_l res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs);
void poly_qiss_vec_l_mul_scalar(poly_qiss_vec_l res, const poly_qiss_vec_l arg, const coeff_qiss fac);
void poly_qiss_vec_l_mul_inner(poly_qiss res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs);
uint64_t poly_qiss_vec_l_norm2(const poly_qiss_vec_l arg);
int poly_qiss_vec_l_equal(const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs);
void poly_qiss_vec_l_dump(const poly_qiss_vec_l arg);
void poly_qiss_vec_l_pack(uint8_t buf[POLYQISS_VECL_PACKEDBYTES], const poly_qiss_vec_l arg);

#endif /* POLY_QISS_VEC_L_H */
