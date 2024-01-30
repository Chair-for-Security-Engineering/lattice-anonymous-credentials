#ifndef POLY_QISS_VEC_256_L_H
#define POLY_QISS_VEC_256_L_H

#include "poly_qiss.h"

#define POLYQISS_VEC256L_PACKEDBYTES (PARAM_ARP_DIV_N_L_ISS * POLYQISS_PACKEDBYTES)

typedef struct {
	poly_qiss entries[PARAM_ARP_DIV_N_L_ISS];
} __poly_qiss_vec_256_l;

typedef __poly_qiss_vec_256_l poly_qiss_vec_256_l[1];

void poly_qiss_vec_256_l_setup(void);
void poly_qiss_vec_256_l_teardown(void);

void poly_qiss_vec_256_l_init(poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_clear(poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_zero(poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_set(poly_qiss_vec_256_l res, const poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_get_poly(poly_qiss res, const poly_qiss_vec_256_l arg, size_t pos);
void poly_qiss_vec_256_l_set_poly(poly_qiss_vec_256_l res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_256_l_neg(poly_qiss_vec_256_l res, const poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_add(poly_qiss_vec_256_l res, const poly_qiss_vec_256_l lhs, const poly_qiss_vec_256_l rhs);
void poly_qiss_vec_256_l_sub(poly_qiss_vec_256_l res, const poly_qiss_vec_256_l lhs, const poly_qiss_vec_256_l rhs);
void poly_qiss_vec_256_l_mul_scalar(poly_qiss_vec_256_l res, const poly_qiss_vec_256_l arg, const coeff_qiss fac);
void poly_qiss_vec_256_l_mul_poly_qiss(poly_qiss_vec_256_l out, const poly_qiss_vec_256_l lhs, const poly_qiss rhs);
void poly_qiss_vec_256_l_mul_inner(poly_qiss res, const poly_qiss_vec_256_l lhs, const poly_qiss_vec_256_l rhs);
uint64_t poly_qiss_vec_256_l_norm2(const poly_qiss_vec_256_l arg);
int poly_qiss_vec_256_l_equal(const poly_qiss_vec_256_l lhs, const poly_qiss_vec_256_l rhs);
void poly_qiss_vec_256_l_dump(const poly_qiss_vec_256_l arg);
void poly_qiss_vec_256_l_pack(uint8_t buf[POLYQISS_VEC256L_PACKEDBYTES], const poly_qiss_vec_256_l arg);

#endif /* POLY_QISS_VEC_256_L_H */
