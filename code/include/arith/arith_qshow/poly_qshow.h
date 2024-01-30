#ifndef POLY_QSHOW_H
#define POLY_QSHOW_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/nmod_poly.h>

#include "params.h"

#ifdef __SIZEOF_INT128__
__extension__ typedef __int128 int128;
__extension__ typedef unsigned __int128 uint128;
#endif

#define COEFFQSHOW_PACKEDBYTES 8
#define POLYQSHOW_PACKEDBYTES (PARAM_N_SHOW * COEFFQSHOW_PACKEDBYTES)

//typedef nmod_poly_t __poly_qshow;

typedef nmod_poly_t poly_qshow;

typedef int64_t coeff_qshow;

void arith_qshow_setup(void);
void arith_qshow_teardown(void);

void poly_qshow_init(poly_qshow res);
void poly_qshow_clear(poly_qshow arg);
void poly_qshow_zero(poly_qshow res);
void poly_qshow_set(poly_qshow res, const poly_qshow arg);
coeff_qshow poly_qshow_get_coeff(const poly_qshow arg, size_t n);
coeff_qshow poly_qshow_get_coeff_centered(const poly_qshow arg, size_t n);
void poly_qshow_set_coeff(poly_qshow arg, size_t n, coeff_qshow c);
void poly_qshow_neg(poly_qshow res, const poly_qshow arg);
void poly_qshow_add(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_neg(poly_qshow res, const poly_qshow arg);
void poly_qshow_sub(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_mul(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_mul_x(poly_qshow res, const poly_qshow arg);
void poly_qshow_mul_scalar(poly_qshow out, const poly_qshow lhs, const coeff_qshow rhs);
void poly_qshow_muladd_constant(poly_qshow arg, const coeff_qshow c0_lhs, const coeff_qshow c0_rhs);
void poly_qshow_shift_left(poly_qshow res, const poly_qshow arg, size_t n);
void poly_qshow_conjugate(poly_qshow out, const poly_qshow arg);
int poly_qshow_equal(const poly_qshow lhs, const poly_qshow rhs);
void poly_qshow_dump(const poly_qshow arg);
uint128 poly_qshow_sq_norm2(const poly_qshow arg);
void poly_qshow_pack(uint8_t buf[POLYQSHOW_PACKEDBYTES], const poly_qshow arg);
void coeff_qshow_pack(uint8_t buf[COEFFQSHOW_PACKEDBYTES], const coeff_qshow arg);

#endif /* POLY_QSHOW_H */
