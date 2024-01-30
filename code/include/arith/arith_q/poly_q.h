#ifndef POLY_Q_H
#define POLY_Q_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/nmod_poly.h>

#include "params.h"

typedef nmod_poly_t poly_q;

typedef int64_t coeff_q;

void arith_q_setup(void);
void arith_q_teardown(void);
void poly_q_init(poly_q res);
void poly_q_clear(poly_q arg);
void poly_q_zero(poly_q res);
void poly_q_set(poly_q res, const poly_q arg);
coeff_q poly_q_get_coeff(const poly_q arg, size_t n);
coeff_q poly_q_get_coeff_centered(const poly_q arg, size_t n);
void poly_q_set_coeff(poly_q arg, size_t n, coeff_q c);
void poly_q_from_bits(poly_q arg, const uint8_t coeffs[PARAM_N / 8]);
void poly_q_neg(poly_q res, const poly_q arg);
void poly_q_add(poly_q res, const poly_q lhs, const poly_q rhs);
void poly_q_sub(poly_q res, const poly_q lhs, const poly_q rhs);
void poly_q_mul(poly_q res, const poly_q lhs, const poly_q rhs);
void poly_q_mul_scalar(poly_q res, const poly_q arg, const coeff_q fac);
void poly_q_shift_left(poly_q res, const poly_q arg, size_t n);
void poly_q_shift_right(poly_q res, const poly_q arg, size_t n);
void poly_q_conjugate(poly_q res, const poly_q arg);
int poly_q_equal(const poly_q lhs, const poly_q rhs);
void poly_q_dump(const poly_q arg);
uint64_t poly_q_sq_norm2(const poly_q arg);
int64_t poly_q_weight(const poly_q arg);
void poly_q_invert(poly_q res, const poly_q arg);

#endif /* POLY_Q_H */
