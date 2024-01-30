#ifndef POLY_Z_H
#define POLY_Z_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "params.h"

typedef fmpz_poly_t poly_z;

typedef fmpz_t coeff_z;

void arith_z_setup(void);
void arith_z_teardown(void);

void poly_z_init(poly_z res);
void poly_z_clear(poly_z arg);
void poly_z_zero(poly_z res);
void poly_z_set(poly_z res, const poly_z arg);
void poly_z_get_coeff(coeff_z c, const poly_z arg, size_t n);
void poly_z_set_coeff(poly_z arg, size_t n, coeff_z c);
void poly_z_set_coeff_si(poly_z arg, size_t n, int64_t c);
void poly_z_neg(poly_z res, const poly_z arg);
void poly_z_add(poly_z res, const poly_z lhs, const poly_z rhs);
void poly_z_sub(poly_z res, const poly_z lhs, const poly_z rhs);
void poly_z_shift_left(poly_z res, const poly_z arg, size_t n);
void poly_z_shift_right(poly_z res, const poly_z arg, size_t n);
void poly_z_mul(poly_z res, const poly_z lhs, const poly_z rhs);
void poly_z_mul_scalar(poly_z res, const poly_z arg, coeff_z fac);
int poly_z_equal(const poly_z lhs, const poly_z rhs);
void poly_z_dump(const poly_z arg);

#endif /* POLY_Q_H */
