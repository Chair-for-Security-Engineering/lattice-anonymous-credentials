#include "arith_z.h"
#include "macros.h"

static fmpz_t TMP;
static fmpz_t AUX;

/*************************************************
* Name:        arith_z_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              multiprecision integers. This is strictly required 
*              and must be called once before any other function 
*              from here is used.
**************************************************/
void arith_z_setup(void) {
  fmpz_init(TMP);
  fmpz_init(AUX);
}

/*************************************************
* Name:        arith_z_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              multiprecision integers. This is strictly required 
*              and must be called once at the very end to release 
*              any resources.
**************************************************/
void arith_z_teardown(void) {
  fmpz_clear(TMP);
  fmpz_clear(AUX);
}

/*************************************************
* Name:        poly_z_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_z res: polynomial to be initialized
**************************************************/
void poly_z_init(poly_z res) {
  fmpz_poly_init(res);
}

/*************************************************
* Name:        poly_z_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_z arg: polynomial to be cleared
**************************************************/
void poly_z_clear(poly_z arg) {
  fmpz_poly_clear(arg);
}

/*************************************************
* Name:        poly_z_zero
*
* Description: Set an initialized polynomial to zero
* 
* Arguments:   - poly_z res: polynomial to be zeroized (initialized)
**************************************************/
void poly_z_zero(poly_z res) {
  fmpz_poly_zero(res);
}

/*************************************************
* Name:        poly_z_set
*
* Description: Set a polynomial equal to another polynomial
* 
* Arguments:   - poly_z res: polynomial to be set (initialized)
*              - const poly_z arg: polynomial to be read
**************************************************/
void poly_z_set(poly_z res, const poly_z arg) {
  fmpz_poly_set(res, arg);
}

/*************************************************
* Name:        poly_z_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - coeff_z c: the fmpz to host the coefficient
*              - const poly_z res: polynomial to be read
*              - size_t n: degree of the coefficient to be read
**************************************************/
void poly_z_get_coeff(coeff_z c, const poly_z arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot get coefficient of poly at given position.");
  fmpz_poly_get_coeff_fmpz(c, arg, (signed long) n);
}

/*************************************************
* Name:        poly_z_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - poly_z arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_z c: the new fmpz coefficient
**************************************************/
void poly_z_set_coeff(poly_z arg, size_t n, coeff_z c) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot set coefficient of poly at given position.");
  fmpz_poly_set_coeff_fmpz(arg, (signed long) n, c);
}

/*************************************************
* Name:        poly_z_set_coeff_si
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - poly_z arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - int64_t c: the new signed 64-bit int coefficient
**************************************************/
void poly_z_set_coeff_si(poly_z arg, size_t n, int64_t c) {
  fmpz_set_si(TMP, c);
  poly_z_set_coeff(arg, n, TMP);
}

/*************************************************
* Name:        poly_z_neg
*
* Description: Negate a polynomial coefficient-wise
* 
* Arguments:   - poly_z res: polynomial to host the negation (initialized)
*              - const poly_z arg: polynomial to be negated
**************************************************/
void poly_z_neg(poly_z res, const poly_z arg) {
  fmpz_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_z_add
*
* Description: Add two polynomials 
* 
* Arguments:   - poly_z res: polynomial to host the sum (initialized)
*              - const poly_z lhs: first polynomial summand
*              - const poly_z rhs: second polynomial summand
**************************************************/
void poly_z_add(poly_z res, const poly_z lhs, const poly_z rhs) {
  fmpz_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_z_sub
*
* Description: Substract two polynomials 
* 
* Arguments:   - poly_z res: polynomial to host the difference (initialized)
*              - const poly_z lhs: first polynomial term
*              - const poly_z rhs: second polynomial term
**************************************************/
void poly_z_sub(poly_z res, const poly_z lhs, const poly_z rhs) {
  fmpz_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_z_shift_left
*
* Description: Shift the coefficients of a polynomial to
*              the left by n places. Corresponds to a 
*              multiplication by x^n (NOT REDUCED MOD X^PARAM_N + 1)
*              Trailing zeros are inserted
* 
* Arguments:   - poly_z res: polynomial to host the shift (initialized)
*              - const poly_z arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_z_shift_left(poly_z res, const poly_z arg, size_t n) {
  fmpz_poly_shift_left(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_z_shift_right
*
* Description: Shift the coefficients of a polynomial to
*              the right by n places. Corresponds to a 
*              division by x^n (NOT REDUCED MOD X^PARAM_N + 1)
*              Leading zeros are inserted
* 
* Arguments:   - poly_z res: polynomial to host the shift (initialized)
*              - const poly_z arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_z_shift_right(poly_z res, const poly_z arg, size_t n) {
  fmpz_poly_shift_right(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_z_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N + 1
* 
* Arguments:   - poly_z res: polynomial to host the multiplication (initialized)
*              - const poly_z lhs: first polynomial factor
*              - const poly_z rhs: second polynomial factor
**************************************************/
void poly_z_mul(poly_z res, const poly_z lhs, const poly_z rhs) {
  // multiplication gives a polynomial of degree < 2*PARAM_N
  fmpz_poly_mul(res, lhs, rhs);
  // modulo X^PARAM_N + 1 simply substracts the coeffs of degree N,...,2N-1
  // to those of degree 0,...,N-1
  for (size_t i = 0; i < PARAM_N; i++) {
    fmpz_poly_get_coeff_fmpz(TMP,  res, (signed long)i);
    fmpz_poly_get_coeff_fmpz(AUX, res, (signed long)i + PARAM_N);
    fmpz_sub(TMP, TMP, AUX);
    fmpz_poly_set_coeff_fmpz(res, (signed long)i, TMP);
    fmpz_poly_set_coeff_ui(res, (signed long)i + PARAM_N, 0);
  }
}

/*************************************************
* Name:        poly_z_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
* 
* Arguments:   - poly_z res: polynomial to host the multiplication (initialized)
*              - const poly_z arg: first polynomial factor
*              - coeff_z fac: second scalar factor
**************************************************/
void poly_z_mul_scalar(poly_z res, const poly_z arg, coeff_z fac) {
  fmpz_poly_scalar_mul_fmpz(res, arg, fac);
}

/*************************************************
* Name:        poly_z_equal
*
* Description: Equality test between two polynomials
* 
* Arguments:   - const poly_z lhs: first polynomial
*              - const poly_z rhs: second polynomial
* 
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_z_equal(const poly_z lhs, const poly_z rhs) {
  return fmpz_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_z_dump
*
* Description: Print a polynomial
* 
* Arguments:   - const poly_z arg: the polynomial to be printed
**************************************************/
void poly_z_dump(const poly_z arg) {
  char* str = fmpz_poly_get_str_pretty(arg, "x");
  printf("%s", str);
  flint_free(str);
}
