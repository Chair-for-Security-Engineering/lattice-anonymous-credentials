#include "arith_q.h"
#include "macros.h"

static nmod_poly_t POLY_F;
static nmod_poly_t POLY_F_REV_INV;

static poly_q TMP;
static poly_q AUX;

/*************************************************
* Name:        nmod_poly_invert
*
* Description: Invert an nmod_poly to later compute faster multiplication
*              using the precomputed inverse of reverse(x^PARAM_N + 1)
*
* Arguments:   - nmod_poly_t res: polynomial to host the inverse
*              - const nmod_poly_t arg: polynomial to invert
*              - const nmod_poly_t mod: polynomial to reduce with
*              - int64_t q: modulus
**************************************************/
static void nmod_poly_invert(nmod_poly_t res, const nmod_poly_t arg, const nmod_poly_t mod, int64_t q) {
  nmod_poly_t G, S;
  nmod_poly_init(G, q);
  nmod_poly_init(S, q);
  nmod_poly_xgcd(G, S, res, mod, arg);
  ASSERT_DEBUG(nmod_poly_is_one(G), "GDC != 1");
  nmod_poly_mulmod(G, res, arg, mod);
  ASSERT_DEBUG(nmod_poly_is_one(G), "arg * res % mod != 1");
  nmod_poly_clear(G);
  nmod_poly_clear(S);
}

/*************************************************
* Name:        arith_q_setup
*
* Description: Initialize and setup the backend for arithmetic
*              modulo PARAM_Q. This is strictly required
*              and must be called once before any other function
*              from here is used.
**************************************************/
void arith_q_setup(void) {
  ASSERT_DEBUG(PARAM_N % 8 == 0, "Illegal parameter: ring degree must be divisible by 8 to allow `poly_from_bits` function.");
  nmod_poly_init(POLY_F, PARAM_Q);
  nmod_poly_set_coeff_ui(POLY_F, 0, 1);
  nmod_poly_set_coeff_ui(POLY_F, PARAM_N, 1);

  nmod_poly_t f_len;
  nmod_poly_init(f_len, PARAM_Q);
  nmod_poly_set_coeff_ui(f_len, PARAM_N + 1, 1);

  nmod_poly_t f_rev;
  nmod_poly_init(f_rev, PARAM_Q);

  nmod_poly_init(POLY_F_REV_INV, PARAM_Q);
  nmod_poly_reverse(f_rev, POLY_F, PARAM_N + 1);
  nmod_poly_invert(POLY_F_REV_INV, f_rev, f_len, PARAM_Q);

  nmod_poly_clear(f_len);
  nmod_poly_clear(f_rev);

  poly_q_init(TMP);
  poly_q_init(AUX);
}

/*************************************************
* Name:        arith_q_teardown
*
* Description: Clean up and teardown the backend for arithmetic
*              modulo PARAM_Q. This is strictly required
*              and must be called once at the very end to release
*              any resources.
**************************************************/
void arith_q_teardown(void) {
  poly_q_clear(TMP);
  poly_q_clear(AUX);

  nmod_poly_clear(POLY_F);
  nmod_poly_clear(POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_q_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations
*              are done with/on the polynomial.
*
* Arguments:   - poly_q res: polynomial to be initialized
**************************************************/
void poly_q_init(poly_q res) {
	nmod_poly_init(res, PARAM_Q);
}

/*************************************************
* Name:        poly_q_clear
*
* Description: Clears a polynomial and releases all associated memory.
*              This is strictly required to avoid memory leaks and the
*              polynomial must not be used again (unless reinitialized).
*
* Arguments:   - poly_q arg: polynomial to be cleared
**************************************************/
void poly_q_clear(poly_q arg) {
	nmod_poly_clear(arg);
}

/*************************************************
* Name:        poly_q_zero
*
* Description: Set an initialized polynomial to zero
*
* Arguments:   - poly_q res: polynomial to be zeroized (initialized)
**************************************************/
void poly_q_zero(poly_q res) {
	nmod_poly_zero(res);
}

/*************************************************
* Name:        poly_q_set
*
* Description: Set a polynomial equal to another polynomial.
*              Coefficients are reduced mod PARAM_Q
*
* Arguments:   - poly_q res: polynomial to be set (initialized)
*              - const poly_q arg: polynomial to be read
**************************************************/
void poly_q_set(poly_q res, const poly_q arg) {
	nmod_poly_set(res, arg);
}

/*************************************************
* Name:        poly_q_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
*
* Arguments:   - const poly_q arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
*
* Returns the coefficients of x^n of arg
**************************************************/
coeff_q poly_q_get_coeff(const poly_q arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot get coefficient of poly at given position.");
  return nmod_poly_get_coeff_ui(arg, n);
}

/*************************************************
* Name:        poly_q_get_coeff_centered
*
* Description: Get coefficient of x^n of a polynomial in centered representation
*              condition: [0 <= n < PARAM_N]
*
* Arguments:   - const poly_q arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
*
* Returns the coefficients of x^n of arg in [-PARAM_Q/2, PARAM_Q/2]
**************************************************/
coeff_q poly_q_get_coeff_centered(const poly_q arg, size_t n) {
  coeff_q tmp = poly_q_get_coeff(arg, n);
  return tmp - ((~((tmp - PARAM_Q/2) >> (sizeof(coeff_q)*8-1))) & PARAM_Q);
}

/*************************************************
* Name:        poly_q_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
*              Coefficient is reduced mod PARAM_Q
*
* Arguments:   - poly_q arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_q c: the new coefficient
**************************************************/
void poly_q_set_coeff(poly_q arg, size_t n, coeff_q c) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot set coefficient of poly at given position.");
  c %= PARAM_Q;
  c += (c >> (sizeof(coeff_q)*8-1)) & PARAM_Q;
  ASSERT_DEBUG(c >= 0, "Don't give me too much negativity!");
  nmod_poly_set_coeff_ui(arg, n, c);
}

/*************************************************
* Name:        poly_q_from_bits
*
* Description: Set polynomial with {0,1} coefficients from byte-array
*              Example: Given an array { 0xB2, 0x08 } where 0xB2 = 0b1011 0010
*                   and 0x08 = 0b0000 1000 and a polynomial p with precisely 16
*                   coefficients, this functions sets p = x^15 + x^13 + x^12 + x^9 + x^3.
*
* Arguments:   - poly_q arg: polynomial to be set from byte-array (initialized)
*              - const uint8_t *coeffs: byte array for coefficients (allocated PARAM_N/8 bytes)
**************************************************/
void poly_q_from_bits(poly_q arg, const uint8_t coeffs[PARAM_N / 8]) {
  poly_q_zero(arg);
  for (size_t i = 0; i < PARAM_N; ++i) {
    coeff_q c = coeffs[i/8];
    c = (c >> (i % 8)) & 1;
    poly_q_set_coeff(arg, i, c);
  }
}

/*************************************************
* Name:        poly_q_neg
*
* Description: Negate a polynomial coefficient-wise
*              Coefficients are reduced mod PARAM_Q
*
* Arguments:   - poly_q res: polynomial to host the negation (initialized)
*              - const poly_q arg: polynomial to be negated
**************************************************/
void poly_q_neg(poly_q res, const poly_q arg) {
  nmod_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_q_add
*
* Description: Add two polynomials. Coefficients are reduced mod PARAM_Q
*
* Arguments:   - poly_q res: polynomial to host the sum (initialized)
*              - const poly_q lhs: first polynomial summand
*              - const poly_q rhs: second polynomial summand
**************************************************/
void poly_q_add(poly_q res, const poly_q lhs, const poly_q rhs) {
	nmod_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_q_sub
*
* Description: Substract two polynomials. Coefficients are reduced mod PARAM_Q
*
* Arguments:   - poly_q res: polynomial to host the difference (initialized)
*              - const poly_q lhs: first polynomial term
*              - const poly_q rhs: second polynomial term
**************************************************/
void poly_q_sub(poly_q res, const poly_q lhs, const poly_q rhs) {
  nmod_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_q_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N + 1.
*              Coefficients are reduced mod PARAM_Q
*
* Arguments:   - poly_q res: polynomial to host the multiplication (initialized)
*              - const poly_q lhs: first polynomial factor
*              - const poly_q rhs: second polynomial factor
**************************************************/
void poly_q_mul(poly_q res, const poly_q lhs, const poly_q rhs) {
  ASSERT_DEBUG(nmod_poly_degree(lhs) < PARAM_N, "Argument to `poly_q_mul` must already be reduced for FLINT.");
  ASSERT_DEBUG(nmod_poly_degree(rhs) < PARAM_N, "Argument to `poly_q_mul` must already be reduced for FLINT.");
	nmod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_q_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
*
* Arguments:   - poly_q res: polynomial to host the multiplication (initialized)
*              - const poly_q arg: first polynomial factor
*              - coeff_q fac: second scalar factor
**************************************************/
void poly_q_mul_scalar(poly_q res, const poly_q arg, coeff_q fac) {
  ASSERT_DEBUG(fac < PARAM_Q, "Scalar must not exceed q for scalar polynomial multiplication.");
  nmod_poly_scalar_mul_nmod(res, arg, fac);
}

/*************************************************
* Name:        poly_q_shift_left
*
* Description: Shift the coefficients of a polynomial to
*              the left by n places. Corresponds to a
*              multiplication by x^n (NOT REDUCED MOD X^PARAM_N + 1)
*              Trailing zeros are inserted
*
* Arguments:   - poly_q res: polynomial to host the shift (initialized)
*              - const poly_q arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_q_shift_left(poly_q res, const poly_q arg, size_t n) {
  nmod_poly_shift_left(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_q_shift_right
*
* Description: Shift the coefficients of a polynomial to
*              the right by n places. Corresponds to a
*              division by x^n (NOT REDUCED MOD X^PARAM_N + 1)
*              Leading zeros are inserted
*
* Arguments:   - poly_q res: polynomial to host the shift (initialized)
*              - const poly_q arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_q_shift_right(poly_q res, const poly_q arg, size_t n) {
  nmod_poly_shift_right(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_q_conjugate
*
* Description: Compute the conjugate of a polynomial (evaluation at x^-1)
*
* Arguments:   - poly_q res: polynomial to host the conjugate (initialized)
*              - const poly_q arg: polynomial to be conjugated
**************************************************/
void poly_q_conjugate(poly_q res, const poly_q arg) {
	coeff_q c = poly_q_get_coeff(arg,0);
	poly_q_set_coeff(res, 0, c);
	for (size_t i = 1; i < PARAM_N; ++i) {
		c = poly_q_get_coeff(arg, PARAM_N - i);
    poly_q_set_coeff(res, i, -c);
  }
}

/*************************************************
* Name:        poly_q_equal
*
* Description: Equality test between two polynomials
*
* Arguments:   - const poly_q lhs: first polynomial
*              - const poly_q rhs: second polynomial
*
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_q_equal(const poly_q lhs, const poly_q rhs) {
  return nmod_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_q_dump
*
* Description: Print a polynomial
*
* Arguments:   - const poly_q arg: the polynomial to be printed
**************************************************/
void poly_q_dump(const poly_q arg) {
  printf("[");
  for (size_t i = 0; i < PARAM_N; i++) {
    printf("%ld, ", poly_q_get_coeff_centered(arg, i));
  }
  printf("]");
}

/*************************************************
* Name:        poly_q_sq_norm2
*
* Description: Compute the square l2 norm of a polynomial
*
* Arguments:   - const poly_q arg: the polynomial
*
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_q_sq_norm2(const poly_q arg) {
	uint64_t sq_norm2 = 0;
	coeff_q c;
	for (size_t i = 0; i < PARAM_N; i++) {
		c = poly_q_get_coeff_centered(arg, i);
		CHK_UI_OVF_ADDITION(sq_norm2, (uint64_t) (c * c));
	}
	return sq_norm2;
}

/*************************************************
* Name:        poly_q_weight
*
* Description: Compute the Hamming weight of a binary
*              polynomial, fails if not binary
*
* Arguments:   - const poly_q arg: the polynomial
*
* Returns the Hamming weight of the binary polynomial, and -1 if not binary
**************************************************/
int64_t poly_q_weight(const poly_q arg) {
  int64_t weight = 0;
  coeff_q tmp;
  for (size_t i = 0; i < PARAM_N; i++) {
  	tmp = poly_q_get_coeff(arg, i);
  	if (tmp >> 1 != 0) {
  		return -1;
  	}
  	weight += tmp;
  }
  return weight;
}

/*************************************************
* Name:        poly_q_invert
*
* Description: Compute the inverse of a polynomial in Z[x]/(q, x^n+1)
*
* Arguments:   - poly_q res: polynomial to host the inverse mod PARAM_Q
*              - const poly_q arg: polynomial to invert
**************************************************/
void poly_q_invert(poly_q res, const poly_q arg) {
  nmod_poly_xgcd(TMP, AUX, res, POLY_F, arg);
  ASSERT_DEBUG(nmod_poly_is_one(TMP), "Polynomial inversion did not work correctly.");
}
