#include "arith_q.h"
#include "arith_qiss.h"
#include "macros.h"

static nmod_poly_t POLY_F;
static nmod_poly_t POLY_F_REV_INV;

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
* Name:        arith_qiss_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              modulo PARAM_Q_ISS. This is strictly required 
*              and must be called once before any other function 
*              from here is used.
**************************************************/
void arith_qiss_setup(void) {
	nmod_poly_init(POLY_F, PARAM_Q_ISS);
	nmod_poly_set_coeff_ui(POLY_F, 0, 1);
	nmod_poly_set_coeff_ui(POLY_F, PARAM_N_ISS, 1);

  nmod_poly_t f_len;
  nmod_poly_init(f_len, PARAM_Q_ISS);
  nmod_poly_set_coeff_ui(f_len, PARAM_N_ISS + 1, 1);

  nmod_poly_t f_rev;
  nmod_poly_init(f_rev, PARAM_Q_ISS);

  nmod_poly_init(POLY_F_REV_INV, PARAM_Q_ISS);
  nmod_poly_reverse(f_rev, POLY_F, PARAM_N_ISS + 1);
  nmod_poly_invert(POLY_F_REV_INV, f_rev, f_len, PARAM_Q_ISS);

  nmod_poly_clear(f_len);
  nmod_poly_clear(f_rev);
}

/*************************************************
* Name:        arith_qiss_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              modulo PARAM_Q_ISS. This is strictly required 
*              and must be called once at the very end to release 
*              any resources.
**************************************************/
void arith_qiss_teardown(void) {
  nmod_poly_clear(POLY_F);
  nmod_poly_clear(POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qiss_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_qiss res: polynomial to be initialized
**************************************************/
void poly_qiss_init(poly_qiss res) {
  nmod_poly_init(res, PARAM_Q_ISS);
}

/*************************************************
* Name:        poly_qiss_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qiss arg: polynomial to be cleared
**************************************************/
void poly_qiss_clear(poly_qiss arg) {
  nmod_poly_clear(arg);
}

/*************************************************
* Name:        poly_qiss_zero
*
* Description: Set an initialized polynomial to zero
* 
* Arguments:   - poly_qiss res: polynomial to be zeroized (initialized)
**************************************************/
void poly_qiss_zero(poly_qiss res) {
  nmod_poly_zero(res);
}

/*************************************************
* Name:        poly_qiss_set
*
* Description: Set a polynomial equal to another polynomial.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to be set (initialized)
*              - const poly_qiss arg: polynomial to be read
**************************************************/
void poly_qiss_set(poly_qiss res, const poly_qiss arg) {
	nmod_poly_set(res, arg);
}

/*************************************************
* Name:        poly_qiss_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_ISS]
* 
* Arguments:   - const poly_qiss arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg
**************************************************/
coeff_qiss poly_qiss_get_coeff(const poly_qiss arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N_ISS, "Illegal argument: cannot get coefficient of poly at given position.");
	return nmod_poly_get_coeff_ui(arg, n);
}

/*************************************************
* Name:        poly_qiss_get_coeff_centered
*
* Description: Get coefficient of x^n of a polynomial in centered representation
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - const poly_qiss arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg in [-PARAM_Q_ISS/2, PARAM_Q_ISS/2]
**************************************************/
coeff_qiss poly_qiss_get_coeff_centered(const poly_qiss arg, size_t n) {
	coeff_qiss tmp = poly_qiss_get_coeff(arg, n);
	return tmp - ((~((tmp - PARAM_Q_ISS/2) >> (sizeof(coeff_qiss)*8-1))) & PARAM_Q_ISS);
}

/*************************************************
* Name:        poly_qiss_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_ISS]
*              Coefficient is reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_qiss c: the new coefficient
**************************************************/
void poly_qiss_set_coeff(poly_qiss arg, size_t n, coeff_qiss c) {
  ASSERT_DEBUG(n < PARAM_N_ISS, "Illegal argument: cannot set coefficient of poly at given position.");
	c += PARAM_Q_ISS;
	c %= PARAM_Q_ISS;
  c += (c >> (sizeof(coeff_qiss)*8 - 1)) & PARAM_Q_ISS;
  ASSERT_DEBUG(c >= 0, "Coefficients of polys must not be negative.");
  nmod_poly_set_coeff_ui(arg, n, c);
}

/*************************************************
* Name:        poly_qiss_neg
*
* Description: Negate a polynomial coefficient-wise
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the negation (initialized)
*              - const poly_qiss arg: polynomial to be negated
**************************************************/
void poly_qiss_neg(poly_qiss res, const poly_qiss arg) {
  nmod_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_qiss_add
*
* Description: Add two polynomials. Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the sum (initialized)
*              - const poly_qiss lhs: first polynomial summand
*              - const poly_qiss rhs: second polynomial summand
**************************************************/
void poly_qiss_add(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
	nmod_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_sub
*
* Description: Substract two polynomials. Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the difference (initialized)
*              - const poly_qiss lhs: first polynomial term
*              - const poly_qiss rhs: second polynomial term
**************************************************/
void poly_qiss_sub(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
	nmod_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N_ISS + 1.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss lhs: first polynomial factor
*              - const poly_qiss rhs: second polynomial factor
**************************************************/
void poly_qiss_mul(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
  ASSERT_DEBUG(nmod_poly_degree(lhs) < PARAM_N_ISS, "Argument to `poly_qiss_mul` must already be reduced for FLINT.");
  ASSERT_DEBUG(nmod_poly_degree(rhs) < PARAM_N_ISS, "Argument to `poly_qiss_mul` must already be reduced for FLINT.");
	nmod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qiss_mul_x
*
* Description: Multiplication of a polynomial by x and reduced
*              modulo x^PARAM_N_ISS + 1.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: polynomial to be multiplied by x
**************************************************/
void poly_qiss_mul_x(poly_qiss res, const poly_qiss arg) {
  coeff_qiss c = poly_qiss_get_coeff(arg, PARAM_N_ISS - 1);
  poly_qiss_shift_left(res, arg, 1);
  poly_qiss_set_coeff(res, 0, -c);
  // removing leading coefficient manually as poly_qiss_set_coeff only authorizes exponent < PARAM_N_ISS
  nmod_poly_set_coeff_ui(res, PARAM_N_ISS, 0);
}

/*************************************************
* Name:        poly_qiss_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: first polynomial factor
*              - coeff_qiss fac: second scalar factor
**************************************************/
void poly_qiss_mul_scalar(poly_qiss out, const poly_qiss lhs, const coeff_qiss rhs) {
  coeff_qiss f = rhs % PARAM_Q_ISS;
  f += (f >> (sizeof(coeff_qiss)*8 - 1)) & PARAM_Q_ISS;
  ASSERT_DEBUG(f >= 0, "Factors for scalar multiplication must not be negative.");
  nmod_poly_scalar_mul_nmod(out, lhs, f);
}

/*************************************************
* Name:        poly_qiss_shift_left
*
* Description: Shift the coefficients of a polynomial to
*              the left by n places. Corresponds to a 
*              multiplication by x^n (NOT REDUCED MOD X^PARAM_N_ISS + 1)
*              Trailing zeros are inserted
* 
* Arguments:   - poly_qiss res: polynomial to host the shift (initialized)
*              - const poly_qiss arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_qiss_shift_left(poly_qiss res, const poly_qiss arg, size_t n) {
  nmod_poly_shift_left(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_qiss_conjugate
*
* Description: Compute the conjugate of a polynomial (evaluation at x^-1)
* 
* Arguments:   - poly_qiss res: polynomial to host the conjugate (initialized)
*              - const poly_qiss arg: polynomial to be conjugated
**************************************************/
void poly_qiss_conjugate(poly_qiss res, const poly_qiss arg) {
  ASSERT_DEBUG(res != arg, "Input and output of conjugation must not be idential.");
  coeff_qiss c = poly_qiss_get_coeff(arg, 0);
  poly_qiss_set_coeff(res, 0, c);
  for (size_t i = 1; i < PARAM_N_ISS; ++i) {
    c = poly_qiss_get_coeff_centered(arg, PARAM_N_ISS - i);
    poly_qiss_set_coeff(res, i, -c);
  }
}

/*************************************************
* Name:        poly_qiss_equal
*
* Description: Equality test between two polynomials
* 
* Arguments:   - const poly_qiss lhs: first polynomial
*              - const poly_qiss rhs: second polynomial
* 
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_qiss_equal(const poly_qiss lhs, const poly_qiss rhs) {
  return nmod_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_dump
*
* Description: Print a polynomial
* 
* Arguments:   - const poly_qiss arg: the polynomial to be printed
**************************************************/
void poly_qiss_dump(const poly_qiss arg) {
	nmod_poly_print(arg);
}

/*************************************************
* Name:        poly_qiss_sq_norm2
*
* Description: Compute the square l2 norm of a polynomial
* 
* Arguments:   - const poly_qiss arg: the polynomial
* 
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_qiss_sq_norm2(const poly_qiss arg) {
	uint64_t sq_norm2 = 0;
	coeff_qiss c;
	for (size_t i = 0; i < PARAM_N_ISS; i++) {
		c = poly_qiss_get_coeff_centered(arg, i);
    CHK_UI_OVF_ADDITION(sq_norm2, (uint64_t) (c * c));
	}
	return sq_norm2;
}

/*************************************************
* Name:        poly_qiss_pack
*
* Description: Pack a polynomial mod PARAM_Q_ISS into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQISS_PACKEDBYTES bytes)
*              - const poly_qiss arg: the polynomial to be packed
**************************************************/
void poly_qiss_pack(uint8_t buf[POLYQISS_PACKEDBYTES], const poly_qiss arg) {
  uint64_t x;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
	  x = poly_qiss_get_coeff(arg, i);
#if PARAM_Q_ISS_BITLEN > 40
#error "poly_qiss_pack assumes that PARAM_Q_ISS_BITLEN <= 40"
#endif
#if PARAM_Q1_ISS_BITLEN > 20
#error "poly_qiss_pack assumes that PARAM_Q1_ISS_BITLEN <= 20"
#endif
    buf[5*i + 0] = x >>  0;
    buf[5*i + 1] = x >>  8;
    buf[5*i + 2] = x >> 16;
    buf[5*i + 3] = x >> 24;
    buf[5*i + 4] = x >> 32;
  }
}

/*************************************************
* Name:        coeff_qiss_pack
*
* Description: Pack a polynomial coefficient mod PARAM_Q_ISS into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated COEFFQISS_PACKEDBYTES bytes)
*              - const coeff_qiss arg: the coefficient to be packed
**************************************************/
void coeff_qiss_pack(uint8_t buf[COEFFQISS_PACKEDBYTES], const coeff_qiss arg) {
  uint64_t x;
  coeff_qiss tmp = arg + ((arg >> (sizeof(coeff_qiss)*8-1)) & PARAM_Q_ISS);
  x = tmp % PARAM_Q1_ISS; // first CRT slot on first 20 bits
  x |= (tmp % PARAM_Q2_ISS) << 20; // second CRT slot on next 20 bits
#if PARAM_Q_ISS_BITLEN > 40
#error "poly_qiss_pack assumes that PARAM_Q_ISS_BITLEN <= 40"
#endif
#if PARAM_Q1_ISS_BITLEN > 20
#error "poly_qiss_pack assumes that PARAM_Q1_ISS_BITLEN <= 20"
#endif
  buf[0] = x >>  0;
  buf[1] = x >>  8;
  buf[2] = x >> 16;
  buf[3] = x >> 24;
  buf[4] = x >> 32;
}
