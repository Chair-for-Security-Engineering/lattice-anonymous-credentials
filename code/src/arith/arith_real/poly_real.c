#include "poly_real.h"
#include "random.h"
#include "macros.h"

#include <flint/acb_dft.h>
#include <flint/arb_poly.h>

#if PARAM_N != 256
#error "SampleFz is implemented statically for n = 256."
#endif

static acb_dft_pre_t pre[8];

static arb_t TMP;
static arb_t AUX;

/*************************************************
* Name:        arith_real_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              multiprecision reals. This is strictly required 
*              and must be called once before any other function 
*              from here is used.
**************************************************/
void arith_real_setup(void) {
  for (size_t i = 1; i <= 8; i++) {
    acb_dft_precomp_init(pre[i-1], 2u<<i, POLY_REAL_PRECISION);
  }
  arb_init(TMP);
  arb_init(AUX);
}

/*************************************************
* Name:        arith_real_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              multiprecision reals. This is strictly required 
*              and must be called once at the very end to release 
*              any resources.
**************************************************/
void arith_real_teardown(void) {
  for (size_t i = 0; i < 8; i++)
  {
    acb_dft_precomp_clear(pre[i]);
  }
  arb_clear(TMP);
  arb_clear(AUX);
}

/*************************************************
* Name:        multiply [static]
*
* Description: Multiplication of two polynomials reduced
*              modulo x^(2^log_dim) + 1. Computation in 
*              a subring as used in SampleFz [GM18]
* 
* Arguments:   - arb_poly_t out: polynomial to host the multiplication (initialized)
*              - const arb_poly_t x: first polynomial factor
*              - const arb_poly_t y: second polynomial factor
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void multiply(arb_poly_t out, const arb_poly_t x, const arb_poly_t y, const size_t log_dim) {
  if (log_dim > 0) {
    signed long i;
    arb_poly_mul(out, x, y, POLY_REAL_PRECISION);
    // reduce
    for (i = 0; i < 1u << log_dim; i++) {
      arb_poly_get_coeff_arb(TMP,  out, i);
      arb_poly_get_coeff_arb(AUX, out, i + (1<<log_dim));
      arb_sub(TMP, TMP, AUX, POLY_REAL_PRECISION);
      arb_poly_set_coeff_arb(out, i, TMP);
      arb_poly_set_coeff_si( out, i + (1<<log_dim), 0);
    }
  } else {
    arb_poly_get_coeff_arb(TMP, x, 0);
    arb_poly_get_coeff_arb(AUX, y, 0);
    arb_mul(TMP, TMP, AUX, POLY_REAL_PRECISION);
    arb_poly_set_coeff_arb(out, 0, TMP);
  }
}

/*************************************************
* Name:        invert [static]
*
* Description: Invert a polynomial modulo x^(2^log_dim) + 1. 
*              Computation in a subring as used in SampleFz [GM18]
* 
* Arguments:   - arb_poly_t out: polynomial to host the inverse (initialized)
*              - const arb_poly_t x: polynomial to be inverted
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void invert(arb_poly_t out, const arb_poly_t x, const size_t log_dim) {
  if (log_dim > 0) {
    signed long i;
    acb_struct vec[512];
    //acb_t debugtmp;
    arb_t resq, imsq;
    arb_init(resq);
    arb_init(imsq);
    //acb_init(debugtmp);
    for (i = 0; i < 2u << log_dim; i++) {
      acb_init(&vec[i]);
    }
    for (i = 0; i < 1u << log_dim; i++) {
      arb_poly_get_coeff_arb(TMP, x, i);
      acb_set_arb(&vec[i], TMP);
    }
    acb_dft_precomp(vec, vec, pre[log_dim-1], POLY_REAL_PRECISION);
    for (i = 0; i < (2u << log_dim); i++) { // invert
      //acb_set(debugtmp, &vec[i]);
      acb_inv(&vec[i], &vec[i], POLY_REAL_PRECISION);
      /*if (!acb_is_finite(&vec[i]))
      {
        printf("\n");
        acb_print(debugtmp);
        fflush(stdout);
        acb_inv(debugtmp, debugtmp, POLY_REAL_PRECISION);
        printf("\n");
        acb_print(debugtmp);
        printf("\n");
        acb_print(&vec[i]);
        fflush(stdout);
      }*/
    }
    acb_dft_inverse_precomp(vec, vec, pre[log_dim-1], POLY_REAL_PRECISION);
    // reduce and write result
    for (i = 0; i < 1u << log_dim; i++) {
      acb_get_real(TMP, &vec[i]);
      acb_get_real(AUX, &vec[i+(1u<<log_dim)]);
      //if(!arb_is_finite(TMP)) {printf("low not finite\n");arb_print(TMP);acb_print(&vec[i]);}
      //if(!arb_is_finite(AUX)) {printf("high not finite\n");arb_print(AUX);acb_print(&vec[i+(1u<<log_dim)]);}
      arb_sub(TMP, TMP, AUX, POLY_REAL_PRECISION);
      arb_poly_set_coeff_arb(out, i, TMP);
    }

#if 0
    ///////////////////// TEST
    if (out != x)
    {
      arb_poly_t res;
      arb_poly_init(res);
      multiply(res, out, x, log_dim);
      printf("\n\nthis should be one: ");
      arb_poly_printd(res, 5);
      printf("\n");
      fflush(stdout);
      arb_poly_clear(res);
    } else {
      printf("out == x.");
    }
    ///////////////////// END OF TEST
#endif
    // clean up flint objects
    for (i = 0; i < 2u << log_dim; i++) {
      acb_clear(&vec[i]);
    }
    arb_clear(resq);
    arb_clear(imsq);
    //acb_clear(debugtmp);
  } else {
    arb_poly_get_coeff_arb(TMP, x, 0);
    ASSERT_DEBUG(!arb_is_zero(TMP) && arb_is_finite(TMP), "Unexpected zero value");
    arb_inv(TMP, TMP, POLY_REAL_PRECISION);
    ASSERT_DEBUG(!arb_is_zero(TMP) && arb_is_finite(TMP), "Unexpected zero value");
    arb_poly_set_coeff_arb(out, 0, TMP);
  }
}

/*************************************************
* Name:        conjugate [static]
*
* Description: Conjugate a polynomial modulo x^(2^log_dim) + 1. 
*              Computation in a subring as used in SampleFz [GM18]
* 
* Arguments:   - arb_poly_t x: polynomial to be conjugated [in place] (initialized)
*              - size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void conjugate(arb_poly_t x, size_t log_dim) {
  if (log_dim > 0) {
    signed long i;
    for (i = 1; i < 1u << (log_dim-1); i++) {
      arb_poly_get_coeff_arb(TMP,  x, i);
      arb_poly_get_coeff_arb(AUX, x, (1<<log_dim)-i);

      arb_neg(TMP, TMP);
      arb_neg(AUX, AUX);

      arb_poly_set_coeff_arb(x, i, AUX);
      arb_poly_set_coeff_arb(x, (1<<log_dim)-i, TMP);
    }
    arb_poly_get_coeff_arb(TMP, x, 1<<(log_dim-1));
    arb_neg(TMP, TMP);
    arb_poly_set_coeff_arb(x, 1<<(log_dim-1), TMP);
  }
}

/*************************************************
* Name:        split [static]
*
* Description: Split a polynomial of degree n into two polynomials
*              of degree n/2 by splitting the even and odd exponents. 
* 
* Arguments:   - arb_poly_t even: polynomial to host the even exponents' coefficients (initialized)
*              - arb_poly_t odd: polynomial to host the odd exponents' coefficients (initialized)
*              - const arb_poly_t x: polynomial to be splitted
**************************************************/
static inline void split(arb_poly_t even, arb_poly_t odd, const arb_poly_t x) {
  for (signed long i = 0; i < arb_poly_length(x); i++) {
    arb_poly_get_coeff_arb(TMP, x, i);
    if ((i%2) == 0) {
      arb_poly_set_coeff_arb(even, i/2, TMP);
    } else {
      arb_poly_set_coeff_arb(odd, i/2, TMP);
    }
  }
}

/*************************************************
* Name:        _samplefz [static]
*
* Description: Sample a Gaussian element of Z[x]/(x^(2^log_dim) + 1)
*              of covariance real-valued polynomial f and center 
*              real-valued polynomial c, using the [GM18] SampleFz 
*              algorithm. Computation in a subring as used in [GM18]
* 
* Arguments:   - arb_poly_t res: polynomial to host the Gaussian sample (initialized)
*              - const arb_poly_t f: covariance polynomial 
*              - const arb_poly_t f: center polynomial 
*              - size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static void _samplefz(arb_poly_t res, const arb_poly_t f, const arb_poly_t c, const size_t log_dim) {
  if (log_dim == 0) {
    arb_t tmpc, tmpf;
    arb_init(tmpc);
    arb_init(tmpf);
    arb_poly_get_coeff_arb(tmpc, c, 0);
    arb_poly_get_coeff_arb(tmpf, f, 0);
    arb_sqrt(tmpf, tmpf, POLY_REAL_PRECISION);
    ASSERT_DEBUG(arb_is_finite(tmpc), "Unexpected infinite value");
    ASSERT_DEBUG(arb_is_finite(tmpf), "Unexpected infinite value");
    arb_set_si(tmpf, SampleZ(
      arf_get_d(arb_midref(tmpc), ARF_RND_DOWN),
      arf_get_d(arb_midref(tmpf), ARF_RND_DOWN)
    ));
    arb_poly_set_coeff_arb(res, 0, tmpf);
    arb_clear(tmpf);
    arb_clear(tmpc);
  } else {
    arb_poly_t f_even, f_odd, c_even, c_odd, q_even, q_odd, tmp;
    arb_t tcf;

    // init flint objects
    arb_poly_init(f_even);
    arb_poly_init(f_odd);
    arb_poly_init(c_even);
    arb_poly_init(c_odd);
    arb_poly_init(q_even);
    arb_poly_init(q_odd);
    arb_poly_init(tmp);
    arb_init(tcf);

    // sampling begins
    split(f_even, f_odd, f);
    split(c_even, c_odd, c);

    _samplefz(q_odd, f_even, c_odd, log_dim-1);

    invert(tmp, f_even, log_dim-1);
    multiply(tmp, tmp, f_odd, log_dim-1);
    arb_poly_sub(c_odd, q_odd, c_odd, POLY_REAL_PRECISION); // last usage of c_odd, this is now a temporary variable
    multiply(c_odd, tmp, c_odd, log_dim-1);
    arb_poly_add(c_even, c_even, c_odd, POLY_REAL_PRECISION);

    conjugate(f_odd, log_dim-1);
    multiply(tmp, tmp, f_odd, log_dim-1);
    arb_poly_sub(tmp, f_even, tmp, POLY_REAL_PRECISION); // last usage of f_even

    _samplefz(q_even, tmp, c_even, log_dim-1);

    // merge
    for (signed long i = 0; i < 1u << log_dim; i++) {
      if ((i%2) == 0) {
        arb_poly_get_coeff_arb(tcf, q_even, i/2);
      } else {
        arb_poly_get_coeff_arb(tcf, q_odd, i/2);
      }
      arb_poly_set_coeff_arb(res, i, tcf);
    }

    // clean up flint objects
    arb_poly_clear(f_even);
    arb_poly_clear(f_odd);
    arb_poly_clear(c_even);
    arb_poly_clear(c_odd);
    arb_poly_clear(q_even);
    arb_poly_clear(q_odd);
    arb_poly_clear(tmp);
    arb_clear(tcf);
  }
}

/*************************************************
* Name:        poly_real_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_real res: polynomial to be initialized
**************************************************/
void poly_real_init(poly_real res) {
  arb_poly_init(res);
}

/*************************************************
* Name:        poly_real_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_real res: polynomial to be cleared
**************************************************/
void poly_real_clear(poly_real res) {
  arb_poly_clear(res);
}

/*************************************************
* Name:        poly_real_mul_scalar
*
* Description: Multiplication of a polynomial by a scalar
* 
* Arguments:   - poly_real res: polynomial to host the multiplication (initialized)
*              - const poly_real arg: first polynomial factor
*              - const coeff_real f: second scalar factor
**************************************************/
void poly_real_mul_scalar(poly_real res, const poly_real arg, const coeff_real f) {
  arb_set_d(TMP, f);
  arb_poly_scalar_mul(res, arg, TMP, POLY_REAL_PRECISION);
}

/*************************************************
* Name:        poly_real_add_constant
*
* Description: Add a scalar to the constant coefficient of a polynomial
* 
* Arguments:   - poly_real res: polynomial to be incremented (initialized)
*              - const coeff_real c: scalar summand
**************************************************/
void poly_real_add_constant(poly_real res, const coeff_real c) {
  arb_poly_get_coeff_arb(TMP, res, 0);
  arb_set_d(AUX, c);
  arb_add(TMP, TMP, AUX, POLY_REAL_PRECISION);
  arb_poly_set_coeff_arb(res, 0, TMP);
}

/*************************************************
* Name:        poly_real_invert
*
* Description: Invert a real-valued polynomial modulo
*              x^PARAM_N + 1
* 
* Arguments:   - poly_real res: polynomial to host the inverse (initialized)
*              - const poly_real arg: polynomial to be inverted
**************************************************/
void poly_real_invert(poly_real res, const poly_real arg) {
  invert(res, arg, 8);
}

/*************************************************
* Name:        poly_real_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N + 1
* 
* Arguments:   - poly_real res: polynomial to host the multiplication (initialized)
*              - const poly_real lhs: first polynomial factor
*              - const poly_real rhs: second polynomial factor
**************************************************/
void poly_real_mul(poly_real res, const poly_real lhs, const poly_real rhs) {
  multiply(res, lhs, rhs, 8);
}

/*************************************************
* Name:        poly_real_add
*
* Description: Add two polynomials 
* 
* Arguments:   - poly_real res: polynomial to host the sum (initialized)
*              - const poly_real lhs: first polynomial summand
*              - const poly_real rhs: second polynomial summand
**************************************************/
void poly_real_add(poly_real res, const poly_real lhs, const poly_real rhs) {
  arb_poly_add(res, lhs, rhs, POLY_REAL_PRECISION);
}

/*************************************************
* Name:        poly_real_sub
*
* Description: Substract two polynomials 
* 
* Arguments:   - poly_real res: polynomial to host the difference (initialized)
*              - const poly_real lhs: first polynomial term
*              - const poly_real rhs: second polynomial term
**************************************************/
void poly_real_sub(poly_real res, const poly_real lhs, const poly_real rhs) {
  arb_poly_sub(res, lhs, rhs, POLY_REAL_PRECISION);
}

/*************************************************
* Name:        poly_real_set_si
*
* Description: Set coefficient of x^n of a polynomial
* 
* Arguments:   - poly_real arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - int64_t c: the new signed 64-bit int coefficient
**************************************************/
void poly_real_set_si(poly_real res, size_t n, int64_t c) {
  arb_set_si(TMP, c);
  arb_poly_set_coeff_arb(res, (signed long)n, TMP);
}

/*************************************************
* Name:        poly_real_get_si
*
* Description: Get rounded coefficient of x^n of a polynomial
* 
* Arguments:   - const poly_real arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*
* Returns the signed 64-bit integer rounded coefficient of x^n
**************************************************/
int64_t poly_real_get_si(const poly_real arg, size_t n) {
  int64_t res;
  arb_poly_get_coeff_arb(TMP, arg, (long int) n);
  res = arf_get_si(arb_midref(TMP), ARF_RND_DOWN);
  return res;
}

/*************************************************
* Name:        poly_real_samplefz
*
* Description: Sample a Gaussian element of Z[x]/(x^PARAM_N + 1)
*              of covariance real-valued polynomial f and center 
*              real-valued polynomial c, using the [GM18] SampleFz 
*              algorithm
* 
* Arguments:   - poly_real res: polynomial to host the Gaussian sample (initialized)
*              - const poly_real f: covariance polynomial 
*              - const poly_real f: center polynomial
**************************************************/
void poly_real_samplefz(poly_real res, const poly_real f, const poly_real c) {
  _samplefz(res, f, c, 8);
}

/*************************************************
* Name:        poly_real_get_coeff_rounded
*
* Description: Get rounded coefficient of x^n of a polynomial
* 
* Arguments:   - const poly_real arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*
* Returns the signed 64-bit integer rounded coefficient of x^n
**************************************************/
int64_t poly_real_get_coeff_rounded(const poly_real arg, size_t n) {
  int64_t ret;
  arb_poly_get_coeff_arb(TMP, arg, (signed long)n);
  ret = arf_get_si(arb_midref(TMP), ARF_RND_DOWN);
  return ret;
}

/*************************************************
* Name:        poly_real_dump
*
* Description: Print a polynomial
* 
* Arguments:   - const poly_real arg: the polynomial to be printed
**************************************************/
void poly_real_dump(const poly_real arg) {
  arb_poly_printd(arg, 5);
  printf("\n");
}
