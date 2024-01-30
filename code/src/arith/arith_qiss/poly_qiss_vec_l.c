#include "arith_qiss.h"
#include "macros.h"

static poly_qiss TMP;

/*************************************************
* Name:        poly_qiss_vec_l_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_ISS integer vectors with PARAM_L_ISS entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_qiss_vec_l_setup(void) {
  poly_qiss_init(TMP);
}

/*************************************************
* Name:        poly_qiss_vec_l_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_ISS integer vectors with PARAM_L_ISS entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_qiss_vec_l_teardown(void) {
  poly_qiss_clear(TMP);
}

/*************************************************
* Name:        poly_qiss_vec_l_init
*
* Description: Initialize polynomial vector with PARAM_L_ISS entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_qiss_vec_l arg: polynomial vector to be initialized
**************************************************/
void poly_qiss_vec_l_init(poly_qiss_vec_l arg) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_clear
*
* Description: Clear polynomial vector with PARAM_L_ISS entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qiss_vec_l arg: polynomial vector to be cleared
**************************************************/
void poly_qiss_vec_l_clear(poly_qiss_vec_l arg) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_zero
*
* Description: Set an initialized polynomial vector with PARAM_L_ISS entries to zero
* 
* Arguments:   - poly_qiss_vec_l arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_qiss_vec_l_zero(poly_qiss_vec_l arg) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_set
*
* Description: Set a polynomial vector with PARAM_L_ISS entries equal to another polynomial vector
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to be set (initialized)
* 			   		 - const poly_qiss_vec_l arg: polynomial vector to be read
**************************************************/
void poly_qiss_vec_l_set(poly_qiss_vec_l res, const poly_qiss_vec_l arg) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_L_ISS]
* 
* Arguments:   - poly_qiss res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_qiss_vec_l arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qiss_vec_l_get_poly(poly_qiss res, const poly_qiss_vec_l arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_L_ISS, "Illegal argument: cannot get entry of vector at given position.");
	poly_qiss_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_qiss_vec_l_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_L_ISS]
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to be set (initialized)
* 			   		 - const poly_qiss arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qiss_vec_l_set_poly(poly_qiss_vec_l res, const poly_qiss arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_L_ISS, "Illegal argument: cannot set entry of vector at given position.");
	poly_qiss_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_qiss_vec_l_neg
*
* Description: Negate a polynomial vector with PARAM_L_ISS entries
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_qiss_vec_l arg: polynomial vector to be negated
**************************************************/
void poly_qiss_vec_l_neg(poly_qiss_vec_l res, const poly_qiss_vec_l arg) {
  for (size_t i = 0; i < PARAM_L_ISS; ++i) {
    poly_qiss_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qiss_vec_l_add
*
* Description: Add two polynomial vectors with PARAM_L_ISS entries
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_qiss_vec_l lhs: first polynomial vector summand
* 			   		 - const poly_qiss_vec_l rhs: second polynomial vector summand
**************************************************/
void poly_qiss_vec_l_add(poly_qiss_vec_l res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_sub
*
* Description: Substract two polynomial vectors with PARAM_L_ISS entries
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_qiss_vec_l lhs: first polynomial vector term
* 			   		 - const poly_qiss_vec_l rhs: second polynomial vector term
**************************************************/
void poly_qiss_vec_l_sub(poly_qiss_vec_l res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_L_ISS entries by a integer scalar
* 
* Arguments:   - poly_qiss_vec_l res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qiss_vec_l arg: polynomial vector factor
* 			   		 - coeff_qiss fac: integer factor
**************************************************/
void poly_qiss_vec_l_mul_scalar(poly_qiss_vec_l res, const poly_qiss_vec_l arg, const coeff_qiss fac) {
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_mul_scalar(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_L_ISS entries
* 
* Arguments:   - poly_qiss res: polynomial to host the inner product (initialized)
* 			   		 - const poly_qiss_vec_l lhs: first polynomial vector
* 			   		 - const poly_qiss_vec_l rhs: second polynomial vector
**************************************************/
void poly_qiss_vec_l_mul_inner(poly_qiss res, const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs) {
	poly_qiss_zero(res);
	for (size_t i = 0; i < PARAM_L_ISS; ++i) {
		poly_qiss_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_qiss_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_qiss_vec_l_norm2
*
* Description: Compute the square l2 norm of a polynomial vector
* 
* Arguments:   - const poly_qiss_vec_l arg: the polynomial vector
* 
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_qiss_vec_l_norm2(const poly_qiss_vec_l arg) {
  uint64_t sq_norm2, tmp;
  size_t i;
  sq_norm2 = 0;
  for (i = 0; i < PARAM_L_ISS; i++) {
  	tmp = poly_qiss_sq_norm2(arg->entries[i]);
		CHK_UI_OVF_ADDITION(sq_norm2, tmp);
  }
  return sq_norm2;
}

/*************************************************
* Name:        poly_qiss_vec_l_equal
*
* Description: Equality test between two polynomial vectors with PARAM_L_ISS entries
* 
* Arguments:   - const poly_qiss_vec_l lhs: first polynomial vector
* 			   		 - const poly_qiss_vec_l rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_qiss_vec_l_equal(const poly_qiss_vec_l lhs, const poly_qiss_vec_l rhs) {
  for (size_t i = 0; i < PARAM_L_ISS; ++i) {
    if (!poly_qiss_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qiss_vec_l_dump
*
* Description: Print a polynomial vector with PARAM_L_ISS entries
* 
* Arguments:   - const poly_qiss_vec_l arg: polynomial vector to be printed
**************************************************/
void poly_qiss_vec_l_dump(const poly_qiss_vec_l arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_L_ISS - 1; ++i) {
		poly_qiss_dump(arg->entries[i]);
		printf(", ");
	}
	poly_qiss_dump(arg->entries[PARAM_L_ISS - 1]);
	printf("]");
}

/*************************************************
* Name:        poly_qiss_vec_l_pack
*
* Description: Pack a polynomial vector mod PARAM_Q_ISS into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQISS_VECD_PACKEDBYTES bytes)
*              - const poly_qiss arg: the polynomial vector to be packed
**************************************************/
void poly_qiss_vec_l_pack(uint8_t buf[POLYQISS_VECL_PACKEDBYTES], const poly_qiss_vec_l arg) {
	for (size_t i = 0; i < PARAM_L_ISS; i++) {
		poly_qiss_pack(&buf[i*POLYQISS_PACKEDBYTES], arg->entries[i]);
	}
}
