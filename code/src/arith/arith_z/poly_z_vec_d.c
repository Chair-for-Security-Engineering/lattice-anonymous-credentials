#include "arith_z.h"
#include "macros.h"

static poly_z TMP;

/*************************************************
* Name:        poly_z_vec_d_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              multiprecision integer vectors with PARAM_D entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_z_vec_d_setup(void) {
  poly_z_init(TMP);
}

/*************************************************
* Name:        poly_z_vec_d_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              multiprecision integer vectors with PARAM_D entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_z_vec_d_teardown(void) {
  poly_z_clear(TMP);
}

/*************************************************
* Name:        poly_z_vec_d_init
*
* Description: Initialize polynomial vector with PARAM_D entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_z_vec_d arg: polynomial vector to be initialized
**************************************************/
void poly_z_vec_d_init(poly_z_vec_d arg) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_clear
*
* Description: Clear polynomial vector with PARAM_D entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_z_vec_d arg: polynomial vector to be cleared
**************************************************/
void poly_z_vec_d_clear(poly_z_vec_d arg) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_zero
*
* Description: Set an initialized polynomial vector with PARAM_D entries to zero
* 
* Arguments:   - poly_z_vec_d arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_z_vec_d_zero(poly_z_vec_d arg) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_set
*
* Description: Set a polynomial vector with PARAM_D entries equal to another polynomial vector
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_z_vec_d arg: polynomial vector to be read
**************************************************/
void poly_z_vec_d_set(poly_z_vec_d res, const poly_z_vec_d arg) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_D]
* 
* Arguments:   - poly_z res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_z_vec_d arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_z_vec_d_get_poly(poly_z res, const poly_z_vec_d arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_D, "Illegal argument: cannot get entry of vector at given position.");
	poly_z_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_z_vec_d_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_D]
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_z arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_z_vec_d_set_poly(poly_z_vec_d res, const poly_z arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_D, "Illegal argument: cannot set entry of vector at given position.");
	poly_z_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_z_vec_d_neg
*
* Description: Negate a polynomial vector with PARAM_D entries
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_z_vec_d arg: polynomial vector to be negated
**************************************************/
void poly_z_vec_d_neg(poly_z_vec_d res, const poly_z_vec_d arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_z_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_z_vec_d_add
*
* Description: Add two polynomial vectors with PARAM_D entries
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_z_vec_d lhs: first polynomial vector summand
* 			   		 - const poly_z_vec_d rhs: second polynomial vector summand
**************************************************/
void poly_z_vec_d_add(poly_z_vec_d res, const poly_z_vec_d lhs, const poly_z_vec_d rhs) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_sub
*
* Description: Substract two polynomial vectors with PARAM_D entries
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_z_vec_d lhs: first polynomial vector term
* 			   		 - const poly_z_vec_d rhs: second polynomial vector term
**************************************************/
void poly_z_vec_d_sub(poly_z_vec_d res, const poly_z_vec_d lhs, const poly_z_vec_d rhs) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_z_vec_d_mul_poly
*
* Description: Multiplication of a polynomial vector with PARAM_D entries by a polynomial
* 
* Arguments:   - poly_z_vec_d res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_z_vec_d arg: first polynomial vector factor
* 			   		 - const poly_z fac: second polynomial factor
**************************************************/
void poly_z_vec_d_mul_poly(poly_z_vec_d res, const poly_z_vec_d arg, const poly_z fac) {
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_mul(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_z_vec_d_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_D entries
* 
* Arguments:   - poly_z res: polynomial to host the inner product (initialized)
* 			   		 - const poly_z_vec_d lhs: first polynomial vector
* 			   		 - const poly_z_vec_d rhs: second polynomial vector
**************************************************/
void poly_z_vec_d_mul_inner(poly_z res, const poly_z_vec_d lhs, const poly_z_vec_d rhs) {
	poly_z_zero(res);
	for (size_t i = 0; i < PARAM_D; ++i) {
		poly_z_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_z_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_z_vec_d_equal
*
* Description: Equality test between two polynomial vectors with PARAM_D entries
* 
* Arguments:   - const poly_z_vec_d lhs: first polynomial vector
* 			   		 - const poly_z_vec_d rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_z_vec_d_equal(const poly_z_vec_d lhs, const poly_z_vec_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    if (!poly_z_equal(lhs->entries[i], rhs->entries[i])) {
        return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_z_vec_d_dump
*
* Description: Print a polynomial vector with PARAM_D entries
* 
* Arguments:   - const poly_z_vec_d arg: polynomial vector to be printed
**************************************************/
void poly_z_vec_d_dump(const poly_z_vec_d arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D - 1; ++i) {
		poly_z_dump(arg->entries[i]);
		printf(", ");
	}
	poly_z_dump(arg->entries[PARAM_D - 1]);
	printf("]");
}
