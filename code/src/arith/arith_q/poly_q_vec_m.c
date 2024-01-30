#include "arith_q.h"
#include "macros.h"

static poly_q TMP;

/*************************************************
* Name:        poly_q_vec_m_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_M entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_q_vec_m_setup(void) {
  poly_q_init(TMP);
}

/*************************************************
* Name:        poly_q_vec_m_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_M entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_q_vec_m_teardown(void) {
  poly_q_clear(TMP);
}

/*************************************************
* Name:        poly_q_vec_m_init
*
* Description: Initialize polynomial vector with PARAM_M entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_q_vec_m arg: polynomial vector to be initialized
**************************************************/
void poly_q_vec_m_init(poly_q_vec_m arg) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_clear
*
* Description: Clear polynomial vector with PARAM_M entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_q_vec_m arg: polynomial vector to be cleared
**************************************************/
void poly_q_vec_m_clear(poly_q_vec_m arg) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_zero
*
* Description: Set an initialized polynomial vector with PARAM_M entries to zero
* 
* Arguments:   - poly_q_vec_m arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_q_vec_m_zero(poly_q_vec_m arg) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_set
*
* Description: Set a polynomial vector with PARAM_M entries equal to another polynomial vector
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to be set (initialized)
* 			   		 - const poly_q_vec_m arg: polynomial vector to be read
**************************************************/
void poly_q_vec_m_set(poly_q_vec_m res, const poly_q_vec_m arg) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_M]
* 
* Arguments:   - poly_q res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_q_vec_m arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_q_vec_m_get_poly(poly_q res, const poly_q_vec_m arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_M, "Illegal argument: cannot get entry of vector at given position.");
	poly_q_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_q_vec_m_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_M]
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to be set (initialized)
* 			   		 - const poly_q arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_q_vec_m_set_poly(poly_q_vec_m res, const poly_q arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_M, "Illegal argument: cannot set entry of vector at given position.");
	poly_q_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_q_vec_m_neg
*
* Description: Negate a polynomial vector with PARAM_M entries
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_q_vec_m arg: polynomial vector to be negated
**************************************************/
void poly_q_vec_m_neg(poly_q_vec_m res, const poly_q_vec_m arg) {
  for (size_t i = 0; i < PARAM_M; ++i) {
    poly_q_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_q_vec_m_add
*
* Description: Add two polynomial vectors with PARAM_M entries
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_q_vec_m lhs: first polynomial vector summand
* 			   		 - const poly_q_vec_m rhs: second polynomial vector summand
**************************************************/
void poly_q_vec_m_add(poly_q_vec_m res, const poly_q_vec_m lhs, const poly_q_vec_m rhs) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_sub
*
* Description: Substract two polynomial vectors with PARAM_M entries
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_q_vec_m lhs: first polynomial vector term
* 			   		 - const poly_q_vec_m rhs: second polynomial vector term
**************************************************/
void poly_q_vec_m_sub(poly_q_vec_m res, const poly_q_vec_m lhs, const poly_q_vec_m rhs) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_q_vec_m_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_M entries by a integer scalar
* 
* Arguments:   - poly_q_vec_m res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_q_vec_m arg: polynomial vector factor
* 			   		 - coeff_q fac: integer factor
**************************************************/
void poly_q_vec_m_mul_scalar(poly_q_vec_m res, const poly_q_vec_m arg, const poly_q fac) {
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_mul(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_q_vec_m_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_M entries
* 
* Arguments:   - poly_q res: polynomial to host the inner product (initialized)
* 			   		 - const poly_q_vec_m lhs: first polynomial vector
* 			   		 - const poly_q_vec_m rhs: second polynomial vector
**************************************************/
void poly_q_vec_m_mul_inner(poly_q res, const poly_q_vec_m lhs, const poly_q_vec_m rhs) {
	poly_q_zero(res);
	for (size_t i = 0; i < PARAM_M; ++i) {
		poly_q_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_q_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_q_vec_m_equal
*
* Description: Equality test between two polynomial vectors with PARAM_M entries
* 
* Arguments:   - const poly_q_vec_m lhs: first polynomial vector
* 			   		 - const poly_q_vec_m rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_q_vec_m_equal(const poly_q_vec_m lhs, const poly_q_vec_m rhs) {
  for (size_t i = 0; i < PARAM_M; ++i) {
    if (!poly_q_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_q_vec_m_dump
*
* Description: Print a polynomial vector with PARAM_M entries
* 
* Arguments:   - const poly_q_vec_m arg: polynomial vector to be printed
**************************************************/
void poly_q_vec_m_dump(const poly_q_vec_m arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_M - 1; ++i) {
		poly_q_dump(arg->entries[i]);
		printf(", ");
	}
	poly_q_dump(arg->entries[PARAM_M - 1]);
	printf("]");
}
