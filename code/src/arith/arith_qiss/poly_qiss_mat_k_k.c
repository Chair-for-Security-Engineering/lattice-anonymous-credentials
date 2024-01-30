#include "arith_qiss.h"
#include "macros.h"

static poly_qiss TMP;

/*************************************************
* Name:        poly_qiss_mat_k_k_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_ISS integer matrices with PARAM_K_ISS x PARAM_K_ISS entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_qiss_mat_k_k_setup(void) {
  poly_qiss_init(TMP);
}

/*************************************************
* Name:        poly_qiss_mat_k_k_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_ISS integer matrices with PARAM_K_ISS x PARAM_K_ISS entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_qiss_mat_k_k_teardown(void) {
  poly_qiss_clear(TMP);
}

/*************************************************
* Name:        poly_qiss_mat_k_k_init
*
* Description: Initialize polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to be initialized
**************************************************/
void poly_qiss_mat_k_k_init(poly_qiss_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_clear
*
* Description: Clear polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to be cleared
**************************************************/
void poly_qiss_mat_k_k_clear(poly_qiss_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_zero
*
* Description: Set an initialized polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries to zero
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_qiss_mat_k_k_zero(poly_qiss_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_set
*
* Description: Set a polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries equal to another polynomial matrix
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to be set (initialized)
*              - const poly_qiss_mat_k_k arg: polynomial matrix to be read
**************************************************/
void poly_qiss_mat_k_k_set(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k arg) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_neg
*
* Description: Negate a polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to host the negation (initialized)
*              - const poly_qiss_mat_k_k arg: polynomial matrix to be negated
**************************************************/
void poly_qiss_mat_k_k_neg(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k arg) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_add
*
* Description: Add two polynomial matrices with PARAM_K_ISS x PARAM_K_ISS entries
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to host the sum (initialized)
*              - const poly_qiss_mat_k_k lhs: first polynomial matrix summand
*              - const poly_qiss_mat_k_k rhs: second polynomial matrix summand
**************************************************/
void poly_qiss_mat_k_k_add(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_sub
*
* Description: Substract two polynomial matrices with PARAM_K_ISS x PARAM_K_ISS entries
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix to host the difference (initialized)
*              - const poly_qiss_mat_k_k lhs: first polynomial matrix term
*              - const poly_qiss_mat_k_k rhs: second polynomial matrix term
**************************************************/
void poly_qiss_mat_k_k_sub(poly_qiss_mat_k_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_mul_vec_k
*
* Description: Product of a polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries
*              with a polynomial vector with PARAM_K_ISS entries
* 
* Arguments:   - poly_qiss_vec_k res: polynomial vector to host the multiplication (initialized)
*              - const poly_qiss_mat_k_k lhs: polynomial matrix to multiply
*              - const poly_qiss_vec_k rhs: polynomial vector to multiply
**************************************************/
void poly_qiss_mat_k_k_mul_vec_k(poly_qiss_vec_k res, const poly_qiss_mat_k_k lhs, const poly_qiss_vec_k rhs) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    poly_qiss_vec_k_mul_inner(TMP, lhs->rows[i], rhs);
    poly_qiss_vec_k_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_qiss_mat_k_k_equal
*
* Description: Equality test between two polynomial matrices with PARAM_K_ISS x PARAM_K_ISS entries
* 
* Arguments:   - const poly_qiss_mat_k_k lhs: first polynomial matrix
*              - const poly_qiss_mat_k_k rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_qiss_mat_k_k_equal(const poly_qiss_mat_k_k lhs, const poly_qiss_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_ISS; ++i) {
    if (!poly_qiss_vec_k_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qiss_mat_k_k_dump
*
* Description: Print a polynomial matrix with PARAM_K_ISS x PARAM_K_ISS entries
* 
* Arguments:   - const poly_qiss_mat_k_k arg: polynomial matrix to be printed
**************************************************/
void poly_qiss_mat_k_k_dump(const poly_qiss_mat_k_k arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_K_ISS - 1; ++i) {
		poly_qiss_vec_k_dump(arg->rows[i]);
		printf(", ");
	}
	poly_qiss_vec_k_dump(arg->rows[PARAM_K_ISS - 1]);
	printf("]");
}
