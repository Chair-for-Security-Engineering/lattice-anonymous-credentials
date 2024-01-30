#include "arith_q.h"
#include "macros.h"

static poly_q TMP;

/*************************************************
* Name:        poly_q_mat_d_k_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_K entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_q_mat_d_k_setup(void) {
  poly_q_init(TMP);
}

/*************************************************
* Name:        poly_q_mat_d_k_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_K entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_q_mat_d_k_teardown(void) {
  poly_q_clear(TMP);
}

/*************************************************
* Name:        poly_q_mat_d_k_init
*
* Description: Initialize polynomial matrix with PARAM_D x PARAM_K entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to be initialized
**************************************************/
void poly_q_mat_d_k_init(poly_q_mat_d_k res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_clear
*
* Description: Clear polynomial matrix with PARAM_D x PARAM_K entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to be cleared
**************************************************/
void poly_q_mat_d_k_clear(poly_q_mat_d_k res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_zero
*
* Description: Set an initialized polynomial matrix with PARAM_D x PARAM_K entries to zero
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_q_mat_d_k_zero(poly_q_mat_d_k res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_set
*
* Description: Set a polynomial matrix with PARAM_D x PARAM_K entries equal to another polynomial matrix
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to be set (initialized)
*              - const poly_q_mat_d_k arg: polynomial matrix to be read
**************************************************/
void poly_q_mat_d_k_set(poly_q_mat_d_k res, const poly_q_mat_d_k arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_neg
*
* Description: Negate a polynomial matrix with PARAM_D x PARAM_K entries
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to host the negation (initialized)
*              - const poly_q_mat_d_k arg: polynomial matrix to be negated
**************************************************/
void poly_q_mat_d_k_neg(poly_q_mat_d_k res, const poly_q_mat_d_k arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_add
*
* Description: Add two polynomial matrices with PARAM_D x PARAM_K entries
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to host the sum (initialized)
*              - const poly_q_mat_d_k lhs: first polynomial matrix summand
*              - const poly_q_mat_d_k rhs: second polynomial matrix summand
**************************************************/
void poly_q_mat_d_k_add(poly_q_mat_d_k res, const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_sub
*
* Description: Substract two polynomial matrices with PARAM_D x PARAM_K entries
* 
* Arguments:   - poly_q_mat_d_k res: polynomial matrix to host the difference (initialized)
*              - const poly_q_mat_d_k lhs: first polynomial matrix term
*              - const poly_q_mat_d_k rhs: second polynomial matrix term
**************************************************/
void poly_q_mat_d_k_sub(poly_q_mat_d_k res, const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_mul_vec_k
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_K entries
*              with a polynomial vector with PARAM_K entries
* 
* Arguments:   - poly_q_vec_d res: polynomial vector to host the multiplication (initialized)
*              - const poly_q_mat_d_k lhs: polynomial matrix to multiply
*              - const poly_q_vec_k rhs: polynomial vector to multiply
**************************************************/
void poly_q_mat_d_k_mul_vec_k(poly_q_vec_d res, const poly_q_mat_d_k lhs, const poly_q_vec_k rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_mul_inner(TMP, lhs->rows[i], rhs);
    poly_q_vec_d_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_equal
*
* Description: Equality test between two polynomial matrices with PARAM_D x PARAM_K entries
* 
* Arguments:   - const poly_q_mat_d_k lhs: first polynomial matrix
*              - const poly_q_mat_d_k rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_q_mat_d_k_equal(const poly_q_mat_d_k lhs, const poly_q_mat_d_k rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    if (!poly_q_vec_k_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_q_mat_d_k_dump
*
* Description: Print a polynomial matrix with PARAM_D x PARAM_K entries
* 
* Arguments:   - const poly_q_mat_d_k arg: polynomial matrix to be printed
**************************************************/
void poly_q_mat_d_k_dump(const poly_q_mat_d_k arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D - 1; ++i) {
		poly_q_vec_k_dump(arg->rows[i]);
		printf(", ");
	}
	poly_q_vec_k_dump(arg->rows[PARAM_D - 1]);
	printf("]");
}
