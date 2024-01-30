#include "arith_qshow.h"
#include "macros.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_mat_d_m1_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_D_SHOW x PARAM_M1_SHOW entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_qshow_mat_d_m1_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_D_SHOW x PARAM_M1_SHOW entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_qshow_mat_d_m1_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_init
*
* Description: Initialize polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to be initialized
**************************************************/
void poly_qshow_mat_d_m1_init(poly_qshow_mat_d_m1 res) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_clear
*
* Description: Clear polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to be cleared
**************************************************/
void poly_qshow_mat_d_m1_clear(poly_qshow_mat_d_m1 res) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_zero
*
* Description: Set an initialized polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries to zero
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_qshow_mat_d_m1_zero(poly_qshow_mat_d_m1 res) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_set
*
* Description: Set a polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries equal to another polynomial matrix
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to be set (initialized)
*              - const poly_qshow_mat_d_m1 arg: polynomial matrix to be read
**************************************************/
void poly_qshow_mat_d_m1_set(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 arg) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_neg
*
* Description: Negate a polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to host the negation (initialized)
*              - const poly_qshow_mat_d_m1 arg: polynomial matrix to be negated
**************************************************/
void poly_qshow_mat_d_m1_neg(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 arg) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_add
*
* Description: Add two polynomial matrices with PARAM_D_SHOW x PARAM_M1_SHOW entries
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to host the sum (initialized)
*              - const poly_qshow_mat_d_m1 lhs: first polynomial matrix summand
*              - const poly_qshow_mat_d_m1 rhs: second polynomial matrix summand
**************************************************/
void poly_qshow_mat_d_m1_add(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_sub
*
* Description: Substract two polynomial matrices with PARAM_D_SHOW x PARAM_M1_SHOW entries
* 
* Arguments:   - poly_qshow_mat_d_m1 res: polynomial matrix to host the difference (initialized)
*              - const poly_qshow_mat_d_m1 lhs: first polynomial matrix term
*              - const poly_qshow_mat_d_m1 rhs: second polynomial matrix term
**************************************************/
void poly_qshow_mat_d_m1_sub(poly_qshow_mat_d_m1 res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_mul_vec_m1
*
* Description: Product of a polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries
*              with a polynomial vector with PARAM_M1_SHOW entries
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the multiplication (initialized)
*              - const poly_qshow_mat_d_m1 lhs: polynomial matrix to multiply
*              - const poly_qshow_vec_m1 rhs: polynomial vector to multiply
**************************************************/
void poly_qshow_mat_d_m1_mul_vec_m1(poly_qshow_vec_d res, const poly_qshow_mat_d_m1 lhs, const poly_qshow_vec_m1 rhs) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_vec_m1_mul_inner(TMP, lhs->rows[i], rhs);
    poly_qshow_vec_d_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_equal
*
* Description: Equality test between two polynomial matrices with PARAM_D_SHOW x PARAM_M1_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_d_m1 lhs: first polynomial matrix
*              - const poly_qshow_mat_d_m1 rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_qshow_mat_d_m1_equal(const poly_qshow_mat_d_m1 lhs, const poly_qshow_mat_d_m1 rhs) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    if (!poly_qshow_vec_m1_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_dump
*
* Description: Print a polynomial matrix with PARAM_D_SHOW x PARAM_M1_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_d_m1 arg: polynomial matrix to be printed
**************************************************/
void poly_qshow_mat_d_m1_dump(const poly_qshow_mat_d_m1 arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D_SHOW - 1; ++i) {
		poly_qshow_vec_m1_dump(arg->rows[i]);
		printf(", ");
	}
	poly_qshow_vec_m1_dump(arg->rows[PARAM_D_SHOW - 1]);
	printf("]");
}
