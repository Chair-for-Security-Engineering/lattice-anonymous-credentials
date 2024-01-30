#include "arith_qiss.h"
#include "macros.h"

static poly_qiss TMP;

/*************************************************
* Name:        poly_qiss_mat_256l_m2_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_ISS integer matrices with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_qiss_mat_256l_m2_setup(void) {
  poly_qiss_init(TMP);
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_ISS integer matrices with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_qiss_mat_256l_m2_teardown(void) {
  poly_qiss_clear(TMP);
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_init
*
* Description: Initialize polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to be initialized
**************************************************/
void poly_qiss_mat_256l_m2_init(poly_qiss_mat_256l_m2 res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_clear
*
* Description: Clear polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to be cleared
**************************************************/
void poly_qiss_mat_256l_m2_clear(poly_qiss_mat_256l_m2 res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_zero
*
* Description: Set an initialized polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries to zero
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_qiss_mat_256l_m2_zero(poly_qiss_mat_256l_m2 res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_set
*
* Description: Set a polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries equal to another polynomial matrix
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to be set (initialized)
*              - const poly_qiss_mat_256l_m2 arg: polynomial matrix to be read
**************************************************/
void poly_qiss_mat_256l_m2_set(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 arg) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_neg
*
* Description: Negate a polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to host the negation (initialized)
*              - const poly_qiss_mat_256l_m2 arg: polynomial matrix to be negated
**************************************************/
void poly_qiss_mat_256l_m2_neg(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 arg) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_add
*
* Description: Add two polynomial matrices with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to host the sum (initialized)
*              - const poly_qiss_mat_256l_m2 lhs: first polynomial matrix summand
*              - const poly_qiss_mat_256l_m2 rhs: second polynomial matrix summand
**************************************************/
void poly_qiss_mat_256l_m2_add(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_sub
*
* Description: Substract two polynomial matrices with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
* 
* Arguments:   - poly_qiss_mat_256l_m2 res: polynomial matrix to host the difference (initialized)
*              - const poly_qiss_mat_256l_m2 lhs: first polynomial matrix term
*              - const poly_qiss_mat_256l_m2 rhs: second polynomial matrix term
**************************************************/
void poly_qiss_mat_256l_m2_sub(poly_qiss_mat_256l_m2 res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_mul_vec_m2
*
* Description: Product of a polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
*              with a polynomial vector with PARAM_M2_ISS entries
* 
* Arguments:   - poly_qiss_vec_256_l res: polynomial vector to host the multiplication (initialized)
*              - const poly_qiss_mat_256l_m2 lhs: polynomial matrix to multiply
*              - const poly_qiss_vec_m2 rhs: polynomial vector to multiply
**************************************************/
void poly_qiss_mat_256l_m2_mul_vec_m2(poly_qiss_vec_256_l res, const poly_qiss_mat_256l_m2 lhs, const poly_qiss_vec_m2 rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    poly_qiss_vec_m2_mul_inner(TMP, lhs->rows[i], rhs);
    poly_qiss_vec_256_l_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_equal
*
* Description: Equality test between two polynomial matrices with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
* 
* Arguments:   - const poly_qiss_mat_256l_m2 lhs: first polynomial matrix
*              - const poly_qiss_mat_256l_m2 rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_qiss_mat_256l_m2_equal(const poly_qiss_mat_256l_m2 lhs, const poly_qiss_mat_256l_m2 rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS; ++i) {
    if (!poly_qiss_vec_m2_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_dump
*
* Description: Print a polynomial matrix with PARAM_ARP_DIV_N_L_ISS x PARAM_M2_ISS entries
* 
* Arguments:   - const poly_qiss_mat_256l_m2 arg: polynomial matrix to be printed
**************************************************/
void poly_qiss_mat_256l_m2_dump(const poly_qiss_mat_256l_m2 arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_ARP_DIV_N_L_ISS - 1; ++i) {
		poly_qiss_vec_m2_dump(arg->rows[i]);
		printf(", ");
	}
	poly_qiss_vec_m2_dump(arg->rows[PARAM_ARP_DIV_N_L_ISS - 1]);
	printf("]");
}
