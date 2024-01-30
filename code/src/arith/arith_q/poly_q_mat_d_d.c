#include "arith_q.h"
#include "macros.h"

static poly_q TMP;
static poly_q ACC;

/*************************************************
* Name:        poly_q_mat_d_d_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_D entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_q_mat_d_d_setup(void) {
  poly_q_init(TMP);
  poly_q_init(ACC);
}

/*************************************************
* Name:        poly_q_mat_d_d_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_D entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_q_mat_d_d_teardown(void) {
  poly_q_clear(TMP);
  poly_q_clear(ACC);
}

/*************************************************
* Name:        poly_q_mat_d_d_init
*
* Description: Initialize polynomial matrix with PARAM_D x PARAM_D entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to be initialized
**************************************************/
void poly_q_mat_d_d_init(poly_q_mat_d_d res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_clear
*
* Description: Clear polynomial matrix with PARAM_D x PARAM_D entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to be cleared
**************************************************/
void poly_q_mat_d_d_clear(poly_q_mat_d_d res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_zero
*
* Description: Set an initialized polynomial matrix with PARAM_D x PARAM_D entries to zero
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_q_mat_d_d_zero(poly_q_mat_d_d res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_set
*
* Description: Set a polynomial matrix with PARAM_D x PARAM_D entries equal to another polynomial matrix
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to be set (initialized)
*              - const poly_q_mat_d_d arg: polynomial matrix to be read
**************************************************/
void poly_q_mat_d_d_set(poly_q_mat_d_d res, const poly_q_mat_d_d arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_neg
*
* Description: Negate a polynomial matrix with PARAM_D x PARAM_D entries
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to host the negation (initialized)
*              - const poly_q_mat_d_d arg: polynomial matrix to be negated
**************************************************/
void poly_q_mat_d_d_neg(poly_q_mat_d_d res, const poly_q_mat_d_d arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_add
*
* Description: Add two polynomial matrices with PARAM_D x PARAM_D entries
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to host the sum (initialized)
*              - const poly_q_mat_d_d lhs: first polynomial matrix summand
*              - const poly_q_mat_d_d rhs: second polynomial matrix summand
**************************************************/
void poly_q_mat_d_d_add(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_sub
*
* Description: Substract two polynomial matrices with PARAM_D x PARAM_D entries
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to host the difference (initialized)
*              - const poly_q_mat_d_d lhs: first polynomial matrix term
*              - const poly_q_mat_d_d rhs: second polynomial matrix term
**************************************************/
void poly_q_mat_d_d_sub(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_mul_vec_d
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_D entries
*              with a polynomial vector with PARAM_D entries
* 
* Arguments:   - poly_q_vec_d res: polynomial vector to host the multiplication (initialized)
*              - const poly_q_mat_d_d lhs: polynomial matrix to multiply
*              - const poly_q_vec_d rhs: polynomial vector to multiply
**************************************************/
void poly_q_mat_d_d_mul_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_mul_inner(TMP, lhs->rows[i], rhs);
    poly_q_vec_d_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_muladd_vec_d
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_D entries
*              with a polynomial vector with PARAM_D entries incremented to res
* 
* Arguments:   - poly_q_vec_d res: polynomial vector to increment with the multiplication (initialized)
*              - const poly_q_mat_d_d lhs: polynomial matrix to multiply
*              - const poly_q_vec_d rhs: polynomial vector to multiply
**************************************************/
void poly_q_mat_d_d_muladd_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_mul_inner(TMP, lhs->rows[i], rhs);
    poly_q_add(res->entries[i], res->entries[i], TMP);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_muladd_vec_d
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_D entries
*              with a polynomial vector with PARAM_D entries decremented to res
* 
* Arguments:   - poly_q_vec_d res: polynomial vector to decrement with the multiplication (initialized)
*              - const poly_q_mat_d_d lhs: polynomial matrix to multiply
*              - const poly_q_vec_d rhs: polynomial vector to multiply
**************************************************/
void poly_q_mat_d_d_mulsub_vec_d(poly_q_vec_d res, const poly_q_mat_d_d lhs, const poly_q_vec_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_d_mul_inner(TMP, lhs->rows[i], rhs);
    poly_q_sub(res->entries[i], res->entries[i], TMP);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_mul_mat_d_d
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_D entries
*              with another polynomial matrix with PARAM_D x PARAM_D entries
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to host the multiplication (initialized)
*              - const poly_q_mat_d_d lhs: first polynomial matrix factor
*              - const poly_q_mat_d_d lhs: second polynomial matrix factor
**************************************************/
void poly_q_mat_d_d_mul_mat_d_d(poly_q_mat_d_d res, const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs) {
  for (size_t r = 0; r < PARAM_D; ++r) {
    for (size_t c = 0; c < PARAM_D; ++c) {
      poly_q_zero(ACC);
		for (size_t i = 0; i < PARAM_D; ++i) {
      poly_q_mul(TMP, lhs->rows[r]->entries[i], rhs->rows[i]->entries[c]);
      poly_q_add(ACC, ACC, TMP);
		}
		poly_q_set(res->rows[r]->entries[c], ACC);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_conjugate
*
* Description: Conjugate a polynomial matrix with PARAM_D x PARAM_D entries
* 
* Arguments:   - poly_q_mat_d_d res: polynomial matrix to host the conjugate (initialized)
*              - const poly_q_mat_d_d arg: polynomial matrix to be conjugated
**************************************************/
void poly_q_mat_d_d_conjugate(poly_q_mat_d_d res, const poly_q_mat_d_d arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    for (size_t j = 0; j < PARAM_D; ++j) {
      poly_q_conjugate(TMP, arg->rows[i]->entries[j]);
		  poly_q_set(res->rows[j]->entries[i], TMP);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_equal
*
* Description: Equality test between two polynomial matrices with PARAM_D x PARAM_D entries
* 
* Arguments:   - const poly_q_mat_d_d lhs: first polynomial matrix
*              - const poly_q_mat_d_d rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_q_mat_d_d_equal(const poly_q_mat_d_d lhs, const poly_q_mat_d_d rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    if (!poly_q_vec_d_equal(lhs->rows[i], rhs->rows[i])) {
     return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_q_mat_d_d_dump
*
* Description: Print a polynomial matrix with PARAM_D x PARAM_D entries
* 
* Arguments:   - const poly_q_mat_d_d arg: polynomial matrix to be printed
**************************************************/
void poly_q_mat_d_d_dump(const poly_q_mat_d_d arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D - 1; ++i) {
		poly_q_vec_d_dump(arg->rows[i]);
		printf(", ");
	}
	poly_q_vec_d_dump(arg->rows[PARAM_D - 1]);
	printf("]");
}


