#include "arith_qshow.h"
#include "macros.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_mat_k_k_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_K_SHOW x PARAM_K_SHOW entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_qshow_mat_k_k_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_k_k_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_K_SHOW x PARAM_K_SHOW entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_qshow_mat_k_k_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_k_k_init
*
* Description: Initialize polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to be initialized
**************************************************/
void poly_qshow_mat_k_k_init(poly_qshow_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_clear
*
* Description: Clear polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to be cleared
**************************************************/
void poly_qshow_mat_k_k_clear(poly_qshow_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_zero
*
* Description: Set an initialized polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries to zero
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_qshow_mat_k_k_zero(poly_qshow_mat_k_k res) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_set
*
* Description: Set a polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries equal to another polynomial matrix
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to be set (initialized)
*              - const poly_qshow_mat_k_k arg: polynomial matrix to be read
**************************************************/
void poly_qshow_mat_k_k_set(poly_qshow_mat_k_k res, const poly_qshow_mat_k_k arg) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_neg
*
* Description: Negate a polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to host the negation (initialized)
*              - const poly_qshow_mat_k_k arg: polynomial matrix to be negated
**************************************************/
void poly_qshow_mat_k_k_neg(poly_qshow_mat_k_k res, const poly_qshow_mat_k_k arg) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_add
*
* Description: Add two polynomial matrices with PARAM_K_SHOW x PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to host the sum (initialized)
*              - const poly_qshow_mat_k_k lhs: first polynomial matrix summand
*              - const poly_qshow_mat_k_k rhs: second polynomial matrix summand
**************************************************/
void poly_qshow_mat_k_k_add(poly_qshow_mat_k_k res, const poly_qshow_mat_k_k lhs, const poly_qshow_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_sub
*
* Description: Substract two polynomial matrices with PARAM_K_SHOW x PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to host the difference (initialized)
*              - const poly_qshow_mat_k_k lhs: first polynomial matrix term
*              - const poly_qshow_mat_k_k rhs: second polynomial matrix term
**************************************************/
void poly_qshow_mat_k_k_sub(poly_qshow_mat_k_k res, const poly_qshow_mat_k_k lhs, const poly_qshow_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_mul_vec_k
*
* Description: Product of a polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries
*              with a polynomial vector with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the multiplication (initialized)
*              - const poly_qshow_mat_k_k lhs: polynomial matrix to multiply
*              - const poly_qshow_vec_k rhs: polynomial vector to multiply
**************************************************/
void poly_qshow_mat_k_k_mul_vec_k(poly_qshow_vec_k res, const poly_qshow_mat_k_k lhs, const poly_qshow_vec_k rhs) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_vec_k_mul_inner(TMP, lhs->rows[i], rhs);
    poly_qshow_vec_k_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_chal_3_embed
*
* Description: Permuted subring multiplication matrix embedding from the 
*              subring embedding vector. Used in computing the quadratic
*              terms sum_i mu_{l + i}G_i" in round 3 and 4 of the proof 
*              of valid credential (anonymous credentials show)
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix to host the embedding
*              - const poly_qshow_vec_k arg: polynomial vector to be embedded
**************************************************/
void poly_qshow_mat_k_k_chal_3_embed(poly_qshow_mat_k_k res, const poly_qshow_vec_k arg) {
  size_t i,j;
  // all except anti-diagonal
  for (i = 0; i < PARAM_K_SHOW - 1; i++) {
    poly_qshow_mul_x(TMP, arg->entries[i]);
    for (j = 0; j < i + 1; j++) {
      poly_qshow_set(res->rows[j]->entries[i-j], arg->entries[i]);
    }
    for (j = i+1; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(res->rows[j]->entries[PARAM_K_SHOW + i - j], TMP);
    }
  }
  // anti-diagonal
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_set(res->rows[i]->entries[PARAM_K_SHOW - 1 - i], arg->entries[PARAM_K_SHOW - 1]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_k_k_equal
*
* Description: Equality test between two polynomial matrices with PARAM_K_SHOW x PARAM_K_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_k_k lhs: first polynomial matrix
*              - const poly_qshow_mat_k_k rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_qshow_mat_k_k_equal(const poly_qshow_mat_k_k lhs, const poly_qshow_mat_k_k rhs) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    if (!poly_qshow_vec_k_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_mat_k_k_dump
*
* Description: Print a polynomial matrix with PARAM_K_SHOW x PARAM_K_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_k_k arg: polynomial matrix to be printed
**************************************************/
void poly_qshow_mat_k_k_dump(const poly_qshow_mat_k_k arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_K_SHOW - 1; ++i) {
		poly_qshow_vec_k_dump(arg->rows[i]);
		printf(", ");
	}
	poly_qshow_vec_k_dump(arg->rows[PARAM_K_SHOW - 1]);
	printf("]");
}
