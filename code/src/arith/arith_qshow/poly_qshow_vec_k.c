#include "arith_qshow.h"
#include "macros.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_vec_k_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_K_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_qshow_vec_k_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_k_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_K_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_qshow_vec_k_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_k_init
*
* Description: Initialize polynomial vector with PARAM_K_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_qshow_vec_k arg: polynomial vector to be initialized
**************************************************/
void poly_qshow_vec_k_init(poly_qshow_vec_k arg) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_clear
*
* Description: Clear polynomial vector with PARAM_K_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_vec_k arg: polynomial vector to be cleared
**************************************************/
void poly_qshow_vec_k_clear(poly_qshow_vec_k arg) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_zero
*
* Description: Set an initialized polynomial vector with PARAM_K_SHOW entries to zero
* 
* Arguments:   - poly_qshow_vec_k arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_qshow_vec_k_zero(poly_qshow_vec_k arg) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_set
*
* Description: Set a polynomial vector with PARAM_K_SHOW entries equal to another polynomial vector
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow_vec_k arg: polynomial vector to be read
**************************************************/
void poly_qshow_vec_k_set(poly_qshow_vec_k res, const poly_qshow_vec_k arg) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_K_SHOW]
* 
* Arguments:   - poly_qshow res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_qshow_vec_k arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_k_get_poly(poly_qshow res, const poly_qshow_vec_k arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_K_SHOW, "Illegal argument: cannot get entry of vector at given position.");
	poly_qshow_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_qshow_vec_k_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_K_SHOW]
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_k_set_poly(poly_qshow_vec_k res, const poly_qshow arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_K_SHOW, "Illegal argument: cannot set entry of vector at given position.");
	poly_qshow_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_qshow_vec_k_neg
*
* Description: Negate a polynomial vector with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_qshow_vec_k arg: polynomial vector to be negated
**************************************************/
void poly_qshow_vec_k_neg(poly_qshow_vec_k res, const poly_qshow_vec_k arg) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    poly_qshow_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_k_add
*
* Description: Add two polynomial vectors with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_qshow_vec_k lhs: first polynomial vector summand
* 			   		 - const poly_qshow_vec_k rhs: second polynomial vector summand
**************************************************/
void poly_qshow_vec_k_add(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_sub
*
* Description: Substract two polynomial vectors with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_qshow_vec_k lhs: first polynomial vector term
* 			   		 - const poly_qshow_vec_k rhs: second polynomial vector term
**************************************************/
void poly_qshow_vec_k_sub(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_K_SHOW entries by a integer scalar
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qshow_vec_k arg: polynomial vector factor
* 			   		 - coeff_qshow fac: integer factor
**************************************************/
void poly_qshow_vec_k_mul_scalar(poly_qshow_vec_k res, const poly_qshow_vec_k arg, const coeff_qshow fac) {
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_mul_scalar(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_mul_poly_qshow
*
* Description: Multiplication of a polynomial vector with PARAM_K_SHOW entries by a polynomial
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qshow_vec_k lhs: first polynomial vector factor
* 			   		 - const poly_qshow rhs: second polynomial factor
**************************************************/
void poly_qshow_vec_k_mul_poly_qshow(poly_qshow_vec_k res, const poly_qshow_vec_k lhs, const poly_qshow rhs) {
	for (size_t i = 0; i < PARAM_K_SHOW; i++) {
		poly_qshow_mul(res->entries[i], lhs->entries[i], rhs);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow res: polynomial to host the inner product (initialized)
* 			   		 - const poly_qshow_vec_k lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_k rhs: second polynomial vector
**************************************************/
void poly_qshow_vec_k_mul_inner(poly_qshow res, const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs) {
	poly_qshow_zero(res);
	for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
		poly_qshow_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_qshow_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_conjugate
*
* Description: Conjugate a polynomial vector with PARAM_K_SHOW entries
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector to host the conjugate (initialized)
* 			   		 - const poly_qshow_vec_k arg: polynomial vector to be conjugated
**************************************************/
void poly_qshow_vec_k_conjugate(poly_qshow_vec_k res, const poly_qshow_vec_k arg) {
	for (size_t i = 0; i < PARAM_K_SHOW; i++) {
		poly_qshow_conjugate(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_k_norm2
*
* Description: Compute the square l2 norm of a polynomial vector
* 
* Arguments:   - const poly_qshow_vec_k arg: the polynomial vector
* 
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_qshow_vec_k_norm2(const poly_qshow_vec_k arg) {
  uint64_t sq_norm2, tmp;
  size_t i;
  sq_norm2 = 0;
  for (i = 0; i < PARAM_K_SHOW; i++) {
  	tmp = poly_qshow_sq_norm2(arg->entries[i]);
		CHK_UI_OVF_ADDITION(sq_norm2, tmp);
  }
  return sq_norm2;
}

/*************************************************
* Name:        poly_qshow_vec_k_equal
*
* Description: Equality test between two polynomial vectors with PARAM_K_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_k lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_k rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_qshow_vec_k_equal(const poly_qshow_vec_k lhs, const poly_qshow_vec_k rhs) {
  for (size_t i = 0; i < PARAM_K_SHOW; ++i) {
    if (!poly_qshow_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_vec_k_dump
*
* Description: Print a polynomial vector with PARAM_K_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_k arg: polynomial vector to be printed
**************************************************/
void poly_qshow_vec_k_dump(const poly_qshow_vec_k arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_K_SHOW - 1; ++i) {
		poly_qshow_dump(arg->entries[i]);
		printf(", ");
	}
	poly_qshow_dump(arg->entries[PARAM_K_SHOW - 1]);
	printf("]");
}

/*************************************************
* Name:        poly_qshow_vec_k_pack
*
* Description: Pack a polynomial vector mod PARAM_Q_SHOW into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQSHOW_VECK_PACKEDBYTES bytes)
*              - const poly_qshow arg: the polynomial vector to be packed
**************************************************/
void poly_qshow_vec_k_pack(uint8_t buf[POLYQSHOW_VECK_PACKEDBYTES], const poly_qshow_vec_k arg) {
  for (size_t i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_pack(&buf[i * POLYQSHOW_PACKEDBYTES], arg->entries[i]);
  }
}
