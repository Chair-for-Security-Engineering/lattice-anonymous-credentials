#ifndef ARITH_H
#define ARITH_H

/**********************************************
* This header serves as an umbrella header for all arithmetic backends.
**********************************************/
#include "arith_z.h"
#include "arith_q.h"
#include "arith_qiss.h"
#include "arith_qshow.h"
#include "arith_real.h"
#include "arith_z.h"

void arith_setup(void);
void arith_teardown(void);

/**********************************************
* All functions (such as conversions) which rely on multiple "arithmetics" go here.
* Function that work on different types (i.e. poly vs d-vector vs k-vector vs matrix) should
* be placed in the corresponding "arithmetics".
**********************************************/
void poly_z_mat_d_d_from_poly_q_mat_d_d(poly_z_mat_d_d res, const poly_q_mat_d_d arg);

void poly_real_from_poly_q(poly_real res, const poly_q arg);

void poly_q_from_poly_real(poly_q res, const poly_real arg);

void poly_q_samplefz(poly_q res, const poly_real f, const poly_real c);

void poly_real_sub_poly_real_poly_q(poly_real res, const poly_q lhs, const poly_real rhs);

void poly_qiss_subring_embed_vec_k(poly_qiss_vec_k res, const poly_q arg, const int64_t fac);
void poly_qiss_subring_embed_mat_k_k(poly_qiss_mat_k_k res, const poly_q arg, const int64_t fac);
uint64_t challenge_size_iss(const poly_qiss c);

void poly_qshow_subring_embed_vec_k(poly_qshow_vec_k res, const poly_q arg, const int64_t fac);
void poly_qshow_subring_embed_mat_k_k(poly_qshow_mat_k_k res, const poly_q arg, const int64_t fac);
uint64_t challenge_size_show(const poly_qshow c);

#endif /* ARITH_H */
