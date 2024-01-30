#include "arith.h"
#include "covariance.h"

/*************************************************
* Name:        arith_setup
*
* Description: Initialize and setup the entire arithmetic backend. 
*              This is strictly required and must be called once
*              before any other function from here is used.
**************************************************/
void arith_setup(void) {
	arith_q_setup();
  arith_qiss_setup();
  arith_qshow_setup();
  arith_real_setup();
  arith_z_setup();
  fft_precomp_setup();

  poly_q_vec_d_setup();
  poly_q_vec_m_setup();
  poly_q_vec_k_setup();
  poly_q_mat_d_m_setup();
  poly_q_mat_d_k_setup();
  poly_q_mat_d_d_setup();

  poly_z_vec_d_setup();
  poly_z_mat_d_d_setup();

  poly_qiss_vec_d_setup();
  poly_qiss_vec_k_setup();
  poly_qiss_vec_m2_setup();
  poly_qiss_vec_l_setup();
  poly_qiss_vec_256_setup();
  poly_qiss_vec_256_l_setup();
  poly_qiss_mat_k_k_setup();
  poly_qiss_mat_d_k_setup();
  poly_qiss_mat_256l_m2_setup();
  poly_qiss_mat_d_m2_setup();

  poly_qshow_vec_d_setup();
  poly_qshow_vec_m1_setup();
  poly_qshow_vec_m2_setup();
  poly_qshow_vec_l_setup();
  poly_qshow_vec_k_setup();
  poly_qshow_vec_256_setup();
  poly_qshow_vec_256_l_setup();
  poly_qshow_mat_k_k_setup();
  poly_qshow_mat_d_k_setup();
  poly_qshow_mat_d_m2_setup();
  poly_qshow_mat_256l_m2_setup();
  poly_qshow_mat_d_m1_setup();
}

/*************************************************
* Name:        arith_teardown
*
* Description: Clean up and teardown the entire arithmetic backend. 
*              This is strictly required and must be called once at 
*              the very end to release any resources.
**************************************************/
void arith_teardown(void) {
	arith_q_teardown();
  arith_qiss_teardown();
  arith_qshow_teardown();
  arith_real_teardown();
  arith_z_teardown();
  fft_precomp_teardown();

  poly_q_vec_d_teardown();
  poly_q_vec_m_teardown();
  poly_q_vec_k_teardown();
  poly_q_mat_d_m_teardown();
  poly_q_mat_d_k_teardown();
  poly_q_mat_d_d_teardown();

  poly_z_vec_d_teardown();
  poly_z_mat_d_d_teardown();

  poly_qiss_vec_d_teardown();
  poly_qiss_vec_k_teardown();
  poly_qiss_vec_m2_teardown();
  poly_qiss_vec_l_teardown();
  poly_qiss_vec_256_teardown();
  poly_qiss_vec_256_l_teardown();
  poly_qiss_mat_k_k_teardown();
  poly_qiss_mat_d_k_teardown();
  poly_qiss_mat_256l_m2_teardown();
  poly_qiss_mat_d_m2_teardown();

  poly_qshow_vec_d_teardown();
  poly_qshow_vec_m1_teardown();
  poly_qshow_vec_m2_teardown();
  poly_qshow_vec_l_teardown();
  poly_qshow_vec_k_teardown();
  poly_qshow_vec_256_teardown();
  poly_qshow_vec_256_l_teardown();
  poly_qshow_mat_k_k_teardown();
  poly_qshow_mat_d_k_teardown();
  poly_qshow_mat_d_m2_teardown();
  poly_qshow_mat_256l_m2_teardown();
  poly_qshow_mat_d_m1_teardown();
}

/*************************************************
* Name:        poly_z_mat_d_d_from_poly_q_mat_d_d
*
* Description: Convert a poly_q_mat_d_d into poly_z_mat_d_d
*              for arithmetic over Z[x]/(x^n + 1) instead
*              of Z[x]/(q, x^n + 1)
* 
* Arguments:   - poly_z_mat_d_d res: matrix of poly_z to host the conversion (initialized)
*              - const poly_q_mat_d_d arg: matrix to be converted
**************************************************/
void poly_z_mat_d_d_from_poly_q_mat_d_d(poly_z_mat_d_d res, const poly_q_mat_d_d arg) {
  size_t i,j,k;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      for (k = 0; k < PARAM_N; k++) {
        poly_z_set_coeff_si(
          res->rows[i]->entries[j],
          k,
          poly_q_get_coeff_centered(arg->rows[i]->entries[j], k)
        );
      }
    }
  }
}

/*************************************************
* Name:        poly_real_from_poly_q
*
* Description: Convert a poly_q into poly_real
*              for arithmetic over R[x]/(x^n + 1) instead
*              of Z[x]/(q, x^n + 1)
* 
* Arguments:   - poly_real res: real-valued polynomial to host the conversion (initialized)
*              - const poly_q arg: polynomial to be converted
**************************************************/
void poly_real_from_poly_q(poly_real res, const poly_q arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_real_set_si(
      res,
      i,
      poly_q_get_coeff_centered(arg, i)
    );
  }
}

/*************************************************
* Name:        poly_q_from_poly_real
*
* Description: Convert a poly_real into poly_q by rounding
*              for arithmetic over Z[x]/(q, x^n + 1) instead
*              of R[x]/(x^n + 1)
* 
* Arguments:   - poly_q res: polynomial to host the conversion (initialized)
*              - const poly_real arg: real-valued polynomial to be rounded/converted
**************************************************/
void poly_q_from_poly_real(poly_q res, const poly_real arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(
      res,
      i,
      poly_real_get_coeff_rounded(arg, i)
    );
  }
}

/*************************************************
* Name:        poly_q_samplefz
*
* Description: Gaussian sampler SampleFz over Z[x]/(x^n + 1) from
*              [GM18] and conversion to poly_q
* 
* Arguments:   - poly_q res: polynomial to host the Gaussian sample (initialized)
*              - const poly_real f: variance polynomial for the Gaussian (Careful: not the stddev)
*              - const poly_real c: center polynomial for the Gaussian
**************************************************/
void poly_q_samplefz(poly_q res, const poly_real f, const poly_real c) {
  poly_real tmp;
  poly_real_init(tmp);
  poly_real_samplefz(tmp, f, c);
  poly_q_from_poly_real(res, tmp);
  poly_real_clear(tmp);
}

/*************************************************
* Name:        poly_real_sub_poly_real_poly_q
*
* Description: Substract a poly_real from a poly_q and store in poly_real
* 
* Arguments:   - poly_real res: real-valued polynomial to host the difference (initialized)
*              - const poly_q lhs: polynomial to substract from
*              - const poly_real rhs: real-valued polynomial to be substracted
**************************************************/
void poly_real_sub_poly_real_poly_q(poly_real res, const poly_q lhs, const poly_real rhs) {
  poly_real tmp;
  poly_real_init(tmp);
  poly_real_from_poly_q(tmp, lhs);
  poly_real_sub(res, tmp, rhs);
  poly_real_clear(tmp);
}

/*************************************************
* Name:        poly_qiss_subring_embed_vec_k
*
* Description: Subring embedding and lifting mod PARAM_Q_ISS. Multiplication by
*              scalar during lifting (paper notation: \theta)
* 
* Arguments:   - poly_qiss_vec_k res: polynomial vector result of fac*theta(arg) modulo PARAM_Q_ISS (initialized)
*              - const poly_q arg: polynomial to be embedded
*              - const int64_t fac: scalar factor during lifting (1 or PARAM_Q1_ISS)
**************************************************/
void poly_qiss_subring_embed_vec_k(poly_qiss_vec_k res, const poly_q arg, const int64_t fac) {
  size_t i,j;
  coeff_q c_q;
  coeff_qiss c_qiss;
  for (i = 0; i < PARAM_K_ISS; i++) {
    for (j = 0; j < PARAM_N_ISS; j++) {
      c_q = poly_q_get_coeff(arg, i + j*PARAM_K_ISS);
      c_qiss = fac * c_q;
      poly_qiss_set_coeff(res->entries[i], j, c_qiss);
    }
  }
}

/*************************************************
* Name:        poly_qiss_subring_embed_mat_k_k
*
* Description: Subring multiplication matrix embedding and lifting mod PARAM_Q_ISS. Multiplication by
*              scalar during lifting (paper notation: M_{\theta})
* 
* Arguments:   - poly_qiss_mat_k_k res: polynomial matrix result of fac*M_theta(arg) modulo PARAM_Q_ISS (initialized)
*              - const poly_q arg: polynomial to be embedded
*              - const int64_t fac: scalar factor during lifting (1 or PARAM_Q1_ISS)
**************************************************/
void poly_qiss_subring_embed_mat_k_k(poly_qiss_mat_k_k res, const poly_q arg, const int64_t fac) {
  size_t i, j;
  poly_qiss_vec_k tmp_vec, tmp_vec_x;
  poly_qiss tmp;

  // init vectors and polynomials
  poly_qiss_vec_k_init(tmp_vec);
  poly_qiss_vec_k_init(tmp_vec_x);
  poly_qiss_init(tmp);

  // subring embedding vector
  poly_qiss_subring_embed_vec_k(tmp_vec, arg, fac);

  // multiplication by x (except the 0 index, not used)
  for (i = 1; i < PARAM_K_ISS; i++) {
    poly_qiss_mul_x(tmp_vec_x->entries[i], tmp_vec->entries[i]);
  }

  // constructing toeplitz matrix
  for (i = 0; i < PARAM_K_ISS; i++) {
    // subdiagonals
    for (j = 0; j < i + 1; j++) {
      poly_qiss_set(res->rows[i]->entries[j], tmp_vec->entries[i-j]);
    }
    // updiagonals
    for (j = i+1; j < PARAM_K_ISS; j++) {
      poly_qiss_set(res->rows[i]->entries[j], tmp_vec_x->entries[PARAM_K_ISS + (i-j)]);
    }
  }

  // clean up vectors and polynomials
  poly_qiss_vec_k_clear(tmp_vec);
  poly_qiss_vec_k_clear(tmp_vec_x);
  poly_qiss_clear(tmp);
}

/*************************************************
* Name:        challenge_size_iss
*
* Description: Compute Manhattan-like bound on ZKP round 4 challenge
*              namely |c^64|_1 ^ (1/64)
* 
* Arguments:   - const poly_qiss c: polynomial challenge
* 
* Returns the 64-th root of the l1 norm of c^64, i.e., |c^64|_1 ^ (1/64)
**************************************************/
uint64_t challenge_size_iss(const poly_qiss c) {
  size_t i, j;
  int exact;
  uint64_t out;
  coeff_qiss tmp_coeff;
  fmpz_t r, norm;

  // init flint objects (does not use poly_z because different ring degree)
  fmpz_init(r);
  fmpz_init(norm);
  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);

  // convert poly_qiss to fmpz_poly_t
  for (i = 0; i < PARAM_N_ISS; i++) {
    tmp_coeff = poly_qiss_get_coeff_centered(c, i);
    fmpz_poly_set_coeff_si(tmp, i, tmp_coeff);
  }

  // computing c^64
  for (i = 0; i < 6; i++) {
    fmpz_poly_mul(tmp, tmp, tmp);
    // modulo X^PARAM_N_ISS + 1 simply substracts the coeffs of degree N_ISS,...,2N_ISS-1
    // to those of degree 0,...,N_ISS-1
    for (j = 0; j < PARAM_N_ISS; j++) {
      // r and norm are used as temporary variables
      fmpz_poly_get_coeff_fmpz(r,  tmp, (signed long)i);
      fmpz_poly_get_coeff_fmpz(norm, tmp, (signed long)i + PARAM_N_ISS);
      fmpz_sub(r, r, norm);
      fmpz_poly_set_coeff_fmpz(tmp, (signed long)i, r);
      fmpz_poly_set_coeff_ui(tmp, (signed long)i + PARAM_N_ISS, 0);
    }
  }
  // tmp = c^64
  fmpz_zero(norm);
  for (i = 0; i < PARAM_N_ISS; i++) {
    // r is used as a temporary variable
    fmpz_poly_get_coeff_fmpz(r,  tmp, (signed long)i);
    fmpz_abs(r, r);
    fmpz_add(norm, norm, r);
  }
  exact = fmpz_root(r, norm, 64);
  out = ((uint64_t) fmpz_get_ui(r)) || (1 - exact); // adds one if not exact

  // clean up flint objects
  fmpz_clear(r);
  fmpz_clear(norm);
  fmpz_poly_clear(tmp);

  return out;
}

/*************************************************
* Name:        poly_qshow_subring_embed_vec_k
*
* Description: Subring embedding and lifting mod PARAM_Q_SHOW. Multiplication by
*              scalar during lifting (paper notation: \theta)
* 
* Arguments:   - poly_qshow_vec_k res: polynomial vector result of fac*theta(arg) modulo PARAM_Q_SHOW (initialized)
*              - const poly_q arg: polynomial to be embedded
*              - const int64_t fac: scalar factor during lifting (1 or PARAM_Q1_SHOW)
**************************************************/
void poly_qshow_subring_embed_vec_k(poly_qshow_vec_k res, const poly_q arg, const int64_t fac) {
  size_t i,j;
  coeff_q c_q;
  coeff_qshow c_qshow;
  for (i = 0; i < PARAM_K_SHOW; i++) {
    for (j = 0; j < PARAM_N_SHOW; j++) {
      c_q = poly_q_get_coeff_centered(arg, i + j*PARAM_K_SHOW);
      c_qshow = fac * c_q;
      poly_qshow_set_coeff(res->entries[i], j, c_qshow);
    }
  }
}

/*************************************************
* Name:        poly_qshow_subring_embed_mat_k_k
*
* Description: Subring multiplication matrix embedding and lifting mod PARAM_Q_SHOW. Multiplication by
*              scalar during lifting (paper notation: M_{\theta})
* 
* Arguments:   - poly_qshow_mat_k_k res: polynomial matrix result of fac*M_theta(arg) modulo PARAM_Q_SHOW (initialized)
*              - const poly_q arg: polynomial to be embedded
*              - const int64_t fac: scalar factor during lifting (1 or PARAM_Q1_SHOW)
**************************************************/
void poly_qshow_subring_embed_mat_k_k(poly_qshow_mat_k_k res, const poly_q arg, const int64_t fac) {
  size_t i, j;
  poly_qshow_vec_k tmp_vec, tmp_vec_x;
  poly_qshow tmp;

  // init vectors and polynomials
  poly_qshow_vec_k_init(tmp_vec);
  poly_qshow_vec_k_init(tmp_vec_x);
  poly_qshow_init(tmp);

  // subring embedding vector
  poly_qshow_subring_embed_vec_k(tmp_vec, arg, fac);

  // multiplication by x (except the 0 index, not used)
  for (i = 1; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul_x(tmp_vec_x->entries[i], tmp_vec->entries[i]);
  }

  // constructing toeplitz matrix
  for (i = 0; i < PARAM_K_SHOW; i++) {
    // subdiagonals
    for (j = 0; j < i + 1; j++) {
      poly_qshow_set(res->rows[i]->entries[j], tmp_vec->entries[i-j]);
    }
    // updiagonals
    for (j = i+1; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(res->rows[i]->entries[j], tmp_vec_x->entries[PARAM_K_SHOW + (i-j)]);
    }
  }

  // clean up vectors and polynomials
  poly_qshow_vec_k_clear(tmp_vec);
  poly_qshow_vec_k_clear(tmp_vec_x);
  poly_qshow_clear(tmp);
}

/*************************************************
* Name:        challenge_size_show
*
* Description: Compute Manhattan-like bound on ZKP round 4 challenge
*              namely |c^64|_1 ^ (1/64)
* 
* Arguments:   - const poly_qshow c: polynomial challenge
* 
* Returns the 64-th root of the l1 norm of c^64, i.e., |c^64|_1 ^ (1/64)
**************************************************/
uint64_t challenge_size_show(const poly_qshow c) {
  size_t i, j;
  int exact;
  uint64_t out;
  coeff_qshow tmp_coeff;
  fmpz_t r, norm;

  // init flint objects
  fmpz_init(r);
  fmpz_init(norm);
  fmpz_poly_t tmp;
  fmpz_poly_init(tmp);

  // convert poly_qshow to fmpz_poly_t
  for (i = 0; i < PARAM_N_SHOW; i++) {
    tmp_coeff = poly_qshow_get_coeff_centered(c, i);
    fmpz_poly_set_coeff_si(tmp, i, tmp_coeff);
  }

  // computing c^64
  for (i = 0; i < 6; i++) {
    fmpz_poly_mul(tmp, tmp, tmp);
    // modulo X^PARAM_N_SHOW + 1 simply substracts the coeffs of degree N_SHOW,...,2N_SHOW-1
    // to those of degree 0,...,N_SHOW-1
    for (j = 0; j < PARAM_N_SHOW; j++) {
      // r and norm are used as temporary variables
      fmpz_poly_get_coeff_fmpz(r,  tmp, (signed long)i);
      fmpz_poly_get_coeff_fmpz(norm, tmp, (signed long)i + PARAM_N_SHOW);
      fmpz_sub(r, r, norm);
      fmpz_poly_set_coeff_fmpz(tmp, (signed long)i, r);
      fmpz_poly_set_coeff_ui(tmp, (signed long)i + PARAM_N_SHOW, 0);
    }
  }
  // tmp = c^64
  fmpz_zero(norm);
  for (i = 0; i < PARAM_N_SHOW; i++) {
    // r is used as a temporary variable
    fmpz_poly_get_coeff_fmpz(r,  tmp, (signed long)i);
    fmpz_abs(r, r);
    fmpz_add(norm, norm, r);
  }
  exact = fmpz_root(r, norm, 64);
  out = ((uint64_t) fmpz_get_ui(r)) || (1 - exact); // adds one if not exact

  // clean up flint objects
  fmpz_clear(r);
  fmpz_clear(norm);
  fmpz_poly_clear(tmp);

  return out;
}
