#include "covariance.h"
#include "arith.h"

#include <flint/arb.h>
#include <flint/acb_dft.h>

/* Precision for FFT in spectral norm estimation */
#define SPECTRAL_NORM_PRECISION 53

static acb_dft_pre_t pre; // for fft precomputation

/*************************************************
* Name:        fft_precomp_setup
*
* Description: Initialize flint fft precomputation
**************************************************/
void fft_precomp_setup(void) {
  acb_dft_precomp_init(pre, 2*PARAM_N, SPECTRAL_NORM_PRECISION);
}

/*************************************************
* Name:        fft_precomp_teardown
*
* Description: Clear flint fft precomputation
**************************************************/
void fft_precomp_teardown(void) {
  acb_dft_precomp_clear(pre);
}

/*************************************************
* Name:        compute_covariance
*
* Description: Computes Schur complements of covariance for perturbation 
*              sampling from the signer's secret key
*
* Arguments:   - poly_real_mat_2d_2d S: polynomial matrix to host the sampling materials
*              - const poly_q_mat_d_d *RRstar: array of polynomial matrices containing R.R*
**************************************************/
void compute_covariance(poly_real_mat_2d_2d S, const poly_q_mat_d_d RRstar[2][2]) {
  size_t i,j,k;
  poly_real finv, tmppoly;

  // init polynomials
  poly_real_init(finv);
  poly_real_init(tmppoly);

  // convert RRstar to poly_real_mat_2d_2d and write in S (RRstar needs to be read in centered representation!)
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_from_poly_q(S->rows[i]->entries[j], RRstar[i/PARAM_D][j/PARAM_D]->rows[i%PARAM_D]->entries[j%PARAM_D]);
    }
  }

  // multiply by constant (s_G^-2 - s_2^-2)^-1, subtract from s_1^2 times I_2d
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_mul_scalar(S->rows[i]->entries[j], S->rows[i]->entries[j], PARAM_SGINVSQ_S2INVSQ); // PARAM_SGINVSQ_S2INVSQ = - (s_G^-2 - s_2^-2)^-1  note the leading minus
      if (i == j) {
        poly_real_add_constant(S->rows[i]->entries[j], PARAM_S1SQ); // PARAM_S1SQ = s_1^2
      }
    }
  }

  for (i = 2*PARAM_D-1; i > 0; i--) {
    // i-th row and i-th column as well as everything beyond stands as it is

    // invert f_i
    poly_real_invert(finv, S->rows[i]->entries[i]);

    // multiply finv to si (i-th column until excluding i-th polynomial, which is f_i)
    for (j = 0; j < i; j++) {
      poly_real_mul(S->rows[j]->entries[i], S->rows[j]->entries[i], finv);
    }

    // update S->rows[j]->entries[k] with j=0..i-1, k=0..i-1 by subtracting S->rows[j]->entries[i] times S->rows[i]->entries[k]
    // note that S->rows[j]->entries[i] has finv already multiplied
    for (j = 0; j < i; j++) {
      for (k = 0; k < i; k++) {
        poly_real_mul(tmppoly, S->rows[j]->entries[i], S->rows[i]->entries[k]);
        poly_real_sub(S->rows[j]->entries[k], S->rows[j]->entries[k], tmppoly);
      }
    }
  }
  
  // clean up polynomials
  poly_real_clear(finv);
  poly_real_clear(tmppoly);
}

/*************************************************
* Name:        sk_sq_spectral_norm
*
* Description: Computes an estimate to the square spectral norm of R
*              based on the iterated power method and computes R.R*
*
* Arguments:   - poly_q_mat_d_d *RRstar: array of polynomial matrices hosting R.R*
*              - const poly_q_mat_d_d *R: polynomial matrix containing the secret key R
* 
* Returns a double corresponding to the estimated spectral norm of R
**************************************************/
double sk_sq_spectral_norm(poly_q_mat_d_d RRstar[2][2], const poly_q_mat_d_d R[2][PARAM_K]) {
  size_t i,j,k;
  acb_struct num_vec[PARAM_N*2], den_vec[PARAM_N*2];
  poly_q_mat_d_d Rstar[PARAM_K][2];
  poly_z_mat_d_d RRstar_z[2][2];
  poly_z_mat_d_d power_mat[2][2];
  poly_z_mat_d_d tmp_mat[2][2];
  poly_z_mat_d_d tmp_z;
  poly_q_mat_d_d tmp_q;
  fmpz_t c;
  arb_t num_embedding, den_embedding;
  double sq_norm_R, test_norm;

  // init flint objects for fft
  for (i = 0; i < PARAM_N*2; i++) {
    acb_init(&num_vec[i]);
    acb_init(&den_vec[i]);
  }
  arb_init(num_embedding);
  arb_init(den_embedding);
  fmpz_init(c);
  
  // init matrices
  poly_z_mat_d_d_init(RRstar_z[0][0]);
  poly_z_mat_d_d_init(RRstar_z[0][1]);
  poly_z_mat_d_d_init(RRstar_z[1][0]);
  poly_z_mat_d_d_init(RRstar_z[1][1]);
  poly_z_mat_d_d_init(power_mat[0][0]);
  poly_z_mat_d_d_init(power_mat[0][1]);
  poly_z_mat_d_d_init(power_mat[1][0]);
  poly_z_mat_d_d_init(power_mat[1][1]);
  poly_z_mat_d_d_init(tmp_mat[0][0]);
  poly_z_mat_d_d_init(tmp_mat[0][1]);
  poly_z_mat_d_d_init(tmp_mat[1][0]);
  poly_z_mat_d_d_init(tmp_mat[1][1]);
  poly_z_mat_d_d_init(tmp_z);
  poly_q_mat_d_d_init(tmp_q);
  // init and compute R*
  for (i = 0; i < 2; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_init(Rstar[j][i]);
      poly_q_mat_d_d_conjugate(Rstar[j][i], R[i][j]);
    }
  }

  // multiply R by R*
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[i][j], R[i][0], Rstar[0][j]);
      for (k = 1; k < PARAM_K; k++) {
        poly_q_mat_d_d_mul_mat_d_d(tmp_q, R[i][k], Rstar[k][j]);
        poly_q_mat_d_d_add(RRstar[i][j], RRstar[i][j], tmp_q);
      }
      // convert RRstar to poly_z_mat_d_d[2][2] and write in RRstar_z (needs to be read in centered representation!)
      poly_z_mat_d_d_from_poly_q_mat_d_d(RRstar_z[i][j], RRstar[i][j]);
      poly_z_mat_d_d_set(power_mat[i][j], RRstar_z[i][j]);
    }
  }

  // Iterated power: compute (R.R*)^32 and store in power_mat and tmp_mat
  for (k = 0; k < PARAM_IT_SPEC_NORM; k++) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        poly_z_mat_d_d_mul_mat_d_d(tmp_mat[i][j], power_mat[i][0], power_mat[0][j]);
        poly_z_mat_d_d_mul_mat_d_d(tmp_z, power_mat[i][1], power_mat[1][j]);
        poly_z_mat_d_d_add(tmp_mat[i][j], tmp_mat[i][j], tmp_z);
      }
    }
    poly_z_mat_d_d_set(power_mat[0][0], tmp_mat[0][0]);
    poly_z_mat_d_d_set(power_mat[0][1], tmp_mat[0][1]);
    poly_z_mat_d_d_set(power_mat[1][0], tmp_mat[1][0]);
    poly_z_mat_d_d_set(power_mat[1][1], tmp_mat[1][1]);
  }

  // multiplying power_mat by R.R*. power_mat contains (R.R*)^33 and tmp_mat contains (R.R*)^32
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_z_mat_d_d_mul_mat_d_d(power_mat[i][j], tmp_mat[i][0], RRstar_z[0][j]);
      poly_z_mat_d_d_mul_mat_d_d(tmp_z, tmp_mat[i][1], RRstar_z[1][j]);
      poly_z_mat_d_d_add(power_mat[i][j], power_mat[i][j], tmp_z);
    }
  }
  
  // fft of (RR^*)^32[0][0] and (RR^*)^33[0][0]
  for (i = 0; i < PARAM_N; i++) {
    poly_z_get_coeff(c, power_mat[0][0]->rows[0]->entries[0], i);
    acb_set_fmpz(&num_vec[i], c);
    poly_z_get_coeff(c, tmp_mat[0][0]->rows[0]->entries[0], i);
    acb_set_fmpz(&den_vec[i], c);
  }
  acb_dft_precomp(num_vec, num_vec, pre, SPECTRAL_NORM_PRECISION);
  acb_dft_precomp(den_vec, den_vec, pre, SPECTRAL_NORM_PRECISION);

  // finding maximal embedding ratio (spectral norm estimate)
  sq_norm_R = 0.0;
  for (i = 0; i < PARAM_N; i++) {
    acb_get_real(num_embedding, &num_vec[2*i+1]);
    acb_get_real(den_embedding, &den_vec[2*i+1]);
    arb_div(num_embedding, num_embedding, den_embedding, SPECTRAL_NORM_PRECISION);
    test_norm = arf_get_d(arb_midref(num_embedding), ARF_RND_DOWN); // round to arb midpoint

    // maximum test
    if (test_norm > sq_norm_R) {
      sq_norm_R = test_norm;
    }
  }

  // cleanup flint objects for fft
  for (i = 0; i < 2*PARAM_N; i++) {
    acb_clear(&num_vec[i]);
    acb_clear(&den_vec[i]);
  }
  arb_clear(num_embedding);
  arb_clear(den_embedding);
  fmpz_clear(c);

  // cleanup matrices
  for (i = 0; i < 2; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_clear(Rstar[j][i]);
    }
  }
  poly_z_mat_d_d_clear(RRstar_z[0][0]);
  poly_z_mat_d_d_clear(RRstar_z[0][1]);
  poly_z_mat_d_d_clear(RRstar_z[1][0]);
  poly_z_mat_d_d_clear(RRstar_z[1][1]);
  poly_z_mat_d_d_clear(power_mat[0][0]);
  poly_z_mat_d_d_clear(power_mat[0][1]);
  poly_z_mat_d_d_clear(power_mat[1][0]);
  poly_z_mat_d_d_clear(power_mat[1][1]);
  poly_z_mat_d_d_clear(tmp_mat[0][0]);
  poly_z_mat_d_d_clear(tmp_mat[0][1]);
  poly_z_mat_d_d_clear(tmp_mat[1][0]);
  poly_z_mat_d_d_clear(tmp_mat[1][1]);
  poly_z_mat_d_d_clear(tmp_z);
  poly_q_mat_d_d_clear(tmp_q);

  return sq_norm_R;
}
