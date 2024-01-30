#include "poly_q_sampling.h"
#include "fips202.h"
#include "arith.h"
#include "poly_real_vec_2d.h"
#include "random.h"

/*************************************************
* Name:        poly_q_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q deterministically from a seed.
* 
* Arguments:   - poly_q pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
**************************************************/
static void poly_q_uniform(poly_q pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, size_t i, size_t j) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, seed, SEED_BYTES);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_absorb(&state, (const uint8_t*) &i, sizeof(i));
  shake128_absorb(&state, (const uint8_t*) &j, sizeof(j));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_N) {
#if PARAM_Q_BITLEN > 20
#error "PARAM_Q_BITLEN too big for uniform sampling."
#else
    // idea: take 5 byte, divide them into two partitions each of 20 bits, potentially ignore the MSBs, perform rejection sampling
    if (bytecnt < 5) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    int64_t tmp5byte = (int64_t) (output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32));
    int64_t tmp = tmp5byte & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp = (tmp5byte >> 20) & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
    }

    off += 5;
    bytecnt -= 5;
#if PARAM_Q_BITLEN < 19
#warning "PARAM_Q_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D x PARAM_D modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_mat_d_d mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - uint8_t offset: offset domain separator for XOF
**************************************************/
void poly_q_mat_d_d_uniform(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_q_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j + offset);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D x PARAM_K modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_mat_d_k mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_mat_d_k_uniform(poly_q_mat_d_k mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_q_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_m_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D x PARAM_M modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_mat_d_m mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_mat_d_m_uniform(poly_q_mat_d_m mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_M; j++) {
      poly_q_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_q_vec_d_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_D modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_vec_d vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_vec_d_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i;
  for (i = 0; i < PARAM_D; i++) {
    poly_q_uniform(vec->entries[i], seed, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_binomial
*
* Description: Sample a centered binomial polynomial matrix of
*              size PARAM_D x PARAM_D with binomial parameter 1 
*              deterministically from a seed.
* 
* Arguments:   - poly_q_mat_d_d mat: output binomial polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t cnt: repetition domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_mat_d_d_binomial(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint64_t output[PARAM_N*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N/64];
  uint64_t coef_sign[PARAM_N/64];
  keccak_state state;
  size_t i,j,k;

  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      shake256_init(&state);
      shake256_absorb(&state, seed, SEED_BYTES);
      shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
      shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
      shake256_absorb(&state, (const uint8_t*) &i, sizeof(i));
      shake256_absorb(&state, (const uint8_t*) &j, sizeof(j));
      shake256_finalize(&state);
      shake256_squeeze((uint8_t*)output, PARAM_N*2/8, &state);

      for (k = 0; k < PARAM_N/64; k++) {
        coef_lsb[k] = output[2*k] ^ output[2*k+1];
        coef_sign[k] = output[2*k] & output[2*k+1];
      }
      for (k = 0; k < PARAM_N; k++) {
        poly_q_set_coeff(mat->rows[i]->entries[j], k, (int32_t)((coef_lsb[k/64] >> (k%64))&1) + (int32_t)(((coef_sign[k/64] >> ((k%64))) << 1)&2) - 1);
        // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
      }
    }
  }
}

/*************************************************
* Name:        poly_q_vec_d_bin_uniform
*
* Description: Sample a uniform binary polynomial vector of
*              size PARAM_D with deterministically from a seed.
* 
* Arguments:   - poly_q_vec_d vec: output binomial polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - uint8_t offset: offset domain separator for XOF
**************************************************/
void poly_q_vec_d_bin_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint8_t output[PARAM_N/8];
  keccak_state state;
  size_t i,j;

  for (i = 0; i < PARAM_D; i++) {
    j = i + offset;
    shake256_init(&state);
    shake256_absorb(&state, seed, SEED_BYTES);
    shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
    shake256_absorb(&state, (const uint8_t*) &j, sizeof(j));
    shake256_finalize(&state);
    shake256_squeeze(output, PARAM_N/8, &state);

    poly_q_from_bits(vec->entries[i], output);
  }
}

/*************************************************
* Name:        poly_q_vec_2d_dk_sample_perturb
*
* Description: Gaussian perturbation sampler (SamplePerturb)
* 
* Arguments:   - poly_q_vec_d *p1: array of polynomial vectors to host top perturbation (initialized)
*              - poly_q_vec_d *p2: array of polynomial vectors to host bottom perturbation (initialized)
*              - const poly_q_mat_d_d *R: array of polynomial matrices, secret key R 
*              - const poly_real_mat_2d_2d *S: array of polynomial matrices, Schur complements for sampling
**************************************************/
static void poly_q_vec_2d_dk_sample_perturb(poly_q_vec_d p1[2], poly_q_vec_d p2[PARAM_K], const poly_q_mat_d_d R[2][PARAM_K], const poly_real_mat_2d_2d S) {
  size_t i,j;
  poly_real tmp_sub, tmp_mul;
  poly_real_vec_2d c;

  // init vectors and polynomials
  poly_real_init(tmp_sub);
  poly_real_init(tmp_mul);
  poly_real_vec_2d_init(c);

  // sample p2 from Gaussian distribution with variance sqrt(s_2^2 - s_G^2)
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_gaussian_sqrt_s2sq_sGsq(p2[i]); // probabilistic
  }

  // p1 = R times p2  (p1 is only used as intermediate variable)
  poly_q_mat_d_d_mul_vec_d(p1[0], R[0][0], p2[0]);
  poly_q_mat_d_d_mul_vec_d(p1[1], R[1][0], p2[0]);
  for (i = 1; i < PARAM_K; i++) {
    poly_q_mat_d_d_muladd_vec_d(p1[0], R[0][i], p2[i]);
    poly_q_mat_d_d_muladd_vec_d(p1[1], R[1][i], p2[i]);
  }

  // convert p1 to c, which is a poly_real_vec_2d, and scale with -s_G^2/(s_2^2 - s_G^2)
  for (i = 0; i < PARAM_D; i++) {
    poly_real_from_poly_q(c->entries[i          ], p1[0]->entries[i]);
    poly_real_from_poly_q(c->entries[i + PARAM_D], p1[1]->entries[i]);
    poly_real_mul_scalar(c->entries[i], c->entries[i], PARAM_NEGSGSQ_DIV_S2SQ_SGSQ);
    poly_real_mul_scalar(c->entries[i + PARAM_D], c->entries[i + PARAM_D], PARAM_NEGSGSQ_DIV_S2SQ_SGSQ);
  }

  for (i = 2*PARAM_D; i > 0; i--) {
    // sample p1[(i-1)/PARAM_D]->entries[(i-1)%PARAM_D] from a discrete Gaussian distribution with
    // center c->entries[i-1] and variance M_tau(S->rows[i-1]->entries[i-1])
    poly_q_samplefz(p1[(i-1)/PARAM_D]->entries[(i-1)%PARAM_D], S->rows[i-1]->entries[i-1], c->entries[i-1]); // probabilistic

    // update the entries from c with indices j=0..i-2 by adding S->rows[j]->entries[i-1] (this already contains f_i^-1)
    poly_real_sub_poly_real_poly_q(tmp_sub, p1[(i-1)/PARAM_D]->entries[(i-1)%PARAM_D], c->entries[i-1]);
    for (j = 0; j < i-1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i-1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  //  clean up vectors and polynomials
  poly_real_vec_2d_clear(c);
  poly_real_clear(tmp_mul);
  poly_real_clear(tmp_sub);
}

#if PARAM_K == 5
// Gaussian widths for Klein Sampler
static const double klein_widths[PARAM_K] = {
  3.42997312037747414948L,
  3.43866739136833965418L,
  3.43871169237754559234L,
  3.43871191840160950193L,
  4.35451735086351643389L
};
// Scaled Gram-Schmidt of the basis of the gadget lattice (must be negated)
static const double neg_scaled_gadget_gso[PARAM_K][PARAM_K] = {
  {-0.07106598984771574090, -0.00507601066998161245, -0.00036257214280532797, -0.00002589801018292152, -0.00000234851491659608},
  { 0.00507614213197969504, -0.07106414937974257773, -0.00507600999927459128, -0.00036257214256090126, -0.00003287920883229431},
  { 0.00000000000000000000,  0.00510190868360396732, -0.07106413998984428826, -0.00507600999585261763, -0.00046030892365212574},
  { 0.00000000000000000000,  0.00000000000000000000,  0.00510204014218007627, -0.07106413994193663819, -0.00644432493112979867},
  { 0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000,  0.00510204081288700810, -0.09022054903581719354}
};
// Gadget lattice basis 
static const int64_t gadget_basis[PARAM_K][PARAM_K] = {
  {14,   0,   0,   0,   5},
  {-1,  14,   0,   0,   6},
  { 0,  -1,  14,   0,   2},
  { 0,   0,  -1,  14,   1},
  { 0,   0,   0,  -1,  11}
};
#else
#error "Some constants are not generated for the given PARAM_K."
#endif

/*************************************************
* Name:        sample_klein
*
* Description: Gaussian Klein sampler on the gadget lattice for SamplePre
* 
* Arguments:   - poly_q_vec_d *v: array of polynomial vectors to host sample (initialized)
*              - const poly_q_vec_d w: polynomial vector, center
**************************************************/
static void sample_klein(poly_q_vec_d v[PARAM_K], const poly_q_vec_d w) {
  int64_t i,j,i_1,i_2,i_2_1,i_2_2;
  double di;
  int64_t zi;
  coeff_q w_coeff, v_coeff;

  for (i = PARAM_N * PARAM_D * PARAM_K; i > 0; i--) {
    i_1   = (i - 1) / (PARAM_N * PARAM_D);
    i_2   = (i - 1) % (PARAM_N * PARAM_D);
    i_2_1 = i_2 / PARAM_N;
    i_2_2 = i_2 % PARAM_N;

    // di = <v - [-w|0|0...], Btilde'[i]>
    w_coeff = poly_q_get_coeff(w->entries[i_2_1], (size_t) i_2_2);
    // coeffs of w do not have to be in centered representation
    if (i == PARAM_N * PARAM_D * PARAM_K) {
      di = neg_scaled_gadget_gso[0][i_1] * ((double) w_coeff);
      // "below" w, everything is zero and does not contribute to the dot product
    } else {
      // coeffs of v have to be in centered representation to compute di
      v_coeff = poly_q_get_coeff_centered(v[0]->entries[i_2_1], (size_t) i_2_2);
      di = neg_scaled_gadget_gso[0][i_1] * ((double) (w_coeff + v_coeff));
      for (j = 1; j < PARAM_K; j++) {
        v_coeff = poly_q_get_coeff_centered(v[j]->entries[i_2_1], (size_t) i_2_2);
        di += neg_scaled_gadget_gso[j][i_1] * ((double) v_coeff);
      }
    }

    // zi following D_{Z, si, di}
    zi = SampleZ(di, klein_widths[i_1]);

    // v = v + zi*B[i]
    for (j = 0; j < PARAM_K; j++) {
      v_coeff = poly_q_get_coeff(v[j]->entries[i_2_1], (size_t) i_2_2);
      v_coeff += zi * gadget_basis[j][i_1];
      poly_q_set_coeff(v[j]->entries[i_2_1], (size_t) i_2_2, v_coeff);
    }
  }
}

/*************************************************
* Name:        poly_q_vec_2d_dk_sample_pre
*
* Description: Elliptic Gaussian preimage sampler (SamplePre)
* 
* Arguments:   - poly_q_vec_d *v1: array of polynomial vectors to host top preimage (initialized)
*              - poly_q_vec_d *v2: array of polynomial vectors to host bottom preimage (initialized)
*              - const poly_q_mat_d_d *R: array of polynomial matrices, secret key R 
*              - const poly_q_mat_d_d A: polynomial matrix, public A'
*              - const poly_q_mat_d_d *B: array of polynomial matrices, public B
*              - const poly_q_vec_d u: polynomial vector, public u
*              - const poly_q tag: polynomial, tag
*              - const poly_real_mat_2d_2d *S: array of polynomial matrices, Schur complements for sampling
**************************************************/
void poly_q_vec_2d_dk_sample_pre(
    poly_q_vec_d v1[2], 
    poly_q_vec_d v2[PARAM_K], 
    const poly_q_mat_d_d R[2][PARAM_K], 
    const poly_q_mat_d_d A, 
    const poly_q_mat_d_d B[PARAM_K], 
    const poly_q_vec_d u, 
    const poly_q tag, 
    const poly_real_mat_2d_2d S) {
  size_t i;
  poly_q_vec_d p1[2];
  poly_q_vec_d p2[PARAM_K];
  poly_q taginv;
  poly_q_vec_d w, tmp, y[PARAM_K];
  int64_t bexpi;

  // init vectors and polynomials
  poly_q_init(taginv);
  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
  poly_q_vec_d_init(p1[0]);
  poly_q_vec_d_init(p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(p2[i]);
    poly_q_vec_d_init(y[i]);
  }

  // sample_perturb
  poly_q_vec_2d_dk_sample_perturb(p1, p2, R, S);

  // w = tag^{-1} * (u - [A | -B]*p)  - G*p_2
  // = tag^-1 * (u - p1[0] - A*p1[1] + B * p2) - G*p2
  poly_q_invert(taginv, tag);
  poly_q_vec_d_sub(w, u, p1[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p2[i]);
  }
  poly_q_vec_d_mul_poly(w, w, taginv);
  poly_q_vec_d_sub(w, w, p2[0]);
  bexpi = 1;
  for (i = 1; i < PARAM_K; i++) {
    bexpi *= PARAM_B;
    poly_q_vec_d_mul_scalar(tmp, p2[i], bexpi);
    poly_q_vec_d_sub(w, w, tmp);
  }

  // c = G^-1(w) = [w|0|0|...|0], but we only store and pass w
  // sample y from Klein sampler
  sample_klein(y, w);

  // z = c + y; z is represented by y now
  poly_q_vec_d_add(y[0], y[0], w);

  // v1 = p1 + Rz
  poly_q_mat_d_d_mul_vec_d(v1[0], R[0][0], y[0]);
  poly_q_mat_d_d_mul_vec_d(v1[1], R[1][0], y[0]);
  poly_q_vec_d_add(v1[0], v1[0], p1[0]);
  poly_q_vec_d_add(v1[1], v1[1], p1[1]);
  for (i = 1; i < PARAM_K; i++) {
    poly_q_mat_d_d_muladd_vec_d(v1[0], R[0][i], y[i]);
    poly_q_mat_d_d_muladd_vec_d(v1[1], R[1][i], y[i]);
  }

  // v2 = p2 + z
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_add(v2[i], p2[i], y[i]);
  }

  // clean up vectors and polynomials
  poly_q_clear(taginv);
  poly_q_vec_d_clear(w);
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(p1[0]);
  poly_q_vec_d_clear(p1[1]);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(p2[i]);
    poly_q_vec_d_clear(y[i]);
  }
}

/*************************************************
* Name:        poly_q_binary_fixed_weight
*
* Description: Sample binary polynomial with fixed Hamming weight PARAM_W
* 
* Arguments:   - poly_q res: output uniform binary polynomial with fixed Hamming weight (initialized)
*              - uint8_t *state_in: state from which to expand the polynomial (allocated STATE_BYTES bytes)
**************************************************/
void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]) {
  unsigned int i, b, pos = 0;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_absorb_once(&state, state_in, STATE_BYTES);
  shake256_squeezeblocks(buf, 1, &state);

  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, 0);
  }
  for (i = PARAM_N - PARAM_W; i < PARAM_N; i++) {
    do {
      if (pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }
      b = buf[pos++];
    } while (b > i);

    poly_q_set_coeff(res, i, poly_q_get_coeff(res, b));
    poly_q_set_coeff(res, b, 1);
  }
}
