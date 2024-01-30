#include "poly_qshow_sampling.h"
#include "fips202.h"
#include "random.h"

/*************************************************
* Name:        poly_qshow_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q_SHOW deterministically from a seed.
* 
* Arguments:   - poly_qshow pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
**************************************************/
static void poly_qshow_uniform(poly_qshow pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, size_t i, size_t j) {
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
  while (cnt < PARAM_N_SHOW) {
#if PARAM_Q_SHOW_BITLEN > 60
#error "PARAM_Q_SHOW_BITLEN too big for uniform sampling."
#else
    // idea: take 15 bytes, ignore MSBs
    if (bytecnt < 15) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp8byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32) | ((uint64_t)output[off+5] << 40) | ((uint64_t)output[off+6] << 48) | ((uint64_t)output[off+7] << 56);
    uint64_t tmp = tmp8byte & ((1UL<<PARAM_Q_SHOW_BITLEN)-1);
    if (tmp < PARAM_Q_SHOW) {
      poly_qshow_set_coeff(pout, cnt++, tmp);
    }
    if (cnt >= PARAM_N_SHOW) {
      break;
    }
    tmp8byte >>= 60;
    tmp8byte |= (output[off+8] | ((uint64_t)output[off+9] << 8) | ((uint64_t)output[off+10] << 16) | ((uint64_t)output[off+11] << 24) | ((uint64_t)output[off+12] << 32) | ((uint64_t)output[off+13] << 40) | ((uint64_t)output[off+14] << 48)) << 4;
    tmp = tmp8byte & ((1UL<<PARAM_Q_SHOW_BITLEN)-1);
    if (tmp < PARAM_Q_SHOW) {
      poly_qshow_set_coeff(pout, cnt++, tmp);
    }
    off += 15;
    bytecnt -= 15;
#if PARAM_Q_SHOW_BITLEN < 58
#warning "PARAM_Q_SHOW_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        vec_qshow_uniform
*
* Description: Sample a uniformly random vector of integer modulo
*              PARAM_Q_SHOW of size PARAM_ARP_SHOW + 6 deterministically 
*              from an input buffer.
* 
* Arguments:   - coeff_qshow *out: output uniform vector (allocated PARAM_ARP_SHOW+6 coeff_qshow)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - const uint32_t domain_separator: domain separator for XOF
*              - const uint32_t counter: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void vec_qshow_uniform(coeff_qshow out[PARAM_ARP_SHOW + 6], const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, buf, buflen);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_absorb(&state, (const uint8_t*)&counter, sizeof(uint32_t));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_ARP_SHOW + 6) {
#if PARAM_Q_SHOW_BITLEN > 60
#error "PARAM_Q_SHOW_BITLEN too big for uniform sampling."
#else
    // idea: take 15 bytes, ignore MSBs
    if (bytecnt < 15) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp8byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32) | ((uint64_t)output[off+5] << 40) | ((uint64_t)output[off+6] << 48) | ((uint64_t)output[off+7] << 56);
    uint64_t tmp = tmp8byte & ((1UL<<PARAM_Q_SHOW_BITLEN)-1);
    if (tmp < PARAM_Q_SHOW) {
      out[cnt++] = tmp;
    }
    if (cnt >= PARAM_ARP_SHOW + 6) {
      break;
    }
    tmp8byte >>= 60;
    tmp8byte |= (output[off+8] | ((uint64_t)output[off+9] << 8) | ((uint64_t)output[off+10] << 16) | ((uint64_t)output[off+11] << 24) | ((uint64_t)output[off+12] << 32) | ((uint64_t)output[off+13] << 40) | ((uint64_t)output[off+14] << 48)) << 4;
    tmp = tmp8byte & ((1UL<<PARAM_Q_SHOW_BITLEN)-1);
    if (tmp < PARAM_Q_SHOW) {
      out[cnt++] = tmp;
    }
    off += 15;
    bytecnt -= 15;
#if PARAM_Q_SHOW_BITLEN < 58
#warning "PARAM_Q_SHOW_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_qshow_uniform_but_zero
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q_SHOW with zero constant coefficient,
*              deterministically from a seed.
* 
* Arguments:   - poly_qshow out: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t kappa: rejection domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qshow_uniform_but_zero(poly_qshow out, const uint8_t seed[SEED_BYTES], uint32_t kappa, uint32_t domain_separator) {
  // TODO implement dedicated function?
  poly_qshow_uniform(out, seed, domain_separator, kappa, 0);
  poly_qshow_set_coeff(out, 0, 0);
}

/*************************************************
* Name:        poly_qshow_mat_d_m1_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D_SHOW x PARAM_M1_SHOW modulo PARAM_Q_SHOW 
*              deterministically from a seed.
* 
* Arguments:   - poly_qshow_mat_d_m1 mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - uint8_t offset: offset domain separator for XOF
**************************************************/
void poly_qshow_mat_d_m1_uniform(poly_qshow_mat_d_m1 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D_SHOW; i++) {
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      poly_qshow_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_m2_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D_SHOW x PARAM_M2_SHOW modulo PARAM_Q_SHOW 
*              deterministically from a seed.
* 
* Arguments:   - poly_qshow_mat_d_m2 mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qshow_mat_d_m2_uniform(poly_qshow_mat_d_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D_SHOW; i++) {
    for (j = 0; j < PARAM_M2_SHOW; j++) {
      poly_qshow_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_256l_m2_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_SHOW modulo PARAM_Q_SHOW 
*              deterministically from a seed.
* 
* Arguments:   - poly_qshow_mat_256l_m2 mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qshow_mat_256l_m2_uniform(poly_qshow_mat_256l_m2 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
    for (j = 0; j < PARAM_M2_SHOW; j++) {
      poly_qshow_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_qshow_vec_m2_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_M2_SHOW modulo PARAM_Q_SHOW 
*              deterministically from a seed.
* 
* Arguments:   - poly_qshow_vec_m2 vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qshow_vec_m2_uniform(poly_qshow_vec_m2 vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i;
  for (i = 0; i < PARAM_M2_SHOW; i++) {
    poly_qshow_uniform(vec->entries[i], seed, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_qshow_vec_l_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_L_SHOW modulo PARAM_Q_SHOW 
*              deterministically from an input buffer.
* 
* Arguments:   - poly_qshow_vec_l vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qshow_vec_l_uniform(poly_qshow_vec_l vec, const uint8_t *buf, uint32_t domain_separator, size_t buflen) {
  size_t i;
  uint8_t lazybuf[SEED_BYTES];
  sha3_256(lazybuf, buf, buflen); // TODO remove this and use buf directly for sampling input
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_uniform(vec->entries[i], lazybuf, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_qshow_vec_k_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_K_SHOW modulo PARAM_Q_SHOW 
*              deterministically from an input buffer.
* 
* Arguments:   - poly_qshow_vec_k vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - uint32_t cnt: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qshow_vec_k_uniform(poly_qshow_vec_k vec, const uint8_t *buf, uint32_t domain_separator, uint32_t cnt, size_t buflen) {
  size_t i;
  uint8_t lazybuf[SEED_BYTES];
  sha3_256(lazybuf, buf, buflen); // TODO remove this and use buf directly for sampling input
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_uniform(vec->entries[i], lazybuf, domain_separator, i, cnt);
  }
}

// TODO streamline the binomial sampling function signatures
/*************************************************
* Name:        poly_qshow_vec_m1_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_M1_SHOW with binomial parameter 1 
*              deterministically from a buffer.
* 
* Arguments:   - poly_qshow_vec_m1 res: output binomial polynomial vector (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
*              - size_t inlen: length of the input buffer
**************************************************/
void poly_qshow_vec_m1_binomial(poly_qshow_vec_m1 res, const uint8_t *buf, uint32_t domain_separator, uint32_t i, size_t inlen) {
#if (PARAM_N_SHOW%64) != 0
#error "PARAM_N_SHOW must be divisible by 64"
#endif
  uint64_t output[PARAM_M1_SHOW*PARAM_N_SHOW*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N_SHOW/64];
  uint64_t coef_sign[PARAM_N_SHOW/64];
  keccak_state state;
  size_t j,k;
  shake256_init(&state);
  shake256_absorb(&state, buf, inlen);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*) &i, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, PARAM_M1_SHOW*PARAM_N_SHOW*2/8, &state);
  for (j = 0; j < PARAM_M1_SHOW; j++) {
    for (k = 0; k < PARAM_N_SHOW/64; k++) {
      coef_lsb[k]  = output[2*k + 2*PARAM_N_SHOW/64*j] ^ output[2*k+1 + 2*PARAM_N_SHOW/64*j];
      coef_sign[k] = output[2*k + 2*PARAM_N_SHOW/64*j] & output[2*k+1 + 2*PARAM_N_SHOW/64*j];
    }
    for (k = 0; k < PARAM_N_SHOW; k++) {
      poly_qshow_set_coeff(res->entries[j], k, (int32_t)((coef_lsb[k/64] >> (k%64))&1) + (int32_t)(((coef_sign[k/64] >> ((k%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
}

/*************************************************
* Name:        poly_qshow_vec_m2_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_M2_SHOW with binomial parameter 1 
*              deterministically from a buffer.
* 
* Arguments:   - poly_qshow_vec_m2 res: output binomial polynomial vector (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input (allocated SEED_BYTES bytes)
*              - const uint32_t cnt: rejection domain separator for XOF
*              - const uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qshow_vec_m2_binomial(poly_qshow_vec_m2 res, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator) {
#if (PARAM_N_SHOW%64) != 0
#error "PARAM_N_SHOW must be divisible by 64"
#endif
  uint64_t output[PARAM_M2_SHOW*PARAM_N_SHOW*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N_SHOW/64];
  uint64_t coef_sign[PARAM_N_SHOW/64];
  keccak_state state;
  size_t k,l;
  shake256_init(&state);
  shake256_absorb(&state, buf, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, PARAM_M2_SHOW*PARAM_N_SHOW*2/8, &state);
  for (k = 0; k < PARAM_M2_SHOW; k++) {
    for (l = 0; l < PARAM_N_SHOW/64; l++) {
      coef_lsb[l]  = output[2*l+2*PARAM_N_SHOW/64*k] ^ output[2*l+1+2*PARAM_N_SHOW/64*k];
      coef_sign[l] = output[2*l+2*PARAM_N_SHOW/64*k] & output[2*l+1+2*PARAM_N_SHOW/64*k];
    }
    for (l = 0; l < PARAM_N_SHOW; l++) {
      poly_qshow_set_coeff(res->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
}

/*************************************************
* Name:        poly_qshow_vec_m1_sample_gaussian_s1
*
* Description: Sample a polynomial vector with PARAM_M1_SHOW entries
*              from the centered spherical Gaussian with parameter
*              PARAM_S1_SHOW
* 
* Arguments:   - poly_qshow_vec_m1 res: the polynomial to host the Gaussian sample
**************************************************/
void poly_qshow_vec_m1_sample_gaussian_s1(poly_qshow_vec_m1 res) {
  for (size_t i = 0; i < PARAM_M1_SHOW; i++) {
    for (size_t j = 0; j < PARAM_N_SHOW; j++) {
      poly_qshow_set_coeff(res->entries[i], j, SampleZ(0, PARAM_S1_SHOW));
    }
  }
}

/*************************************************
* Name:        poly_qshow_vec_m2_sample_gaussian_s2
*
* Description: Sample a polynomial vector with PARAM_M2_SHOW entries
*              from the centered spherical Gaussian with parameter
*              PARAM_S2_SHOW
* 
* Arguments:   - poly_qshow_vec_m2 res: the polynomial to host the Gaussian sample
**************************************************/
void poly_qshow_vec_m2_sample_gaussian_s2(poly_qshow_vec_m2 res) {
  for (size_t i = 0; i < PARAM_M2_SHOW; i++) {
    for (size_t j = 0; j < PARAM_N_SHOW; j++) {
      poly_qshow_set_coeff(res->entries[i], j, SampleZ(0, PARAM_S2_SHOW));
    }
  }
}

/*************************************************
* Name:        poly_qshow_sample_challenge
*
* Description: Sample a uniformly random polynomial with 
*              coefficients between [-PARAM_RHO_SHOW, PARAM_RHO_SHOW]
*              and self-adjoint deterministically from an input buffer.
* 
* Arguments:   - poly_qshow out: output uniform bounded self-adjoint polynomial (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - const uint32_t domain_separator: domain separator for XOF
*              - const uint32_t counter: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qshow_sample_challenge(poly_qshow out, const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen) {
  uint8_t outbuf[SHAKE256_RATE];
  size_t outcnt = 0, bytecnt, pos = 0;
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, buf, buflen);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&counter, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeezeblocks(outbuf, 1, &state);
  bytecnt = SHAKE256_RATE;
  
#if PARAM_RHO_SHOW != 8
#error "poly_qshow_sample_challenge is implemented specifically for PARAM_RHO_SHOW = 8"
#endif
  // the idea: for each coefficient sample a byte and reject the byte iff it is 255
  // if accepted, the coefficient is (byte%17) - 8
  while (outcnt < PARAM_N_SHOW/2) {
    if (bytecnt == 0) {
      shake256_squeezeblocks(outbuf, 1, &state);
      bytecnt = SHAKE256_RATE;
      pos = 0;
    }
    if (outbuf[pos] < 255) {
      coeff_qshow tmp = ((coeff_qshow)(outbuf[pos]%17)) - 8;
      poly_qshow_set_coeff(out, outcnt, tmp);
      if (outcnt > 0) {
        poly_qshow_set_coeff(out, PARAM_N_SHOW-outcnt, -tmp);
      }
      outcnt += 1;
    }
    bytecnt -= 1;
    pos += 1;
  }
}
