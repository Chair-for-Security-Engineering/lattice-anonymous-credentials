#include "arith.h"
#include "randombytes.h"
#include "sep.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "covariance.h"

/*************************************************
* Name:        sep_keys_init
*
* Description: Initialize keys of signature scheme (SEP)
*              by calling Flint initialization
*
* Arguments:   - sep_pk_t *pk: pointer to public key structure
*              - sep_sk_t *sk: pointer to secret key structure
**************************************************/
void sep_keys_init(sep_pk_t *pk, sep_sk_t *sk) {
  size_t i, j;
  for (i = 0; i < PARAM_K; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_init(sk->R[j][i]);
    }
    poly_q_mat_d_d_init(pk->B[i]);
  }
  poly_real_mat_2d_2d_init(sk->S);
}

/*************************************************
* Name:        sep_keys_clear
*
* Description: Clear keys of signature scheme (SEP)
*              by calling Flint clean up
*
* Arguments:   - sep_pk_t *pk: pointer to public key structure
*              - sep_sk_t *sk: pointer to secret key structure
**************************************************/
void sep_keys_clear(sep_pk_t *pk, sep_sk_t *sk) {
  size_t i, j;
  for (i = 0; i < PARAM_K; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_clear(sk->R[j][i]);
    }
    poly_q_mat_d_d_clear(pk->B[i]);
  }
  poly_real_mat_2d_2d_clear(sk->S);
}

/*************************************************
* Name:        sep_sig_init
*
* Description: Initialize structure to host the signature (SEP)
*              by calling Flint initialization
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure
**************************************************/
void sep_sig_init(sep_sig_t *sig) {
  size_t i;
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(sig->v2[i]);
  }
  poly_q_vec_d_init(sig->v12);
  poly_q_vec_k_init(sig->v3);
  poly_q_init(sig->tag);
}

/*************************************************
* Name:        sep_sig_clear
*
* Description: Clear structure that hosts the signature (SEP)
*              by calling Flint clean up
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure
**************************************************/
void sep_sig_clear(sep_sig_t *sig) {
  size_t i;
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(sig->v2[i]);
  }
  poly_q_vec_d_clear(sig->v12);
  poly_q_vec_k_clear(sig->v3);
  poly_q_clear(sig->tag);
}

/*************************************************
* Name:        sep_keygen
*
* Description: Generates public and private key
*
* Arguments:   - sep_pk_t *pk: pointer to public key structure (initialized)
*              - sep_sk_t *sk: pointer to secret key structure (initialized)
**************************************************/
void sep_keygen(sep_pk_t *pk, sep_sk_t *sk) {
  uint8_t root_seed[SEED_BYTES], seeds[SEED_BYTES*2];
  uint8_t *public_seed = seeds, *secret_seed = &seeds[SEED_BYTES];
  poly_q_mat_d_d A;
  poly_q_mat_d_d RRstar[2][2];
  uint32_t kappa;

  // generate random seed(s)
  randombytes(root_seed, SEED_BYTES);
  sha3_512(seeds, root_seed, SEED_BYTES);
#if SEED_BYTES != 32
#error "SEED_BYTES must be 32."
#endif

  // init matrices
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]);
  poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]);
  poly_q_mat_d_d_init(RRstar[1][1]);

  // expand uniform A' from seed
  poly_q_mat_d_d_uniform(A, public_seed, DOMAIN_SEPARATOR_A, 0);

  // sample R from B_1 (binomial)
  kappa = 0;
  do {
    for (size_t i = 0; i < PARAM_K; i++) {
      for (size_t j = 0; j < 2; j++) {
        poly_q_mat_d_d_binomial(sk->R[j][i], secret_seed, kappa++, DOMAIN_SEPARATOR_R);
      }
    }
  } while(sk_sq_spectral_norm(RRstar, sk->R) > PARAM_R_MAX_SQ_SPECTRAL_NORM);

  // compute B = (I | A')R
  for (size_t i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_mul_mat_d_d(pk->B[i], A, sk->R[1][i]);
    poly_q_mat_d_d_add(pk->B[i], pk->B[i], sk->R[0][i]);
  }

  // append public seed to pk for extending A', u, D, A3
  for (size_t i = 0; i < SEED_BYTES; i++) {
    pk->seed[i] = public_seed[i];
  }

  // compute S (for perturbation sampling)
  compute_covariance(sk->S, RRstar);

  // clean up
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(RRstar[0][0]);
  poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]);
  poly_q_mat_d_d_clear(RRstar[1][1]);
}

/*************************************************
* Name:        _sep_sign_commitment
*
* Description: Computes the signature from the commitment cmt = D.m (standalone
*              signature) or cmt = (I|A')r + Ds.usk + D.m (anonymous credentials issuance)
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure (initialized)
*              - uint8_t *state: pointer to signer's state byte array (allocated STATE_BYTES bytes)
*              - const sep_sk_t *sk: pointer to secret key structure
*              - const sep_pk_t *pk: pointer to public key structure
*              - const poly_q_vec_d cmt: polynomial vector hosting the commitment to be signed
**************************************************/
void _sep_sign_commitment(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const poly_q_vec_d cmt) {
  size_t i;
  poly_q_mat_d_k A3;
  poly_q_mat_d_d A;
  poly_q_vec_d u, tmp;
  poly_q_vec_d v1[2];
  uint64_t norm2sq_v1, norm2sq_v2, norm2sq_v3;

  // init matrices and vectors
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(tmp);
  poly_q_vec_d_init(v1[0]);
  poly_q_vec_d_init(v1[1]);

  // compute tag from state
  poly_q_binary_fixed_weight(sig->tag, state);

  // increment state
#if (STATE_BYTES%4) != 0
#error "STATE_BYTES must be multiple of 4"
#endif
  uint64_t stateinc = ((uint64_t)*(uint32_t*)state) + 1; // cast byte array to uint32, then to uint64, then increment
  *(uint32_t*)state = (uint32_t) stateinc;
  uint64_t carry = stateinc >> 32;
  for (i = 1; i < STATE_BYTES/4; i++)
  {
    stateinc = ((uint64_t)*(uint32_t*)&state[i*4]) + carry; // cast byte array to uint32, then to uint64, then add carry
    *(uint32_t*)&state[i*4] = (uint32_t) stateinc;
    carry = stateinc >> 32;
  }

  // expand uniform matrices A', A3, and uniform vector u
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

reject_signature:
  // sample v3 from discrete Gaussian
  poly_q_vec_k_sample_gaussian_s2(sig->v3); // probabilistic

  // compute u + cmt - A3.v3
  poly_q_mat_d_k_mul_vec_k(tmp, A3, sig->v3);
  poly_q_vec_d_sub(tmp, cmt, tmp);
  poly_q_vec_d_add(tmp, u, tmp);

  // call SamplePre, output in v1, sig->v2
  poly_q_vec_2d_dk_sample_pre(v1, sig->v2, sk->R, A, pk->B, tmp, sig->tag, sk->S); // probabilistic

  // check l2 norms before outputting signature
  norm2sq_v1 = poly_q_vec_d_norm2(v1[0]);
  norm2sq_v1 += poly_q_vec_d_norm2(v1[1]);
  norm2sq_v2 = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_K; i++) {
    norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
  }
  norm2sq_v3 = poly_q_vec_k_norm2(sig->v3);
  if((norm2sq_v1 > PARAM_B1SQ) || (norm2sq_v2 > PARAM_B2SQ) || (norm2sq_v3 > PARAM_B3SQ)) {
    goto reject_signature;
  }

  // extract v12 from v1
  poly_q_vec_d_set(sig->v12, v1[1]);

  // clean up matrices and vectors
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(tmp);
  poly_q_vec_d_clear(v1[0]);
  poly_q_vec_d_clear(v1[1]);
}

/*************************************************
* Name:        sep_sign
*
* Description: Computes the signature on the message (standalone signature)
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure (initialized)
*              - uint8_t *state: pointer to signer's state byte array (allocated STATE_BYTES bytes)
*              - const sep_sk_t *sk: pointer to secret key structure
*              - const sep_pk_t *pk: pointer to public key structure
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
**************************************************/
void sep_sign(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const uint8_t msg[PARAM_M*PARAM_N/8]) {
  size_t i;
  poly_q_vec_d cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;

  // init matrices and vectors
  poly_q_vec_d_init(cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  // commitment cmt = D.m (standalone signature)
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_m_mul_vec_m(cmt, D, m);

  // sign commitment
  _sep_sign_commitment(sig, state, sk, pk, cmt);

  // clean up matrices and vectors
  poly_q_vec_d_clear(cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);
}

/*************************************************
* Name:        _sep_verify_from_commitment
*
* Description: Verifies signature from the precomputed commitment cmt = D.m (standalone
*              signature) or cmt = Ds.usk + D.m (anonymous credentials issuance)
*
* Arguments:   - const sep_sig_t *sig: pointer to the input signature structure 
*              - const poly_q_vec_d cmt: binding commitment used in the verification
*              - const sep_pk_t *pk: pointer to public key structure
* 
* Returns 1 if signature could be verified correctly and 0 otherwise
**************************************************/
int _sep_verify_from_commitment(const sep_sig_t *sig, const poly_q_vec_d cmt, const sep_pk_t *pk) {
  size_t i, j;
  int64_t bexpi;
  poly_q_vec_d u, v11;
  poly_q_mat_d_k A3;
  poly_q_mat_d_d A;
  poly_q_mat_d_d Btmp;
  poly_q tag_times_bexpi;
  uint64_t norm2sq_v1, norm2sq_v2, norm2sq_v3;
  int64_t tag_weight;

  // init matrices, vectors and polynomials
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Btmp);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(v11);
  poly_q_vec_d_init(u);
  poly_q_init(tag_times_bexpi);

  // expand uniform matrices A', A3, and uniform vector u
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

  // compute v11 = u + cmt - A'.v12 + (B - tG).v2 - A3.v3
  poly_q_vec_d_add(v11, u, cmt);
  // now we can use u as a temp variable
  poly_q_mat_d_k_mul_vec_k(u, A3, sig->v3);
  poly_q_mat_d_d_muladd_vec_d(u, A, sig->v12);
  poly_q_vec_d_sub(v11, v11, u);
  // now we have v11 = u + cmt - A'.v12 - A3.v3
  bexpi = 1;
  for (i = 0; i < PARAM_K; i++) {
    // copy B to Btmp
    poly_q_mat_d_d_set(Btmp, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      if (i == 0) {
        poly_q_sub(Btmp->rows[j]->entries[j], Btmp->rows[j]->entries[j], sig->tag);
      } else {
        poly_q_mul_scalar(tag_times_bexpi, sig->tag, bexpi);
        poly_q_sub(Btmp->rows[j]->entries[j], Btmp->rows[j]->entries[j], tag_times_bexpi);
      }
    }
    poly_q_mat_d_d_muladd_vec_d(v11, Btmp, sig->v2[i]);
    bexpi *= PARAM_B;
  }

  // compute l2 norms of v11||v12, v2, and v3
  norm2sq_v1 = poly_q_vec_d_norm2(v11);
  norm2sq_v1 += poly_q_vec_d_norm2(sig->v12);
  norm2sq_v2 = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_K; i++) {
    norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
  }
  norm2sq_v3 = poly_q_vec_k_norm2(sig->v3);
  tag_weight = poly_q_weight(sig->tag); // returns -1 if polynomial is non-binary, else the number of ones

  // clean up matrices, vectors and polynomials
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(Btmp);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(v11);
  poly_q_vec_d_clear(u);
  poly_q_clear(tag_times_bexpi);

  // return check
  return (norm2sq_v1 <= PARAM_B1SQ) && (norm2sq_v2 <= PARAM_B2SQ) && (norm2sq_v3 <= PARAM_B3SQ) && (tag_weight == PARAM_W);
}

/*************************************************
* Name:        sep_verify
*
* Description: Verifies signature on the message (standalone signature)
*
* Arguments:   - const sep_sig_t *sig: pointer to the input signature structure 
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
*              - const sep_pk_t *pk: pointer to public key structure
* 
* Returns 1 if signature could be verified correctly and 0 otherwise
**************************************************/
int sep_verify(const sep_sig_t *sig, const uint8_t msg[PARAM_M*PARAM_N/8], const sep_pk_t *pk) {
  size_t i;
  int is_valid;
  poly_q_vec_d cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;

  // init matrices and vectors
  poly_q_vec_d_init(cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  // commitment cmt = D.m (standalone signature)
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_m_mul_vec_m(cmt, D, m);

  // verify from commitment
  is_valid = _sep_verify_from_commitment(sig, cmt, pk);

  // clean up matrices and vectors
  poly_q_vec_d_clear(cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);

  return is_valid;
}
