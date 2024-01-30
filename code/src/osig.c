#include "arith.h"
#include "randombytes.h"
#include "poly_qiss_sampling.h"
#include "poly_q_sampling.h"
#include "sep.h"
#include "osig.h"

/*************************************************
* Name:        user_keys_init
*
* Description: Initialize user keys for the anonymous credentials
*              by calling Flint initialization
*
* Arguments:   - user_pk_t *upk: pointer to user public key structure
*              - user_sk_t *usk: pointer to user secret key structure
**************************************************/
void user_keys_init(user_pk_t *upk, user_sk_t *usk) {
  poly_q_vec_d_init(usk->s[0]);
  poly_q_vec_d_init(usk->s[1]);
  poly_q_vec_d_init(upk->t);
}

/*************************************************
* Name:        user_keys_clear
*
* Description: Clear user keys for the anonymous credentials
*              by calling Flint clean up
*
* Arguments:   - user_pk_t *upk: pointer to user public key structure
*              - user_sk_t *usk: pointer to user secret key structure
**************************************************/
void user_keys_clear(user_pk_t *upk, user_sk_t *usk) {
  poly_q_vec_d_clear(usk->s[0]);
  poly_q_vec_d_clear(usk->s[1]);
  poly_q_vec_d_clear(upk->t);
}

/*************************************************
* Name:        osig_proof_init
*
* Description: Initialize structure to host the issuance proof
*              by calling Flint initialization
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
**************************************************/
void osig_proof_init(osig_proof_t *proof) {
  size_t i;
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_init(proof->z1[i]);
  }
  poly_qiss_vec_d_init(proof->tA);
  poly_qiss_vec_256_l_init(proof->tB);
  poly_qiss_vec_l_init(proof->h);
  poly_qiss_init(proof->t1);
  poly_qiss_init(proof->c);
  poly_qiss_vec_m2_init(proof->z2);
}

/*************************************************
* Name:        osig_proof_clear
*
* Description: Clear structure that hosts the issuance proof
*              by calling Flint clean up
*
* Arguments:   - osig_proof_t *proof: pointer to issuance proof structure
**************************************************/
void osig_proof_clear(osig_proof_t *proof) {
  size_t i;
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_clear(proof->z1[i]);
  }
  poly_qiss_vec_d_clear(proof->tA);
  poly_qiss_vec_256_l_clear(proof->tB);
  poly_qiss_vec_l_clear(proof->h);
  poly_qiss_clear(proof->t1);
  poly_qiss_clear(proof->c);
  poly_qiss_vec_m2_clear(proof->z2);
}

/*************************************************
* Name:        osig_user_keygen
*
* Description: Generates user public and private key
*
* Arguments:   - user_pk_t *upk: pointer to user public key structure (initialized) (contains seed)
*              - user_sk_t *usk: pointer to user secret key structure (initialized)
*              - const uint8_t *seed: pointer to byte array containing the seed 
*                   for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void osig_user_keygen(user_pk_t *upk, user_sk_t *usk, const uint8_t seed[SEED_BYTES]) {
  size_t i;
  uint8_t secret_seed[SEED_BYTES];
  poly_q_mat_d_d Ds[2];
  poly_q_vec_d tmp;

  // init matrices and vectors
  poly_q_mat_d_d_init(Ds[0]);
  poly_q_mat_d_d_init(Ds[1]);
  poly_q_vec_d_init(tmp);

  // generate random secret seed
  randombytes(secret_seed, SEED_BYTES);

  // expand uniform Ds from seed
  poly_q_mat_d_d_uniform(Ds[0], seed, DOMAIN_SEPARATOR_DS, 0);
  poly_q_mat_d_d_uniform(Ds[1], seed, DOMAIN_SEPARATOR_DS, PARAM_D);

  // sample s from U({0,1})
  // TODO usk->s could be just an uint8_t buffer and this uniform distribution could be generated here, just with shake256
  poly_q_vec_d_bin_uniform(usk->s[0], secret_seed, DOMAIN_SEPARATOR_S, 0);
  poly_q_vec_d_bin_uniform(usk->s[1], secret_seed, DOMAIN_SEPARATOR_S, PARAM_D);

  // compute t = Ds.s
  poly_q_mat_d_d_mul_vec_d(upk->t, Ds[0], usk->s[0]);
  poly_q_mat_d_d_mul_vec_d(tmp, Ds[1], usk->s[1]);
  poly_q_vec_d_add(upk->t, upk->t, tmp);

  /**********************************************
  * In the single signer setting, the seed for A', Ds, D 
  * is that of the signer. In the multiple signers setting, 
  * the seed is that of the public parameters shared across 
  * many signers.
  **********************************************/
  // adding seed to upk to derive A', Ds, D
  for (i = 0; i < SEED_BYTES; i++) {
    upk->seed[i] = seed[i];
  }

  // clean up matrices and vectors
  poly_q_mat_d_d_clear(Ds[0]);
  poly_q_mat_d_d_clear(Ds[1]);
  poly_q_vec_d_clear(tmp);
}

/*************************************************
* Name:        osig_user_commit
*
* Description: Compute hiding commitment to (usk | msg) by cmt = (I|A')r + Ds.usk + D.msg
*
* Arguments:   - poly_q_vec_d *r: array of polynomial vectors for user commitment randomness (initialized)
*              - poly_q_vec_d cmt: polynomial vector for user commitment (initialized)
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
*              - user_pk_t *upk: pointer to user public key structure
**************************************************/
void osig_user_commit(poly_q_vec_d r[2], poly_q_vec_d cmt, const uint8_t msg[PARAM_M*PARAM_N/8], const user_pk_t *upk) {
  size_t i;
  uint8_t randomness_seed[SEED_BYTES];
  poly_q_mat_d_d A;
  poly_q_mat_d_m D;
  poly_q_vec_m m;
  poly_q_vec_d tmp;

  // init matrices and vectors
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);
  poly_q_vec_d_init(tmp);

  // generate random secret seed for commitment randomness
  randombytes(randomness_seed, SEED_BYTES);

  // expand uniform A', D from seed
  poly_q_mat_d_d_uniform(A, upk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_m_uniform(D, upk->seed, DOMAIN_SEPARATOR_D);

  // storing message as polynomial vector
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }

  // sample r from U({0,1})
  poly_q_vec_d_bin_uniform(r[0], randomness_seed, DOMAIN_SEPARATOR_RAND, 0);
  poly_q_vec_d_bin_uniform(r[1], randomness_seed, DOMAIN_SEPARATOR_RAND, PARAM_D);

  // cmt = r[0] + A'.r[1] + Ds.usk + D.msg (but Ds.usk is already computed in upk)
  poly_q_mat_d_d_mul_vec_d(cmt, A, r[1]);
  poly_q_vec_d_add(cmt, cmt, r[0]);
  poly_q_vec_d_add(cmt, cmt, upk->t);
  poly_q_mat_d_m_mul_vec_m(tmp, D, m);
  poly_q_vec_d_add(cmt, cmt, tmp);

  // clean up matrices and vectors
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);
  poly_q_vec_d_clear(tmp);
}

/*************************************************
* Name:        osig_user_embed
*
* Description: Embedding the user relation for the issuance proof
*              of commitment opening and user registration
*
* Arguments:   - poly_qiss_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.A'
*              - poly_qiss_mat_k_k *Ds_embed: array of polynomial matrices to host subring embedding of q1.Ds
*              - poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.D
*              - poly_qiss_vec_k *u: array of polynomial vectors to host subring embedding of q1.(cmt-upk | upk)
*              - poly_qiss_vec_k *s1: array of polynomial vectors to host subring embedding of (r|usk|msg)
*              - const user_pk_t *upk: pointer to user public key structure
*              - const user_sk_t *usk: pointer to user secret key structure
*              - const poly_q_vec_d cmt: polynomial vector for commitment
*              - const poly_q_vec_d *r: array of polynomial vectors for commitment randomness
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
**************************************************/
void osig_user_embed(
    poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qiss_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    poly_qiss_mat_k_k D_embed[PARAM_D][PARAM_M], 
    poly_qiss_vec_k u[2*PARAM_D], 
    poly_qiss_vec_k s1[PARAM_M1_K_ISS], 
    const user_pk_t *upk, 
    const user_sk_t *usk, 
    const poly_q_vec_d cmt, 
    const poly_q_vec_d r[2], 
    const uint8_t msg[PARAM_M*PARAM_N/8]) {
  size_t i,j;
  poly_q tmp;
  poly_q_mat_d_d A;
  poly_q_mat_d_d Ds[2];
  poly_q_mat_d_m D;

  // init matrices and polynomials
  poly_q_init(tmp);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Ds[0]);
  poly_q_mat_d_d_init(Ds[1]);
  poly_q_mat_d_m_init(D);

  // embedding witness vector s1 = [theta(r) | theta(s) | theta(m)]
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_subring_embed_vec_k(s1[i            ], r[i/PARAM_D]->entries[i%PARAM_D], 1); // r
    poly_qiss_subring_embed_vec_k(s1[i + 2*PARAM_D], usk->s[i/PARAM_D]->entries[i%PARAM_D], 1); // s
  }
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(tmp, &msg[i * PARAM_N/8]);
    poly_qiss_subring_embed_vec_k(s1[i + 4*PARAM_D], tmp, 1); // msg
  }

  // embedding syndrome u = q1 * [theta(cmt-upk) | theta(upk)]
  for (i = 0; i < PARAM_D; i++) {
    poly_q_sub(tmp, cmt->entries[i], upk->t->entries[i]);
    poly_qiss_subring_embed_vec_k(u[i          ], tmp, PARAM_Q1_ISS);                // cmt - upk 
    poly_qiss_subring_embed_vec_k(u[i + PARAM_D], upk->t->entries[i], PARAM_Q1_ISS); // upk
  }

  // expanding uniform A', D, Ds
  poly_q_mat_d_d_uniform(A, upk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_m_uniform(D, upk->seed, DOMAIN_SEPARATOR_D);
  poly_q_mat_d_d_uniform(Ds[0], upk->seed, DOMAIN_SEPARATOR_DS, 0);
  poly_q_mat_d_d_uniform(Ds[1], upk->seed, DOMAIN_SEPARATOR_DS, PARAM_D);

  // embedding A, D, Ds
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_ISS); // A'
      poly_qiss_subring_embed_mat_k_k(Ds_embed[i][j          ], Ds[0]->rows[i]->entries[j], PARAM_Q1_ISS); // Ds[:,0:PARAM_D]
      poly_qiss_subring_embed_mat_k_k(Ds_embed[i][j + PARAM_D], Ds[1]->rows[i]->entries[j], PARAM_Q1_ISS); // Ds[:,PARAM_D:]
    }
    for (j = 0; j < PARAM_M; j++) {
      poly_qiss_subring_embed_mat_k_k(D_embed[i][j], D->rows[i]->entries[j], PARAM_Q1_ISS); // D
    }
  }

  // clean up matrices and polynomials
  poly_q_clear(tmp);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(Ds[0]);
  poly_q_mat_d_d_clear(Ds[1]);
  poly_q_mat_d_m_clear(D);
}

/*************************************************
* Name:        osig_signer_sign_commitment
*
* Description: Computes the signature from the commitment cmt = (I|A')r + Ds.usk + D.m 
*              (anonymous credentials issuance)
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure (initialized)
*              - uint8_t *state: pointer to signer's state byte array (allocated STATE_BYTES bytes)
*              - const sep_sk_t *sk: pointer to secret key structure
*              - const sep_pk_t *pk: pointer to public key structure
*              - const poly_q_vec_d cmt: polynomial vector hosting the commitment to be signed
**************************************************/
void osig_signer_sign_commitment(sep_sig_t *sig, uint8_t state[STATE_BYTES], const sep_sk_t *sk, const sep_pk_t *pk, const poly_q_vec_d cmt) {
  _sep_sign_commitment(sig, state, sk, pk, cmt);
}

/*************************************************
* Name:        osig_user_sig_complete
*
* Description: Merge user commitment randomness to obtain signature (anonymous credentials issuance)
*
* Arguments:   - sep_sig_t *sig: pointer to signature structure (initialized)
*              - const poly_q_vec_d *r: array of polynomial vectors for commitment randomness
**************************************************/
void osig_user_sig_complete(sep_sig_t *sig, const poly_q_vec_d r[2]) {
  poly_q_vec_d_sub(sig->v12, sig->v12, r[1]);
}

/*************************************************
* Name:        osig_user_verify
*
* Description: User verification of the obtained signature (anonymous credentials issuance)
*
* Arguments:   - const sep_sig_t *sig: pointer to the input signature structure 
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
*              - const sep_pk_t *pk: pointer to public key structure
* 
* Returns 1 if signature could be verified correctly and 0 otherwise
**************************************************/
int osig_user_verify(const sep_sig_t *sig, const sep_pk_t *pk, const user_pk_t *upk, const uint8_t msg[PARAM_M*PARAM_N/8]) {
  size_t i;
  int is_valid;
  poly_q_vec_d bind_cmt;
  poly_q_mat_d_m D;
  poly_q_vec_m m;

  // init matrices and vectors
  poly_q_vec_d_init(bind_cmt);
  poly_q_mat_d_m_init(D);
  poly_q_vec_m_init(m);

  // storing message as polynomial vector
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }

  // expand uniform D from seed
  poly_q_mat_d_m_uniform(D, upk->seed, DOMAIN_SEPARATOR_D);

  // computing binding commitment bind_cmt = Ds.usk + D.msg (randomness already included in sig)
  poly_q_mat_d_m_mul_vec_m(bind_cmt, D, m);
  poly_q_vec_d_add(bind_cmt, bind_cmt, upk->t);

  // Signature verification
  is_valid = _sep_verify_from_commitment(sig, bind_cmt, pk);

  // clean up matrices and vectors
  poly_q_vec_d_clear(bind_cmt);
  poly_q_mat_d_m_clear(D);
  poly_q_vec_m_clear(m);

  return is_valid;
}
