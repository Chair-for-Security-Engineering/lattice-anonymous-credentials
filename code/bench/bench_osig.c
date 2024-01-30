#include "bench_osig.h"

#include "sep.h"
#include "osig.h"
#include "randombytes.h"
#include "random.h"

double osig_user_commit_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t msg[PARAM_M*PARAM_N/8];

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(msg, PARAM_M*PARAM_N/8);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  osig_user_commit(r, cmt, msg, &upk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_signer_sign_commitment_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sep_sig_init(&sig);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_sign_complete_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  osig_user_sig_complete(&sig, r);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_verify_valid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  osig_user_sig_complete(&sig, r);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = osig_user_verify(&sig, &pk, &upk, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
      printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_verify_invalid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
  osig_user_sig_complete(&sig, r);
  msg[0] ^= 1;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = osig_user_verify(&sig, &pk, &upk, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
      printf("FATAL ERROR: benchmarked signature is valid\n");
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_keygen_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  sep_keygen(&pk, &sk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  osig_user_keygen(&upk, &usk, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_embed_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(u[i]);
  }
  for (size_t i= 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);

  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (size_t i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_user_prove_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(u[i]);
  }
  for (size_t i= 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);

  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (size_t i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_signer_verify_valid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(u[i]);
  }
  for (size_t i= 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);

  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
  osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (size_t i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return time;
}

double osig_signer_verify_invalid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  coeff_qiss coeff;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_init(u[i]);
  }
  for (size_t i= 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }
  sep_keygen(&pk, &sk);
  osig_user_keygen(&upk, &usk, pk.seed);
  randombytes(crs_seed, CRS_SEED_BYTES);

  randombytes(msg, PARAM_M*PARAM_N/8);
  osig_user_commit(r, cmt, msg, &upk);
  osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
  osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);
  coeff = poly_qiss_get_coeff(u[0]->entries[0], 0) + 1;
  poly_qiss_set_coeff(u[0]->entries[0], 0, coeff);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked signature is valid\n");
  }
  for (size_t i = 0; i < 2*PARAM_D; i++) {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (size_t i = 0; i < PARAM_M1_K_ISS; i++) {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (size_t i = 0; i < PARAM_D; i++) {
    for (size_t j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (size_t j = 0; j < PARAM_M; j++) {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return time;
}

