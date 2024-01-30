#include "bench_sep.h"

#include "sep.h"
#include "osig.h"
#include "randombytes.h"
#include "random.h"

double sep_keygen_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sep_keys_init(&pk, &sk);
  sep_keygen(&pk, &sk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  sep_keys_clear(&pk, &sk);
  return time;
}

double sep_sign_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  sep_keys_init(&pk, &sk);
  sep_keygen(&pk, &sk);
  randombytes(msg, PARAM_M*PARAM_N/8);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sep_sig_init(&sig);
  sep_sign(&sig, state, &sk, &pk, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!sep_verify(&sig, msg, &pk)) {
    printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double sep_verify_valid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  sep_keys_init(&pk, &sk);
  sep_keygen(&pk, &sk);
  randombytes(msg, PARAM_M*PARAM_N/8);
  sep_sig_init(&sig);
  sep_sign(&sig, state, &sk, &pk, msg);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = sep_verify(&sig, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
      printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

double sep_verify_invalid_bench(timer* t) {
  double time;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  sep_keys_init(&pk, &sk);
  sep_keygen(&pk, &sk);
  randombytes(msg, PARAM_M*PARAM_N/8);
  sep_sig_init(&sig);
  sep_sign(&sig, state, &sk, &pk, msg);
  msg[0] ^= 1;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = sep_verify(&sig, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
      printf("FATAL ERROR: benchmarked signature is valid\n");
  }
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return time;
}

