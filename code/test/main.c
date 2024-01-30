#include <stdio.h>
#include "sep.h"
#include "osig.h"
#include "show.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 5

static int sep_test(void)
{
  int rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nsep_test\n");

  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      sep_sign(&sig, state, &sk, &pk, msg);
      if (!sep_verify(&sig, msg, &pk))
      {
        printf("sep_verify returned zero for a valid signature.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      msg[0] ^= 1;
      if (sep_verify(&sig, msg, &pk))
      {
        printf("sep_verify returned non-zero for a valid signature on the wrong message.\n");
        rval = 0;
        goto sep_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

  sep_test_cleanup:
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int osig_signing_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  sep_sig_t sig;
  user_sk_t usk;
  user_pk_t upk;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nosig_signing_test\n");

  // init
  sep_keys_init(&pk, &sk);
  sep_sig_init(&sig);
  user_keys_init(&upk, &usk);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
      osig_user_sig_complete(&sig, r);
      if (!osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned zero for a valid signature.\n");
        rval = 0;
        goto osig_signing_test_cleanup;
      }
      msg[0] ^= 1;
      if (osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned non-zero for a valid signature on the wrong message.\n");
        rval = 0;
        goto osig_signing_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

osig_signing_test_cleanup:
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int osig_proof_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  osig_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  coeff_qiss coeff;
  poly_qiss_vec_k u[2*PARAM_D], s1[PARAM_M1_K_ISS];
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  printf("\nosig_proof_test\n");

  // init
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  osig_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  for (i = 0; i < 2*PARAM_D; i++)
  {
    poly_qiss_vec_k_init(u[i]);
  }
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_init(s1[i]);
  }
  for (i = 0; i < PARAM_D; i++)
  {
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qiss_mat_k_k_init(D_embed[i][j]);
    }
  }

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    randombytes(crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_user_embed(A_embed, Ds_embed, D_embed, u, s1, &upk, &usk, cmt, r, msg);
      osig_user_prove(&proof, A_embed, Ds_embed, D_embed, u, s1, crs_seed, upk.seed);
      if (!osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed))
      {
        printf("osig_signer_verify returned zero for a valid proof.\n");
        rval = 0;
        goto osig_proof_test_cleanup;
      }
      coeff = poly_qiss_get_coeff(u[0]->entries[0], 0) + 1;
      poly_qiss_set_coeff(u[0]->entries[0], 0, coeff);
      if (osig_signer_verify(&proof, A_embed, Ds_embed, D_embed, u, crs_seed, upk.seed))
      {
        printf("osig_signer_verify returned non-zero for a valid proof but tampered commitment.\n");
        rval = 0;
        goto osig_proof_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

osig_proof_test_cleanup:
  for (i = 0; i < 2*PARAM_D; i++)
  {
    poly_qiss_vec_k_clear(u[i]);
  }
  for (i = 0; i < PARAM_M1_K_ISS; i++)
  {
    poly_qiss_vec_k_clear(s1[i]);
  }
  for (i = 0; i < PARAM_D; i++)
  {
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qiss_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qiss_mat_k_k_clear(D_embed[i][j]);
    }
  }
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(r[0]);
  osig_proof_clear(&proof);
  user_keys_clear(&upk, &usk);
  sep_keys_clear(&pk, &sk);
  return rval;
}

static int show_proof_test(void)
{
  int i,j,rval = 1;
  sep_sk_t sk;
  sep_pk_t pk;
  user_sk_t usk;
  user_pk_t upk;
  sep_sig_t sig;
  show_proof_t proof;
  poly_q_vec_d r[2];
  poly_q_vec_d cmt;
  coeff_qshow coeff;
  poly_qshow_vec_m1 s1;
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K], A3_embed[PARAM_D][PARAM_K];
  poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
  uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  randombytes(state, STATE_BYTES);

  printf("\nshow_proof_test\n");

  // init
  sep_keys_init(&pk, &sk);
  user_keys_init(&upk, &usk);
  sep_sig_init(&sig);
  show_proof_init(&proof);
  poly_q_vec_d_init(r[0]);
  poly_q_vec_d_init(r[1]);
  poly_q_vec_d_init(cmt);
  poly_qshow_vec_m1_init(s1);
  for (i = 0; i < PARAM_D; i++)
  {
    poly_qshow_vec_k_init(u_embed[i]);
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
      poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D*PARAM_K; j++)
    {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++)
    {
      poly_qshow_mat_k_k_init(A3_embed[i][j]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qshow_mat_k_k_init(D_embed[i][j]);
    }
  }

  for (i = 0; i < NSUBTESTS; i++)
  {
    sep_keygen(&pk, &sk);
    osig_user_keygen(&upk, &usk, pk.seed);
    randombytes(crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_M*PARAM_N/8);
      osig_user_commit(r, cmt, msg, &upk);
      osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
      osig_user_sig_complete(&sig, r);
      if (!osig_user_verify(&sig, &pk, &upk, msg))
      {
        printf("osig_user_verify returned zero for a valid signature.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
      show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, s1, crs_seed, upk.seed);
      if (!show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed))
      {
        printf("show_verify returned zero for a valid proof.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      coeff = poly_qshow_get_coeff(u_embed[0]->entries[0], 0) + 1;
      poly_qshow_set_coeff(u_embed[0]->entries[0], 0, coeff);
      if (show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed))
      {
        printf("show_verify returned non-zero for a valid proof but wrong statement.\n");
        rval = 0;
        goto show_proof_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

show_proof_test_cleanup:
  sep_keys_clear(&pk, &sk);
  user_keys_clear(&upk, &usk);
  sep_sig_clear(&sig);
  show_proof_clear(&proof);
  poly_q_vec_d_clear(r[0]);
  poly_q_vec_d_clear(r[1]);
  poly_q_vec_d_clear(cmt);
  poly_qshow_vec_m1_clear(s1);
  for (i = 0; i < PARAM_D; i++)
  {
    poly_qshow_vec_k_clear(u_embed[i]);
    for (j = 0; j < PARAM_D; j++)
    {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
      poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
    }
    for (j = 0; j < PARAM_D*PARAM_K; j++)
    {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++)
    {
      poly_qshow_mat_k_k_clear(A3_embed[i][j]);
    }
    for (j = 0; j < PARAM_M; j++)
    {
      poly_qshow_mat_k_k_clear(D_embed[i][j]);
    }
  }
  return rval;
}

int main(void) {
  int pass = 1;
  arith_setup();
  random_init();
  printf("Hello from the unit tests.\n");
  for (int i = 0; i < NTESTS; i++)
  {
    pass &= show_proof_test();
    pass &= osig_proof_test();
    pass &= osig_signing_test();
    pass &= sep_test();
    // TODO add more tests

    if (!pass)
    {
      printf("FAILED!\n");
      break;
    } else {
      printf(".");
    }
  }
  if (pass)
  {
    printf("passed.\n");
  }
  arith_teardown();
  return 0;
}
