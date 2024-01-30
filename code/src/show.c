#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "poly_q_sampling.h"
#include "sep.h"
#include "show.h"
#include "four_squares.h"

/*************************************************
* Name:        show_proof_init
*
* Description: Initialize structure to host the show proof
*              by calling Flint initialization
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
**************************************************/
void show_proof_init(show_proof_t *proof) {
  poly_qshow_vec_m1_init(proof->z1);
  poly_qshow_vec_d_init(proof->tA);
  poly_qshow_vec_256_l_init(proof->tB);
  poly_qshow_vec_l_init(proof->h);
  poly_qshow_init(proof->t1);
  poly_qshow_init(proof->c);
  poly_qshow_vec_m2_init(proof->z2);
}

/*************************************************
* Name:        show_proof_init
*
* Description: Clear structure that hosts the show proof
*              by calling Flint clean up
*
* Arguments:   - show_proof_t *proof: pointer to show proof structure
**************************************************/
void show_proof_clear(show_proof_t *proof) {
  poly_qshow_vec_m1_clear(proof->z1);
  poly_qshow_vec_d_clear(proof->tA);
  poly_qshow_vec_256_l_clear(proof->tB);
  poly_qshow_vec_l_clear(proof->h);
  poly_qshow_clear(proof->t1);
  poly_qshow_clear(proof->c);
  poly_qshow_vec_m2_clear(proof->z2);
}

/*************************************************
* Name:        show_user_embed
*
* Description: Embedding the user relation for the show proof
*              of valid credentials
*
* Arguments:   - poly_qshow_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.A'
*              - poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.B
*              - poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.A3
*              - poly_qshow_mat_k_k *Ds_embed: array of polynomial matrices to host subring embedding of q1.Ds
*              - poly_qshow_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.D
*              - poly_qshow_vec_k *u_embed: array of polynomial vectors to host subring embedding of q1.u
*              - poly_qshow_vec_k *s1: array of polynomial vectors to host subring embedding of witness
*                   s1 = [theta(v1)|a1 | theta(v2)|a2 | theta(v3)|a3 | theta(tag) | theta(usk) | theta(msg)]
*              - const user_pk_t *upk: pointer to user public key structure
*              - const user_sk_t *usk: pointer to user secret key structure
*              - const sep_pk_t *pk: pointer to signer public key structure
*              - const sep_sig_t *sig: pointer to signature structure
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_M*PARAM_N/8 bytes)
**************************************************/
void show_user_embed(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_D*PARAM_K], 
    poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    poly_qshow_mat_k_k Ds_embed[PARAM_D][2*PARAM_D], 
    poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], 
    poly_qshow_vec_k   u_embed[PARAM_D], 
    poly_qshow_vec_m1  s1, 
    const user_pk_t    *upk, 
    const user_sk_t    *usk, 
    const sep_pk_t     *pk, 
    const sep_sig_t    *sig,
    const uint8_t      msg[PARAM_M*PARAM_N/8]) {
  size_t i,j;
  int64_t bexpi;
  uint64_t norm2sq_v1, norm2sq_v2, norm2sq_v3, four_squares_res[4];
  poly_q tmp_poly;
  poly_q_mat_d_d A, Ds[2], Btmp;
  poly_q_mat_d_k A3;
  poly_q_mat_d_m D;
  poly_q_vec_d u, v11;
  poly_q_vec_m m;
  poly_qshow_vec_k tmp_vec_k;

  // init matrices, vectors and polynomials
  poly_q_init(tmp_poly);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(Btmp);
  poly_q_mat_d_k_init(A3);
  poly_q_mat_d_m_init(D);
  poly_q_mat_d_d_init(Ds[0]);
  poly_q_mat_d_d_init(Ds[1]);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(v11);
  poly_q_vec_m_init(m);
  poly_qshow_vec_k_init(tmp_vec_k);

  // expanding uniform public parameters
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_mat_d_d_uniform(Ds[0], pk->seed, DOMAIN_SEPARATOR_DS, 0);
  poly_q_mat_d_d_uniform(Ds[1], pk->seed, DOMAIN_SEPARATOR_DS, PARAM_D);
  poly_q_mat_d_m_uniform(D, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

  // embedding u (u_embed = theta(q1.u))
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(u_embed[i], u->entries[i], PARAM_Q1_SHOW);
  }

  // recovering v_11 = u + Ds.usk + D.m - A'.v12 - (tG - B).v2 - A3.v3
  // storing message as polynomial vector
  for (i = 0; i < PARAM_M; i++) {
    poly_q_from_bits(m->entries[i], &msg[i * PARAM_N/8]);
  }
  poly_q_mat_d_m_mul_vec_m(v11, D, m);
  poly_q_vec_d_add(v11, v11, u); // u can be used as a temp variable now
  poly_q_vec_d_add(v11, v11, upk->t); // Ds.usk = upk
  poly_q_mat_d_k_mul_vec_k(u, A3, sig->v3);
  poly_q_mat_d_d_muladd_vec_d(u, A, sig->v12);
  poly_q_vec_d_sub(v11, v11, u);

  // now we have v11 = u + Ds.usk + D.m - A.v12 - A3.v3
  bexpi = 1;
  for (i = 0; i < PARAM_K; i++) {
    // copy B to Btmp
    poly_q_mat_d_d_set(Btmp, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      if (i == 0) {
        poly_q_sub(Btmp->rows[j]->entries[j], Btmp->rows[j]->entries[j], sig->tag);
      } else {
        poly_q_mul_scalar(tmp_poly, sig->tag, bexpi);
        poly_q_sub(Btmp->rows[j]->entries[j], Btmp->rows[j]->entries[j], tmp_poly);
      }
    }
    poly_q_mat_d_d_muladd_vec_d(v11, Btmp, sig->v2[i]);
    bexpi *= PARAM_B;
  }

  // computing square norms and four-square decompositions
  // compute l2 norms of v11||v12, v2, and v3
  norm2sq_v1 = poly_q_vec_d_norm2(v11);
  norm2sq_v1 += poly_q_vec_d_norm2(sig->v12);
  norm2sq_v2 = poly_q_vec_d_norm2(sig->v2[0]);
  for (i = 1; i < PARAM_K; i++) {
    norm2sq_v2 += poly_q_vec_d_norm2(sig->v2[i]);
  }
  norm2sq_v3 = poly_q_vec_k_norm2(sig->v3);
  // B1^2 - |v1|^2
  assert(norm2sq_v1 <= PARAM_B1SQ);
  four_squares(four_squares_res, PARAM_B1SQ - norm2sq_v1);
  for (i = 0; i < 4; i++) {
    poly_qshow_set_coeff(s1->entries[IDX_V2_SHOW - 1], i, (coeff_qshow)(four_squares_res[i]));
  }
  // B2^2 - |v2|^2
  assert(norm2sq_v2 <= PARAM_B2SQ);
  four_squares(four_squares_res, PARAM_B2SQ - norm2sq_v2);
  for (i = 0; i < 4; i++) {
    poly_qshow_set_coeff(s1->entries[IDX_V3_SHOW - 1], i, (coeff_qshow)(four_squares_res[i]));
  }
  // B3^2 - |v3|^2
  assert(norm2sq_v3 <= PARAM_B3SQ);
  four_squares(four_squares_res, PARAM_B3SQ - norm2sq_v3);
  for (i = 0; i < 4; i++) {
    poly_qshow_set_coeff(s1->entries[IDX_TAG_SHOW - 1], i, (coeff_qshow)(four_squares_res[i]));
  }

  // embedding witness vector s1 = [theta(v1)|a1 | theta(v2)|a2 | theta(v3)|a3 | theta(tag) | theta(usk) | theta(msg)]
  // v1
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, v11->entries[i], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_V1_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
    poly_qshow_subring_embed_vec_k(tmp_vec_k, sig->v12->entries[i], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_V1_SHOW + (i + PARAM_D)*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // v2
  for (i = 0; i < PARAM_D*PARAM_K; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, sig->v2[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_V2_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // v3
  for (i = 0; i < PARAM_K; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, sig->v3->entries[i], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_V3_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // tag
  poly_qshow_subring_embed_vec_k(tmp_vec_k, sig->tag, 1);
  for (j = 0; j < PARAM_K_SHOW; j++) {
    poly_qshow_set(s1->entries[IDX_TAG_SHOW + j], tmp_vec_k->entries[j]);
  }
  // usk
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, usk->s[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_USK_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // msg
  for (i = 0; i < PARAM_M; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, m->entries[i], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_M_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }

  // embedding A, B, A3, D, Ds
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_SHOW); // A'
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j          ], Ds[0]->rows[i]->entries[j], PARAM_Q1_SHOW); // Ds[:,0:PARAM_D]
      poly_qshow_subring_embed_mat_k_k(Ds_embed[i][j + PARAM_D], Ds[1]->rows[i]->entries[j], PARAM_Q1_SHOW); // Ds[:,PARAM_D:]
    }
    for (j = 0; j < PARAM_D*PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(B_embed[i][j], pk->B[j / PARAM_D]->rows[i]->entries[j % PARAM_D], PARAM_Q1_SHOW); // B
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(A3_embed[i][j], A3->rows[i]->entries[j], PARAM_Q1_SHOW); // A3
    }
    for (j = 0; j < PARAM_M; j++) {
      poly_qshow_subring_embed_mat_k_k(D_embed[i][j], D->rows[i]->entries[j], PARAM_Q1_SHOW); // D
    }
  }

  // clean up matrices, vectors and polynomials
  poly_q_clear(tmp_poly);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(Btmp);
  poly_q_mat_d_k_clear(A3);
  poly_q_mat_d_m_clear(D);
  poly_q_mat_d_d_clear(Ds[0]);
  poly_q_mat_d_d_clear(Ds[1]);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(v11);
  poly_q_vec_m_clear(m);
  poly_qshow_vec_k_clear(tmp_vec_k);
}
