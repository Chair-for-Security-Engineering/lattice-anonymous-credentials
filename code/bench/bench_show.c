#include "bench_show.h"

#include "sep.h"
#include "osig.h"
#include "show.h"
#include "randombytes.h"
#include "random.h"

double show_user_embed_bench(timer* t) {
	double time;
	int i,j;
	sep_sk_t sk;
	sep_pk_t pk;
	user_sk_t usk;
	user_pk_t upk;
	sep_sig_t sig;
	poly_q_vec_d r[2];
	poly_q_vec_d cmt;
	poly_qshow_vec_m1 s1;
	poly_qshow_vec_k u_embed[PARAM_D];
	poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K], A3_embed[PARAM_D][PARAM_K];
	poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
	uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
	randombytes(state, STATE_BYTES);

	sep_keys_init(&pk, &sk);
	user_keys_init(&upk, &usk);
	sep_sig_init(&sig);
	poly_q_vec_d_init(r[0]);
	poly_q_vec_d_init(r[1]);
	poly_q_vec_d_init(cmt);
	poly_qshow_vec_m1_init(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_init(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_init(A_embed[i][j]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_init(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_init(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_init(D_embed[i][j]);
		}
	}
	sep_keygen(&pk, &sk);
	osig_user_keygen(&upk, &usk, pk.seed);
	randombytes(crs_seed, CRS_SEED_BYTES);
	randombytes(msg, PARAM_M*PARAM_N/8);
	osig_user_commit(r, cmt, msg, &upk);
	osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
	osig_user_sig_complete(&sig, r);
	/* ----------- BEGIN: Code under measurement --------- */
	start_timer(t);
	show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
	time = stop_timer(t);
	/* ----------- END: Code under measurement ----------- */
	sep_keys_clear(&pk, &sk);
	user_keys_clear(&upk, &usk);
	sep_sig_clear(&sig);
	poly_q_vec_d_clear(r[0]);
	poly_q_vec_d_clear(r[1]);
	poly_q_vec_d_clear(cmt);
	poly_qshow_vec_m1_clear(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_clear(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_clear(A_embed[i][j]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_clear(D_embed[i][j]);
		}
	}
	return time;
}

double show_user_prove_bench(timer* t) {
	double time;
	int i,j;
	sep_sk_t sk;
	sep_pk_t pk;
	user_sk_t usk;
	user_pk_t upk;
	sep_sig_t sig;
	show_proof_t proof;
	poly_q_vec_d r[2];
	poly_q_vec_d cmt;
	poly_qshow_vec_m1 s1;
	poly_qshow_vec_k u_embed[PARAM_D];
	poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K], A3_embed[PARAM_D][PARAM_K];
	poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
	uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
	randombytes(state, STATE_BYTES);

	sep_keys_init(&pk, &sk);
	user_keys_init(&upk, &usk);
	sep_sig_init(&sig);
	show_proof_init(&proof);
	poly_q_vec_d_init(r[0]);
	poly_q_vec_d_init(r[1]);
	poly_q_vec_d_init(cmt);
	poly_qshow_vec_m1_init(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_init(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_init(A_embed[i][j]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_init(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_init(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_init(D_embed[i][j]);
		}
	}
	sep_keygen(&pk, &sk);
	osig_user_keygen(&upk, &usk, pk.seed);
	randombytes(crs_seed, CRS_SEED_BYTES);
	randombytes(msg, PARAM_M*PARAM_N/8);
	osig_user_commit(r, cmt, msg, &upk);
	osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
	osig_user_sig_complete(&sig, r);
	show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
	/* ----------- BEGIN: Code under measurement --------- */
	start_timer(t);
	show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, s1, crs_seed, upk.seed);
	time = stop_timer(t);
	/* ----------- END: Code under measurement ----------- */
	sep_keys_clear(&pk, &sk);
	user_keys_clear(&upk, &usk);
	sep_sig_clear(&sig);
	show_proof_clear(&proof);
	poly_q_vec_d_clear(r[0]);
	poly_q_vec_d_clear(r[1]);
	poly_q_vec_d_clear(cmt);
	poly_qshow_vec_m1_clear(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_clear(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_clear(A_embed[i][j]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_clear(D_embed[i][j]);
		}
	}
	return time;
}

double show_user_verify_valid_bench(timer* t) {
	double time;
	int i,j;
	sep_sk_t sk;
	sep_pk_t pk;
	user_sk_t usk;
	user_pk_t upk;
	sep_sig_t sig;
	show_proof_t proof;
	poly_q_vec_d r[2];
	poly_q_vec_d cmt;
	poly_qshow_vec_m1 s1;
	poly_qshow_vec_k u_embed[PARAM_D];
	poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_D*PARAM_K], A3_embed[PARAM_D][PARAM_K];
	poly_qshow_mat_k_k D_embed[PARAM_D][PARAM_M], Ds_embed[PARAM_D][2*PARAM_D];
	uint8_t state[STATE_BYTES], msg[PARAM_M*PARAM_N/8], crs_seed[CRS_SEED_BYTES];
	randombytes(state, STATE_BYTES);

	sep_keys_init(&pk, &sk);
	user_keys_init(&upk, &usk);
	sep_sig_init(&sig);
	show_proof_init(&proof);
	poly_q_vec_d_init(r[0]);
	poly_q_vec_d_init(r[1]);
	poly_q_vec_d_init(cmt);
	poly_qshow_vec_m1_init(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_init(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_init(A_embed[i][j]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_init(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_init(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_init(D_embed[i][j]);
		}
	}
	sep_keygen(&pk, &sk);
	osig_user_keygen(&upk, &usk, pk.seed);
	randombytes(crs_seed, CRS_SEED_BYTES);
	randombytes(msg, PARAM_M*PARAM_N/8);
	osig_user_commit(r, cmt, msg, &upk);
	osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
	osig_user_sig_complete(&sig, r);
	show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
	show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, s1, crs_seed, upk.seed);
	/* ----------- BEGIN: Code under measurement --------- */
	start_timer(t);
	int is_valid = show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed);
	time = stop_timer(t);
	/* ----------- END: Code under measurement ----------- */
	if (!is_valid) {
		printf("FATAL ERROR: benchmarked proof is not valid\n");
	}
	sep_keys_clear(&pk, &sk);
	user_keys_clear(&upk, &usk);
	sep_sig_clear(&sig);
	show_proof_clear(&proof);
	poly_q_vec_d_clear(r[0]);
	poly_q_vec_d_clear(r[1]);
	poly_q_vec_d_clear(cmt);
	poly_qshow_vec_m1_clear(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_clear(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_clear(A_embed[i][j]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_clear(D_embed[i][j]);
		}
	}
	return time;
}

double show_user_verify_invalid_bench(timer* t) {
	double time;
	int i,j;
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

	sep_keys_init(&pk, &sk);
	user_keys_init(&upk, &usk);
	sep_sig_init(&sig);
	show_proof_init(&proof);
	poly_q_vec_d_init(r[0]);
	poly_q_vec_d_init(r[1]);
	poly_q_vec_d_init(cmt);
	poly_qshow_vec_m1_init(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_init(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_init(A_embed[i][j]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_init(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_init(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_init(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_init(D_embed[i][j]);
		}
	}
	sep_keygen(&pk, &sk);
	osig_user_keygen(&upk, &usk, pk.seed);
	randombytes(crs_seed, CRS_SEED_BYTES);
	randombytes(msg, PARAM_M*PARAM_N/8);
	osig_user_commit(r, cmt, msg, &upk);
	osig_signer_sign_commitment(&sig, state, &sk, &pk, cmt);
	osig_user_sig_complete(&sig, r);
	show_user_embed(A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, s1, &upk, &usk, &pk, &sig, msg);
	show_user_prove(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, s1, crs_seed, upk.seed);
	coeff = poly_qshow_get_coeff(u_embed[0]->entries[0], 0) + 1;
	poly_qshow_set_coeff(u_embed[0]->entries[0], 0, coeff);
	/* ----------- BEGIN: Code under measurement --------- */
	start_timer(t);
	int is_valid = show_verify(&proof, A_embed, B_embed, A3_embed, Ds_embed, D_embed, u_embed, crs_seed, upk.seed);
	time = stop_timer(t);
	/* ----------- END: Code under measurement ----------- */
	if (is_valid) {
		printf("FATAL ERROR: benchmarked proof is valid\n");
	}
	sep_keys_clear(&pk, &sk);
	user_keys_clear(&upk, &usk);
	sep_sig_clear(&sig);
	show_proof_clear(&proof);
	poly_q_vec_d_clear(r[0]);
	poly_q_vec_d_clear(r[1]);
	poly_q_vec_d_clear(cmt);
	poly_qshow_vec_m1_clear(s1);
	for (i = 0; i < PARAM_D; i++) {
		poly_qshow_vec_k_clear(u_embed[i]);
		for (j = 0; j < PARAM_D; j++) {
			poly_qshow_mat_k_k_clear(A_embed[i][j]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + 0      ]);
			poly_qshow_mat_k_k_clear(Ds_embed[i][j + PARAM_D]);
		}
		for (j = 0; j < PARAM_D*PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(B_embed[i][j]);
		}
		for (j = 0; j < PARAM_K; j++) {
			poly_qshow_mat_k_k_clear(A3_embed[i][j]);
		}
		for (j = 0; j < PARAM_M; j++) {
			poly_qshow_mat_k_k_clear(D_embed[i][j]);
		}
	}
	return time;
}
