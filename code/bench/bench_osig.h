#ifndef BENCH_OSIG_H
#define BENCH_OSIG_H

#include "benchmark.h"

double osig_user_commit_bench(timer* t);
double osig_signer_sign_commitment_bench(timer* t);
double osig_user_sign_complete_bench(timer* t);
double osig_user_verify_valid_bench(timer* t);
double osig_user_verify_invalid_bench(timer* t);
double osig_user_keygen_bench(timer* t);
double osig_user_embed_bench(timer* t);
double osig_user_prove_bench(timer* t);
double osig_signer_verify_valid_bench(timer* t);
double osig_signer_verify_invalid_bench(timer* t);

#endif /* BENCH_OSIG */

