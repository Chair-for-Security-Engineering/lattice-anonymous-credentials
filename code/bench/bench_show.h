#ifndef BENCH_SHOW_H
#define BENCH_SHOW_H

#include "benchmark.h"

double show_user_embed_bench(timer* t);
double show_user_prove_bench(timer* t);
double show_user_verify_valid_bench(timer* t);
double show_user_verify_invalid_bench(timer* t);

#endif /* BENCH_SHOW_H */
