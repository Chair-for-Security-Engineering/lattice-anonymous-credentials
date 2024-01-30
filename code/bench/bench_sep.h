#ifndef BENCH_SEP_H
#define BENCH_SEP_H

#include "benchmark.h"

double sep_keygen_bench(timer* t);

double sep_sign_bench(timer* t);

double sep_verify_valid_bench(timer* t);

double sep_verify_invalid_bench(timer* t);

#endif /* BENCH_SEP_H */

