#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <stdint.h>
#include <time.h>

typedef struct {
  struct timespec start;
  struct timespec end;
  uint64_t start_cycle;
  uint64_t end_cycle;
} timer;

void start_timer(timer* t);

double stop_timer(timer* t);

void benchmark(const char* name, size_t num_iterations, double (*fun) (timer*));

#endif /* BENCHMARK_H */

