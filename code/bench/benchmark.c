#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <flint/flint.h>

#include "cpucycles.h"
#include "benchmark.h"

void start_timer(timer* t) {
  clock_gettime(CLOCK_MONOTONIC, &t->start);
  t->start_cycle = cpucycles();
}

double stop_timer(timer* t) {
  t->end_cycle = cpucycles();
  clock_gettime(CLOCK_MONOTONIC, &t->end);
  double start = (double) t->start.tv_sec * 1e9 + (double) t->start.tv_nsec;
  double end = (double) t->end.tv_sec * 1e9 + (double) t->end.tv_nsec;
  return (end - start) / 1e6;
}

uint64_t cpucycles_overhead(void) {
  uint64_t t0, t1, overhead = -1LL;
  unsigned int i;

  for(i=0;i<100000;i++) {
    t0 = cpucycles();
    __asm__ volatile ("");
    t1 = cpucycles();
    if(t1 - t0 < overhead)
      overhead = t1 - t0;
  }

  return overhead;
}

double clock_gettime_overhead(void) {
  unsigned int i;
  struct timespec start, end;
  double overhead = 1e9;

  for(i=0;i<100000;i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    __asm__ volatile ("");
    clock_gettime(CLOCK_MONOTONIC, &end);
    double p = ((double) start.tv_sec * 1e9 + (double) start.tv_nsec - ((double) end.tv_sec * 1e9 + (double) end.tv_nsec)) / 1e6;
    if(p < overhead)
      overhead = p;
  }

  return overhead;
}

static int cmp_double(const void *lhs, const void *rhs) {
	double x = *((const double*) lhs);
	double y = *((const double*) rhs);
	if (x < y) {
		return -1;
	} else if (x > y) {
		return 1;
	} else {
		return 0;
	}
}

typedef struct {
	double mean;
	double med;
	double min;
	double max;
	double var;
} stat_t;

static void get_stats(stat_t *stats, double *data, size_t len) {
  qsort(data, len, sizeof(data[0]), cmp_double);
  stats->min = data[0];
  stats->max = data[len - 1];
  stats->med = (data[len / 2] + data[(len - 1) / 2]) / 2.0;

  // Welford's online algorithm (from Wikipedia) to avoid critical floating point errors
  stats->mean = 0;
  stats->var = 0.0;
  for (size_t i = 0; i < len; ++i) {
	  double delta = data[i] - stats->mean;
	  stats->mean += delta / (i + 1);
	  //double delta2 = data[i] - stats->mean;
	  //stats->var += delta * delta2;
  }
  assert(len >= 2 && "len cannot be < 2 in for online mean-variance computation");
  stats->var /= (len - 1);
}

void benchmark(const char* name, size_t num_iterations, double (*fun) (timer*)) {
  double timings[num_iterations];
  double cycles[num_iterations];
  timer t;
  stat_t stats_time;
  stat_t stats_cycles;

  static double ovh_ms = -1;
  static uint64_t ovh_cyc = -1UL;

  if (ovh_ms == -1)
  {
    ovh_ms = clock_gettime_overhead();
  }
  if (ovh_cyc == -1UL)
  {
    ovh_cyc = cpucycles_overhead();
  }

  for (size_t i = 0; i < num_iterations; ++i) {
    timings[i] = fun(&t) - ovh_ms;
    cycles[i] = (double) (t.end_cycle - t.start_cycle - ovh_cyc);
    flint_cleanup();
  }
  printf("%s:\n", name);
  get_stats(&stats_time, timings, num_iterations);
  get_stats(&stats_cycles, cycles, num_iterations);
  printf("\tmean = %.8f ms, med = %.8f ms, min = %.8f ms, max = %.8f ms\n\tmean = %.2f cycles, med = %.0f cycles, min = %.0f cycles, max = %.0f cycles\n",
	 stats_time.mean,
	 stats_time.med,
	 stats_time.min,
	 stats_time.max,
	 stats_cycles.mean,
	 stats_cycles.med,
	 stats_cycles.min,
	 stats_cycles.max);
}
