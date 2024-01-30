/* This file is an adapted version from https://github.com/lucasprabel/module_gaussian_lattice/blob/main/Module_BFRS/random.h */

#include <inttypes.h>
#include <x86intrin.h>

/*
	Code from random_aesni.c
*/

//public API
void random_init(void);

int64_t SampleZ(double c, double sigma);