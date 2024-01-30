/* This file is an adapted version from https://github.com/lucasprabel/module_gaussian_lattice/blob/main/Module_BFRS/random.c */

#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <x86intrin.h>

#include "random.h"
#include "randombytes.h"

static void random_bytes(uint8_t * restrict);

/*************************************************
* Name:        uniform_int_distribution
*
* Description: Sample an integer uniformly in [0,n] using the uniform
*              distribution over [0, 2^32 - 1] provided by random_bytes
*
* Arguments:   - uint32_t n: upper bound on uniform element
* 
* Returns an integer uniformly distributed in [0,n]
**************************************************/
static uint32_t uniform_int_distribution(uint32_t n) {
  uint32_t scaling = (UINT32_MAX) / (n + 1);
  uint32_t past = (n + 1) * scaling;
  uint32_t r_data[4];

  while(true) {
    random_bytes((uint8_t *) r_data);
    for(int i = 0 ; i < 4 ; ++i) {
      uint32_t r = r_data[i];
      if(r < past) {
        return r / scaling;
      }
    }
  }
}

/*
	Code from random_aesni.c
*/

//macros
#define DO_ENC_BLOCK(m,k) \
	do{\
        m = _mm_xor_si128       (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[ 1]); \
        m = _mm_aesenc_si128    (m, k[ 2]); \
        m = _mm_aesenc_si128    (m, k[ 3]); \
        m = _mm_aesenc_si128    (m, k[ 4]); \
        m = _mm_aesenc_si128    (m, k[ 5]); \
        __auto_type m5 = m;\
        m = _mm_aesenc_si128    (m, k[ 6]); \
        m = _mm_aesenc_si128    (m, k[ 7]); \
        m = _mm_aesenc_si128    (m, k[ 8]); \
        m = _mm_aesenc_si128    (m, k[ 9]); \
        m = _mm_aesenclast_si128(m, k[10]);\
        m = _mm_xor_si128(m, m5);\
    }while(0)

#define DO_DEC_BLOCK(m,k) \
	do{\
        m = _mm_xor_si128       (m, k[10+0]); \
        m = _mm_aesdec_si128    (m, k[10+1]); \
        m = _mm_aesdec_si128    (m, k[10+2]); \
        m = _mm_aesdec_si128    (m, k[10+3]); \
        m = _mm_aesdec_si128    (m, k[10+4]); \
        m = _mm_aesdec_si128    (m, k[10+5]); \
        m = _mm_aesdec_si128    (m, k[10+6]); \
        m = _mm_aesdec_si128    (m, k[10+7]); \
        m = _mm_aesdec_si128    (m, k[10+8]); \
        m = _mm_aesdec_si128    (m, k[10+9]); \
        m = _mm_aesdeclast_si128(m, k[0]);\
    }while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

//the expanded key
static __m128i key_schedule[20];
static uint64_t ctr = 0;

static __m128i aes_128_key_expansion(__m128i key, __m128i keygened){
	keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	return _mm_xor_si128(key, keygened);
}

static void aes128_load_key(const int8_t * enc_key){
    key_schedule[0] = _mm_loadu_si128((const __m128i*) enc_key);
	key_schedule[1]  = AES_128_key_exp(key_schedule[0], 0x01);
	key_schedule[2]  = AES_128_key_exp(key_schedule[1], 0x02);
	key_schedule[3]  = AES_128_key_exp(key_schedule[2], 0x04);
	key_schedule[4]  = AES_128_key_exp(key_schedule[3], 0x08);
	key_schedule[5]  = AES_128_key_exp(key_schedule[4], 0x10);
	key_schedule[6]  = AES_128_key_exp(key_schedule[5], 0x20);
	key_schedule[7]  = AES_128_key_exp(key_schedule[6], 0x40);
	key_schedule[8]  = AES_128_key_exp(key_schedule[7], 0x80);
	key_schedule[9]  = AES_128_key_exp(key_schedule[8], 0x1B);
	key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);

	// generate decryption keys in reverse order.
    // k[10] is shared by last encryption and first decryption rounds
    // k[0] is shared by first encryption round and last decryption round (and is the original user key)
    // For some implementation reasons, decryption key schedule is NOT the encryption key schedule in reverse order
	key_schedule[11] = _mm_aesimc_si128(key_schedule[9]);
	key_schedule[12] = _mm_aesimc_si128(key_schedule[8]);
	key_schedule[13] = _mm_aesimc_si128(key_schedule[7]);
	key_schedule[14] = _mm_aesimc_si128(key_schedule[6]);
	key_schedule[15] = _mm_aesimc_si128(key_schedule[5]);
	key_schedule[16] = _mm_aesimc_si128(key_schedule[4]);
	key_schedule[17] = _mm_aesimc_si128(key_schedule[3]);
	key_schedule[18] = _mm_aesimc_si128(key_schedule[2]);
	key_schedule[19] = _mm_aesimc_si128(key_schedule[1]);
}

static void aes128_enc(const int8_t * restrict plainText, int8_t * restrict cipherText){
    __m128i m = _mm_loadu_si128((const __m128i *) plainText);

    DO_ENC_BLOCK(m,key_schedule);

    _mm_storeu_si128((__m128i *) cipherText, m);
}

//public API
void random_init(void) {
	// Seed the AES key
    uint8_t seed[16];
    randombytes(seed, sizeof(uint8_t)*16);
    aes128_load_key((int8_t *) seed); // don't know why the AES key is signed but whatever
    ctr = 0;
}

static void random_bytes(uint8_t * restrict data) {
    uint8_t plaintext[16] = {0};
    *((uint64_t *) plaintext) = ctr++;
    aes128_enc((int8_t *) plaintext, (int8_t *) data); // let's cast them to signed because why not
}

/*
	Code from exp_aes.cpp
*/
static double algorithm_EA(uint64_t * n) {
  const double ln2 = 0.6931471805599453; // log(2);
  const double a = 5.713363152645423; //(4+3*sqrt(2))*log(2);
  const double b = 3.414213562373095; // 2+sqrt(2);
  const double c = -1.6734053240284923; // -(1+sqrt(2))*log(2);
  const double p = 0.9802581434685472; //sqrt(2)*log(2);
  const double A = 5.6005707569738075; //a*p;
  const double B = 3.3468106480569846; // b*p;
  //const double z = c + 2; // unused ?
  //const double r = 1.010089582001449; // 8*exp(-z)*b*(b-1)*log(2)/(a*a); // unused ?
  //const double h = r - p; // unused ?
  const double D = 0.08578643762690495; // 1/(b*b);
  const double H = 0.0026106723602095233;// h*D/p;

  uint8_t data[16];
  random_bytes(data);

  const unsigned long long int II = *((uint64_t *) data);
  const int j = __builtin_ctzll(II);
  const double G = c + j*ln2;

  const double U = *((uint16_t *) (data + 8)) / (UINT16_MAX*1.0);
  random_bytes(data);

  (*n)++;
  if (U <= p) {
    return G + A/(B - U);
  }
  //unsigned long long k = 0; // unused ?
  while(true) {
    for(unsigned i = 0; i < 4; i++) {
      (*n)++;
      const double U1 = *((uint16_t *) (data + i*2)) / (UINT16_MAX*1.0);
      const double U2 = *((uint16_t *) (data + i*2 + 8)) / (UINT16_MAX*1.0);

      const double bU1 = b - U1;
      const double Y = a/bU1;
      const double L = (U2*H + D)*(bU1)*(bU1);
      const double Z = Y + c;

      const double LZ  = L - 1 + Z;
      const double ZZ  = Z*Z;
      const bool cond1 = LZ <= 0;
      const bool cond2 = LZ - ZZ/2 + ZZ*Z/6 <= 0;
      const bool cond3 = LZ - ZZ <= 0;

      if (cond1 || cond2 || (cond3 && (L <= exp(-Z)))) {
        return G + Y;
      }
//      if(!cond1 && !cond2 && cond3) {
//        (*n)++;
//      }
    }
    random_bytes(data);
  }
}

/*
	Code from algoF_aes.cpp
*/

static int algorithmF(double mu, double sigma) {
	const double sigma_floor = floor(sigma);
	const unsigned sigma_max = (unsigned) (((sigma_floor == sigma) && (sigma_floor != 0)) ? sigma - 1 : (sigma_floor));


	uint8_t randomData[16];
	random_bytes(randomData);
	while(true) {
		for(unsigned i = 0; i < 16; i++) {
			uint64_t n = 0;
			unsigned k = (unsigned) ceil(2*algorithm_EA(&n)) - 1;

			int s = (randomData[i] > 127) ? -1 : 1;
			//printf("\ts = %d\n", s);
			unsigned j = uniform_int_distribution(sigma_max);

			int64_t i0 = (int64_t) ceil(sigma*k + s*mu);
   			double x0 = (i0 - (sigma*k + s*mu))/sigma;
    		double x = x0 + j/sigma;

    		unsigned k1 = (unsigned) ceil(2*algorithm_EA(&n)) - 1;
    		double z = algorithm_EA(&n);

    		bool cond1 = k1 >= k*(k-1);
    		bool cond2 = x < 1;
    		bool cond3a = x != 0;
    		bool cond3b = k != 0;
    		bool cond3c = s == 1;
    		bool cond3 = cond3a || cond3b || cond3c;
    		bool cond4 = z > 0.5*x*(2*k+x);

    		if (cond1 && cond2 && cond3 && cond4) {
    			return s*(i0 + ((int) j));
    		}
		}
		random_bytes(randomData);
	}
}

/*************************************************
* Name:        sampleZ
*
* Description: Sample an integer from the discrete Gaussian
*              distribution with density proportional to 
*              exp(-pi * |x-c|^2 / sigma^2)
*
* Arguments:   - double c: center for the Gaussian distribution
*              - double sigma: Gaussian parameter (std dev = sigma / sqrt(2*pi))
* 
* Returns a Gaussian sample from D_{Z, sigma, c}
**************************************************/
int64_t SampleZ(double c, double sigma) {
  return algorithmF(c, sigma/2.50662827463100); // 2.50662827463100 = sqrt(2*pi)
}
