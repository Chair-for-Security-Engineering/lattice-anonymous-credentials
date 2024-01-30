#include "four_squares.h"

/*************************************************
* Name:        two_squares [static]
*
* Description: Finds a two square decomposition of
*              the positive integer n if it exists
*
* Arguments:   - uint64_t *res: pointer to the output array (allocated)
*              - uint64_t n: integer to be decomposed
*
* Returns 1 if decomposition succeeds (0 if it fails)
**************************************************/
static int _two_squares(uint64_t res[2], uint64_t n) {
  if (n == 0) {
    res[0] = 0;
    res[1] = 0;
    return 1;
  }
  
  uint32_t fac = 0;
  while (n % 4  == 0) { // Division by the largest power of 4
    n >>= 2;
    fac += 1;
  }

  if (n % 4 == 3) {
    return 0; // n not a sum of two squares
  } 

  uint64_t i,ii,j,jj,nn;
  // n = 1 mod 4 ==> different parity
  if (n % 4 == 1) {
    j = (uint64_t) sqrt((double) n);
    jj = j * j;
    i = 0;
    ii = 0;
    while (ii <= jj) {
      nn = n - ii;
      while (jj > nn) {
        j -= 1;
        jj = j * j;
      }
      if (jj == nn) {
        res[0] = i<<fac;
        res[1] = j<<fac;
        return 1;
      }
      i += 1;
      ii = i*i;
    } 
  }
  // n = 2 mod 4 ==> both odd
  else {
    i = 1;
    ii = 1;
    j = (uint64_t) sqrt((double) n);
    j += 1 - j%2;
    jj = j * j;
    while (ii <= jj) {
      nn = n - ii;
      while (jj > nn) {
        j -= 2;
        jj = j * j;
      }
      if (jj == nn) {
        res[0] = i<<fac;
        res[1] = j<<fac;
        return 1;
      }
      i += 2;
      ii = i*i;
    } 
  }
  return 0;
}

/*************************************************
* Name:        three_squares [static]
*
* Description: Finds a three square decomposition of
*              the positive integer n if it exists
*
* Arguments:   - uint64_t *res: pointer to the output array (allocated)
*              - uint64_t n: integer to be decomposed
*
* Returns 1 if decomposition succeeds (0 if it fails)
**************************************************/
static int _three_squares(uint64_t res[3], uint64_t n) {
  if (n == 0) {
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    return 1;
  }

  uint32_t fac = 0;
  while (n % 4  == 0) {
    n >>= 2;
    fac += 1;
  }

  if (n % 8 == 7) {
    return 0;
  }

  uint64_t j;
  j = (uint64_t) sqrt((double) n);
  while (!_two_squares(res, n - j*j)) {
    j -= 1;
  }
  res[0] <<= fac;
  res[1] <<= fac;
  res[2] = j<<fac;
  return 1;
}

/*************************************************
* Name:        four_squares
*
* Description: Finds a four square decomposition of
*              the positive integer n
*
* Arguments:   - uint64_t *res: pointer to the output array (allocated)
*              - uint64_t n: integer to be decomposed
*
* Returns 1 if decomposition succeeds (0 if it fails)
**************************************************/
void four_squares(uint64_t res[4], uint64_t n) {
  if (n == 0) {
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else {
    uint32_t fac = 0;
    while (n % 4  == 0) {
      n >>= 2;
      fac += 1;
    }

    uint64_t j;
    j = (uint64_t) sqrt((double) n);
    while (!_three_squares(res, n - j*j)) {
      j -= 1;
    }
    res[0] <<= fac;
    res[1] <<= fac;
    res[2] <<= fac;
    res[3] = j<<fac;
  }
}