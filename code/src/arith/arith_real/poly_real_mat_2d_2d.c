#include "poly_real_mat_2d_2d.h"

/*************************************************
* Name:        poly_real_mat_2d_2d_init
*
* Description: Initialize polynomial matrix with 2*PARAM_D x 2*PARAM_D entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_real_mat_2d_2d res: polynomial matrix to be initialized
**************************************************/
void poly_real_mat_2d_2d_init(poly_real_mat_2d_2d res) {
  size_t i,j;
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_init(res->rows[i]->entries[j]);
    }
  }
}

/*************************************************
* Name:        poly_real_mat_2d_2d_clear
*
* Description: Clear polynomial matrix with 2*PARAM_D x 2*PARAM_D entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_real_mat_2d_2d res: polynomial matrix to be cleared
**************************************************/
void poly_real_mat_2d_2d_clear(poly_real_mat_2d_2d res) {
  size_t i,j;
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_clear(res->rows[i]->entries[j]);
    }
  }
}
