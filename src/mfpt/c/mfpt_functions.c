#include <stdio.h>
#include <stdlib.h>
#include "mfpt_functions.h"
#include "shared/libklp_matrix_header.h"
#include "shared/constants.h"

double compute_mfpt(const TRANSITION_MATRIX transition_matrix, const KLP_PARAMS klp_params, const MFPT_PARAMS parameters) {
  int i, j, x, y, start_pointer;
  double mfpt_from_start;
  TRANSITION_MATRIX inversion_matrix;
  double* mfpt;
  
  if (klp_params.start_state < 0 || klp_params.end_state < 0) {
    if (klp_params.start_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
#endif
      return -1;
    }
    
    if (klp_params.end_state < 0) {
#ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the stopping state.\n");
#endif
      return -2;
    }
  }
  
  // If start_index > end_index, we need to shift to the left by one because the end_index row / column is being removed.
  start_pointer    = klp_params.start_state - (klp_params.start_state > klp_params.end_state ? 1 : 0);
  inversion_matrix = init_transition_matrix(transition_matrix.row_length - 1, transition_matrix.type);
  mfpt             = calloc(inversion_matrix.row_length, sizeof(double));
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      if (i != klp_params.end_state && j != klp_params.end_state) {
        x = (i > klp_params.end_state ? i - 1 : i);
        y = (j > klp_params.end_state ? j - 1 : j);
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (end_index) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversion_matrix and i, j are used for indexing into transition_matrix.
        T_ROW_ORDER(inversion_matrix, x, y) = \
                                              (i == j ? 1 - T_ROW_ORDER(transition_matrix, i, j) : -T_ROW_ORDER(transition_matrix, i, j));
      }
    }
  }
  
  inversion_matrix = inverse(inversion_matrix);
  
  for (i = 0; i < inversion_matrix.row_length; ++i) {
    for (j = 0; j < inversion_matrix.row_length; ++j) {
      mfpt[i] += T_ROW_ORDER(inversion_matrix, i, j);
    }
    
    if (parameters.all_mfpt) {
      // The business with this i < end_index stuff is inorder to ensure that the output MFPT indices are representative of the input data, since we reduce the dimension of the matrix by 1.
      printf("%d\t%+.8f\n", i < klp_params.end_state ? i : i + 1, mfpt[i]);
    }
  }
  
  mfpt_from_start = mfpt[start_pointer];
  
  free(mfpt);
  free_transition_matrix(inversion_matrix);
  return mfpt_from_start;
}
