#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "klp_matrix_initializers.h"
#include "klp_matrix_functions.h"

TRANSITION_MATRIX transpose_matrix(TRANSITION_MATRIX matrix) {
  int i, j;
  TRANSITION_MATRIX transposed_matrix = init_transition_matrix(matrix.row_length, matrix.type);
  
  for (i = 0; i < matrix.row_length; ++i) {
    for (j = 0; j < matrix.row_length; ++j) {
      T_COL_ORDER(transposed_matrix, i, j) = T_ROW_ORDER(matrix, i, j);
    }
  }
  
  free_transition_matrix(matrix);
  return transposed_matrix;
}

TRANSITION_MATRIX inverse(TRANSITION_MATRIX transition_matrix) {
  int i, j, signum;
  gsl_matrix* matrix_to_invert = gsl_matrix_alloc(transition_matrix.row_length, transition_matrix.row_length);
  gsl_matrix* inversion_matrix = gsl_matrix_alloc(transition_matrix.row_length, transition_matrix.row_length);
  gsl_permutation* permutation = gsl_permutation_alloc(transition_matrix.row_length);
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      gsl_matrix_set(matrix_to_invert, i, j, T_ROW_ORDER(transition_matrix, i, j));
    }
  }
  
  gsl_linalg_LU_decomp(matrix_to_invert, permutation, &signum);
  gsl_linalg_LU_invert(matrix_to_invert, permutation, inversion_matrix);
  
  for (i = 0; i < transition_matrix.row_length; ++i) {
    for (j = 0; j < transition_matrix.row_length; ++j) {
      T_ROW_ORDER(transition_matrix, i, j) = gsl_matrix_get(inversion_matrix, i, j);
    }
  }
  
  gsl_matrix_free(matrix_to_invert);
  gsl_matrix_free(inversion_matrix);
  gsl_permutation_free(permutation);
  
  return transition_matrix;
}
