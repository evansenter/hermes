#include <stdio.h>
#include <stdlib.h>
#include "initializers.h"

KLP_MATRIX init_klp_matrix(int length) {
  KLP_MATRIX klp_matrix = {
    .k          = (int*)malloc(length * sizeof(int)),
    .l          = (int*)malloc(length * sizeof(int)),
    .p          = (double*)malloc(length * sizeof(double)),
    .length     = length
  };
  return klp_matrix;
}

void free_klp_matrix(KLP_MATRIX klp_matrix) {
  free(klp_matrix.k);
  free(klp_matrix.l);
  free(klp_matrix.p);
}

void print_klp_matrix(const KLP_MATRIX klp_matrix) {
  int i;
  printf("\nk/l/p matrix (index, k, l, p):\n");
  
  for (i = 0; i < klp_matrix.length; ++i) {
    printf("%d\t%d\t%d\t%+.8f\n", i, klp_matrix.k[i], klp_matrix.l[i], klp_matrix.p[i]);
  }
  
  printf("\n");
}

TRANSITION_MATRIX init_transition_matrix(int length, char type) {
  if (!(type == 'P' || type == 'R')) {
    fprintf(stderr, "Error: init_transition_matrix type was %c, must be one of (P)robability, (R)ate.\n", type);
    abort();
  }
  
  TRANSITION_MATRIX transition_matrix = {
    .matrix     = calloc(length * length, sizeof(double)),
    .row_length = length,
    .type       = type
  };
  
  return transition_matrix;
}

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

void free_transition_matrix(TRANSITION_MATRIX transition_matrix) {
  free(transition_matrix.matrix);
}

void print_transition_matrix(const KLP_MATRIX klp_matrix, const TRANSITION_MATRIX transition_matrix, const MFPT_PARAMS parameters) {
  int i, j;
  printf("Transition matrix:\n");
  printf("(from)\t(to)\tp(to | from)\n");
  
  for (i = 0; i < klp_matrix.length; ++i) {
    for (j = 0; j < klp_matrix.length; ++j) {
      if (RUN_TYPE(parameters, TRANSITION_INPUT_FLAG)) {
        printf(
          "%d\t=>\t%d\t%.8f\n",
          i,
          j,
          T_ROW_ORDER(transition_matrix, i, j)
        );
      } else {
        printf(
          "(%d, %d)\t=>\t(%d, %d)\t%.8f\n",
          klp_matrix.k[i],
          klp_matrix.l[i],
          klp_matrix.k[j],
          klp_matrix.l[j],
          T_ROW_ORDER(transition_matrix, i, j)
        );
      }
    }
  }
  
  printf("\n");
}
