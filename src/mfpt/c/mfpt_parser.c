#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mfpt_params.h"
#include "mfpt_parser.h"
#include "mfpt_initializers.h"

KLP_MATRIX klp_matrix_from_file(const MFPT_PARAMS parameters, void (*usage)(int)) {
  int line_count;
  KLP_MATRIX klp_matrix;

  line_count = count_lines(parameters.input_file);

  if (!line_count) {
    (*usage)(0);
  }

  klp_matrix = init_klp_matrix(line_count);
  populate_arrays(&klp_matrix, parameters);

  return klp_matrix;
}

TRANSITION_MATRIX transition_matrix_from_klp_matrix(KLP_MATRIX* klp_matrix, const MFPT_PARAMS parameters) {
  int i, row_length = 0;
  TRANSITION_MATRIX transition_matrix;

  // We need to infer the dimensions of the transition matrix.
  for (i = 0; i < klp_matrix->length; ++i) {
    row_length = klp_matrix->k[i] > row_length ? klp_matrix->k[i] : row_length;
    row_length = klp_matrix->l[i] > row_length ? klp_matrix->l[i] : row_length;
  }

  // The transition matrix is 0-ordered, so we looked for the highest k, l position above and then we add one for the row length.
  // i.e. if the largest value we saw was 10, then we have rows going from 0..10 so row_length == 10 + 1, or 11.
  transition_matrix = init_transition_matrix(row_length + 1, MATRIX_TYPE(parameters));

  for (i = 0; i < klp_matrix->length; ++i) {
    T_ROW_ORDER(transition_matrix, klp_matrix->k[i], klp_matrix->l[i]) = klp_matrix->p[i];
  }

  return transition_matrix;
}

int count_lines(char* file_path) {
  FILE* file = fopen(file_path, "r");
  int c;
  int line_count = 0;

  if (file == NULL) {
    return 0;
  }

  while ((c = fgetc(file)) != EOF) {
    if (c == '\n') {
      line_count++;
    }
  }

  fclose(file);
  return line_count;
}

void populate_arrays(KLP_MATRIX* klp_matrix, const MFPT_PARAMS parameters) {
  int i = 0;
  FILE* file = fopen(parameters.input_file, "r");
  char* token;
  char line[1024];

  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    klp_matrix->k[i] = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix->l[i] = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix->p[i] = atof(token);

    if (parameters.input && !parameters.energy_based && (klp_matrix->p[i] < 0 || klp_matrix->p[i] > 1)) {
      fprintf(stderr, "Error: line number %d (0-indexed) in the input doesn't satisfy 0 <= probability (%+1.2f) <= 1. Did you forget the -e flag?\n\n", i, klp_matrix->p[i]);
      mfpt_usage(0);
    }

    i++;
  }

  fclose(file);
}
