#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "klp_matrix_parser.h"
#include "klp_matrix_initializers.h"
#include "klp_matrix_functions.h"

KLP_MATRIX klp_matrix_from_file(const char* input_file, const KLP_PARAMS klp_params, void (*usage)()) {
  int line_count;
  KLP_MATRIX klp_matrix;
  
  line_count = count_lines(input_file);
  
  if (!line_count) {
    (*usage)();
  }
  
  klp_matrix = init_klp_matrix(line_count);
  populate_arrays(&klp_matrix, input_file, klp_params);
  
  return klp_matrix;
}

TRANSITION_MATRIX transition_matrix_from_klp_matrix_and_params(KLP_PARAMS* klp_params, KLP_MATRIX* klp_matrix) {
  TRANSITION_MATRIX transition_matrix;
  
  if (RUN_TYPE(klp_params->run_type, TRANSITION_INPUT_FLAG)) {
    // We already have a transition matrix, this is the easy case.
    transition_matrix = transition_matrix_from_klp_matrix(klp_matrix, MATRIX_TYPE(klp_params->rate_matrix));
  } else {
    // We have an energy grid, this requires converting the energy grid into a transition matrix data structure.
    transition_matrix = convert_klp_matrix_to_transition_matrix(klp_matrix, klp_params);
  }
  
  if (klp_params->output_only) {
    print_transition_matrix(transition_matrix);
    exit(0);
  } else {
    return transition_matrix;
  }
}

TRANSITION_MATRIX transition_matrix_from_klp_matrix(KLP_MATRIX* klp_matrix, char matrix_type) {
  int i, row_length = 0;
  TRANSITION_MATRIX transition_matrix;
  
  // We need to infer the dimensions of the transition matrix.
  for (i = 0; i < klp_matrix->length; ++i) {
    row_length = klp_matrix->k[i] > row_length ? klp_matrix->k[i] : row_length;
    row_length = klp_matrix->l[i] > row_length ? klp_matrix->l[i] : row_length;
  }
  
  // The transition matrix is 0-ordered, so we looked for the highest k, l position above and then we add one for the row length.
  // i.e. if the largest value we saw was 10, then we have rows going from 0..10 so row_length == 10 + 1, or 11.
  transition_matrix = init_transition_matrix(row_length + 1, matrix_type);
  
  for (i = 0; i < klp_matrix->length; ++i) {
    T_ROW_ORDER(transition_matrix, klp_matrix->k[i], klp_matrix->l[i]) = klp_matrix->p[i];
  }
  
  return transition_matrix;
}

int count_lines(const char* file_path) {
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

void populate_arrays(KLP_MATRIX* klp_matrix, const char* input_file, const KLP_PARAMS klp_params) {
  int i = 0;
  FILE* file = fopen(input_file, "r");
  char* token;
  char line[1024];
  
  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    klp_matrix->k[i] = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix->l[i] = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix->p[i] = atof(token);
    
    if (
      (klp_params.energy_based && !(RUN_TYPE(klp_params.run_type, TRANSITION_INPUT_FLAG) && klp_params.rate_matrix)) && 
      (klp_matrix->p[i] < 0 || klp_matrix->p[i] > 1)
    ) {
      fprintf(stderr, "Error: line number %d (0-indexed) in the input doesn't satisfy 0 <= probability (%+1.2f) <= 1. Did you forget the -e flag?\n\n", i, klp_matrix->p[i]);
      abort();
    }
    
    i++;
  }
  
  fclose(file);
}
