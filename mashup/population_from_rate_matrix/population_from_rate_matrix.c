#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libmulti_param_header.h"
#include "shared/libmfpt_header.h"
#include "shared/libpopulation_header.h"

void population_from_rate_matrix_usage(int);

char* subparams[] = { "mfpt", "population" };

int main(int argc, char** argv) {
  PARAM_CONTAINER* params;
  KLP_MATRIX klp_matrix;
  MFPT_PARAMS mfpt_params;
  POPULATION_PARAMS population_params;
  TRANSITION_MATRIX transition_matrix;


  params = split_args(argc, argv, subparams, 2);

  mfpt_params       = init_mfpt_params();
  mfpt_params.input = 0;

  parse_mfpt_args(&mfpt_params, params[0].argc, params[0].argv, &population_from_rate_matrix_usage);

  klp_matrix        = klp_matrix_from_file(mfpt_params, &population_from_rate_matrix_usage);
  transition_matrix = transition_matrix_from_klp_matrix(&klp_matrix, mfpt_params);

  population_params       = init_population_params();
  population_params.input = 0;

  parse_population_args(&population_params, params[1].argc, params[1].argv, &population_from_rate_matrix_usage);

  population_from_row_ordered_transition_matrix(population_params, transition_matrix);

  return 0;
}

void population_from_rate_matrix_usage(int _) {
  fprintf(stderr, "RateEq --mfpt-c <path/to/rate_matrix.csv> --population-a <starting index in matrix> --population-z <ending index in matrix> [additional options]\n\n");
  print_valid_multi_params_prefixes(subparams, 2);
  print_available_subflag_header("RNAmfpt", subparams[0], &mfpt_usage);
  print_available_subflag_header("RNAeq", subparams[1], &population_usage);
  abort();
}
