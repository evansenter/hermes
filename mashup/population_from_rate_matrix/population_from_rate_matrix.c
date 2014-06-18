#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include "vienna/externs.h"
#include "shared/libmulti_param_header.h"
#include "shared/libmfpt_header.h"
#include "shared/libpopulation_header.h"

int main(int argc, char** argv) {
  PARAM_CONTAINER* params;
  KLP_MATRIX klp_matrix;
  MFPT_PARAMS mfpt_params;
  POPULATION_PARAMS population_params;
  TRANSITION_MATRIX transition_matrix;

  char* subparams[] = { "mfpt", "population" };
  params            = split_args(argc, argv, subparams, 2);

  mfpt_params       = init_mfpt_params();
  mfpt_params.input = 0;

  parse_mfpt_args(&mfpt_params, params[0].argc, params[0].argv);

  klp_matrix        = klp_matrix_from_file(mfpt_params);
  transition_matrix = transition_matrix_from_klp_matrix(&klp_matrix, mfpt_params);

  population_params       = init_population_params();
  population_params.input = 0;

  parse_population_args(&population_params, params[1].argc, params[1].argv);

  population_from_row_ordered_transition_matrix(population_params, transition_matrix);

  return 0;
}
