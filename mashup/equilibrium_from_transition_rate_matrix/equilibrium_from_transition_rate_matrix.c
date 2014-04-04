#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libmulti_param_header.h"
#include "shared/libmfpt_header.h"
#include "shared/libspectral_header.h"
#include "shared/libtpl_header.h"

int main(int argc, char** argv) {
  PARAM_CONTAINER* params;
  KLP_MATRIX klp_matrix;
  MFPT_PARAMS mfpt_params;
  SPECTRAL_PARAMS spectral_params;
  TRANSITION_MATRIX transition_matrix;

  char* subparams[] = { "mfpt", "spectral" };
  params            = split_args(argc, argv, subparams, 2);

  mfpt_params       = init_mfpt_params();
  mfpt_params.input = 0;

  parse_mfpt_args(&mfpt_params, params[0].argc, params[0].argv);

  klp_matrix        = klp_matrix_from_file(mfpt_params);
  transition_matrix = transition_matrix_from_klp_matrix(&klp_matrix, mfpt_params);

  spectral_params       = init_spectral_params();
  spectral_params.input = 0;

  parse_spectral_args(&spectral_params, params[1].argc, params[1].argv);

  population_proportion_from_row_ordered_transition_matrix(spectral_params, transition_matrix);

  return 0;
}
