#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libmulti_param_header.h"
#include "shared/libmfpt_header.h"
#include "shared/libspectral_header.h"

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

  spectral_params             = init_spectral_params();
  spectral_params.input       = 0;
  spectral_params.step_size   = 1e-2;
  spectral_params.equilibrium = 1e-4;

  parse_spectral_args(&spectral_params, params[1].argc, params[1].argv);

  // This will print everything:
  // EIGENSYSTEM eigensystem;
  // eigensystem = eigensystem_from_row_ordered_transition_matrix(transition_matrix);
  // print_population_proportion(eigensystem, spectral_params);
  // print_equilibrium(eigensystem, spectral_params);

  equilibrium_from_row_ordered_transition_matrix(spectral_params, transition_matrix);

  return 0;
}
