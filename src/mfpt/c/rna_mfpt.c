#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "parser.h"
#include "initializers.h"
#include "functions.h"
#include "constants.h"

int main(int argc, char** argv) {
  double mfpt;
  MFPT_PARAMS parameters;
  KLP_MATRIX klp_matrix;
  TRANSITION_MATRIX transition_matrix;
  parameters = init_mfpt_params();
  parse_mfpt_args(&parameters, argc, argv);

  klp_matrix = klp_matrix_from_file(parameters);
#ifdef SUPER_HEAVY_DEBUG
  print_klp_matrix(klp_matrix);
#endif

  if (RUN_TYPE(parameters, TRANSITION_INPUT_FLAG)) {
    // We already have a transition matrix, this is the easy case. Just need to find MFPT.
    transition_matrix = transition_matrix_from_klp_matrix(&klp_matrix, parameters);
  } else {
    // We have an energy grid, this requires converting the energy grid into a transition matrix data structure before finding MFPT.
    transition_matrix = convert_klp_matrix_to_transition_matrix(&klp_matrix, &parameters);
  }

#if SUPER_HEAVY_DEBUG
  print_transition_matrix(klp_matrix, transition_matrix, parameters);
#endif
  mfpt = compute_mfpt(&klp_matrix, parameters, transition_matrix);
  printf("%+.8f\n", mfpt);
  free_transition_matrix(transition_matrix);
  free_klp_matrix(klp_matrix);
  return 0;
}
