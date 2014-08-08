#include <stdio.h>
#include <stdlib.h>
#include "mfpt_params.h"
#include "mfpt_functions.h"
#include "shared/libklp_matrix_header.h"
#include "shared/constants.h"
#include "shared/timers.h"

int main(int argc, char** argv) {
  double mfpt;
  KLP_PARAMS klp_params;
  MFPT_PARAMS parameters;
  KLP_MATRIX klp_matrix;
  TRANSITION_MATRIX transition_matrix;
  klp_params = init_klp_matrix_params();
  parameters = init_mfpt_params();
  parse_klp_matrix_args(&klp_params, argc, argv, &mfpt_usage);
  parse_mfpt_args(&klp_params, &parameters, argc, argv, &mfpt_usage);
  
  START_ALL_TIMERS
  
  klp_matrix = klp_matrix_from_file(parameters.input_file, parameters.input && !klp_params.energy_based, &mfpt_usage);
  TIMING("klp_matrix_from_file")
    
#ifdef SUPER_HEAVY_DEBUG
  print_klp_matrix(klp_matrix);
#endif
  
  if (RUN_TYPE(klp_params, TRANSITION_INPUT_FLAG)) {
    // We already have a transition matrix, this is the easy case. Just need to find MFPT.
    transition_matrix = transition_matrix_from_klp_matrix(&klp_matrix, MATRIX_TYPE(klp_params));
    TIMING("transition_matrix_from_klp_matrix")
  } else {
    // We have an energy grid, this requires converting the energy grid into a transition matrix data structure before finding MFPT.
    transition_matrix = convert_klp_matrix_to_transition_matrix(&klp_matrix, &klp_params);
    TIMING("convert_klp_matrix_to_transition_matrix")
  }
  
#if SUPER_HEAVY_DEBUG
  print_transition_matrix_with_klp_positions(klp_matrix, transition_matrix);
#endif
  
  mfpt = compute_mfpt(transition_matrix, klp_params, parameters);
  TIMING("compute_mfpt")
    
  printf("%+.8f\n", mfpt);
  free_transition_matrix(transition_matrix);
  free_klp_matrix(klp_matrix);
  
  STOP_REMAINING_TIMERS
  
  return 0;
}
