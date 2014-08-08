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
  
  klp_matrix = klp_matrix_from_file(parameters.input_file, klp_params, &mfpt_usage);
  TIMING("klp_matrix_from_file")
    
#ifdef SUPER_HEAVY_DEBUG
  print_klp_matrix(klp_matrix);
#endif
  
  transition_matrix = transition_matrix_from_klp_matrix_and_params(&klp_params, &klp_matrix);
  TIMING("transition_matrix_from_klp_matrix_and_params")
  
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
