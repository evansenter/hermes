#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libmulti_param_header.h"
#include "shared/libfftbor2d_header.h"
#include "shared/libklp_matrix_header.h"
#include "shared/libpopulation_header.h"

KLP_MATRIX convert_fftbor2d_output_to_klp_matrix(const FFTBOR2D_DATA);
void population_from_fftbor2d_usage();

char* subparams[] = { "fftbor2d", "population" };

int main(int argc, char** argv) {
  PARAM_CONTAINER*  params;
  FFTBOR2D_PARAMS   fftbor2d_params;
  KLP_PARAMS        klp_params;
  POPULATION_PARAMS population_params;
  FFTBOR2D_DATA     fftbor2d_data;
  KLP_MATRIX        klp_matrix;
  TRANSITION_MATRIX transition_matrix;
  
  params = split_args(argc, argv, subparams, 2);
  
  fftbor2d_params = init_fftbor2d_params();
  parse_fftbor2d_args(fftbor2d_params, params[0].argc, params[0].argv, &population_from_fftbor2d_usage);
  fftbor2d_data = fftbor2d_from_params(fftbor2d_params);
  
  klp_params             = init_klp_matrix_params();
  klp_params.rate_matrix = 1;
  klp_params.epsilon     = 1e-8;
  klp_params.max_dist    = fftbor2d_data.row_length;
  klp_params.bp_dist     = fftbor2d_data.bp_dist;
  parse_klp_matrix_args(&klp_params, params[1].argc, params[1].argv, &population_from_fftbor2d_usage);
  
  klp_matrix        = convert_fftbor2d_output_to_klp_matrix(fftbor2d_data);
  transition_matrix = convert_klp_matrix_to_transition_matrix(&klp_matrix, &klp_params);
  
  population_params                 = init_population_params();
  population_params.sequence        = fftbor2d_data.sequence;
  population_params.start_structure = fftbor2d_data.structure_1;
  population_params.end_structure   = fftbor2d_data.structure_2;
  population_params.temperature     = temperature;
  parse_population_args(&klp_params, &population_params, params[1].argc, params[1].argv, &population_from_fftbor2d_usage);
  
  population_from_row_ordered_transition_matrix(klp_params, population_params, transition_matrix);
  
  return 0;
}

KLP_MATRIX convert_fftbor2d_output_to_klp_matrix(const FFTBOR2D_DATA fftbor2d_data) {
  int i;
  KLP_MATRIX klp_matrix;
  
  klp_matrix = init_klp_matrix(fftbor2d_data.non_zero_count);
  
  for (i = 0; i < klp_matrix.length; ++i) {
    klp_matrix.k[i] = fftbor2d_data.non_zero_indices[i] / fftbor2d_data.row_length;
    klp_matrix.l[i] = fftbor2d_data.non_zero_indices[i] % fftbor2d_data.row_length;
    klp_matrix.p[i] = fftbor2d_data.probabilities[fftbor2d_data.non_zero_indices[i]];
  }
  
  return klp_matrix;
}

void population_from_fftbor2d_usage() {
  fprintf(stderr, "FFTeq --fftbor2d-i <sequence> --fftbor2d-j <structure_1> --fftbor2d-k <structure_2> [additional options]\n\n");
  print_valid_multi_params_prefixes(subparams, 2);
  print_available_subflags("FFTbor2D",          subparams[0], &fftbor2d_flags);
  print_available_subflags("transition matrix", subparams[1], &klp_matrix_flags);
  print_available_subflags("RNAeq",             subparams[1], &population_flags);
  abort();
}
