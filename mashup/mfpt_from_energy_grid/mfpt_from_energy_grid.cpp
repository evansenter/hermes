#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libmulti_param_header.h"
#include "shared/libfftbor2d_header.h"
#include "shared/libmfpt_header.h"

KLP_MATRIX convert_fftbor2d_output_to_klp_matrix(const FFTBOR2D_DATA);
void convert_klp_matrix_to_use_epsilon(KLP_MATRIX*, MFPT_PARAMS*, double);

int main(int argc, char** argv) {
  PARAM_CONTAINER* params;
  FFTBOR2D_PARAMS fftbor2d_params;
  FFTBOR2D_DATA fftbor2d_data;
  KLP_MATRIX klp_matrix;
  MFPT_PARAMS mfpt_params;
  double* transition_matrix;
  double mfpt;

  char* subparams[] = { "fftbor2d", "mfpt" };
  params            = split_args(argc, argv, subparams, 2);

  fftbor2d_params = init_fftbor2d_params();
  parse_fftbor2d_args(fftbor2d_params, params[0].argc, params[0].argv);
  fftbor2d_data = fftbor2d_from_params(fftbor2d_params);

  mfpt_params          = init_mfpt_params();
  mfpt_params.input    = 0;
  mfpt_params.epsilon  = 1e-8;
  mfpt_params.max_dist = fftbor2d_data.row_length;
  mfpt_params.bp_dist  = fftbor2d_data.bp_dist;

  parse_mfpt_args(&mfpt_params, params[1].argc, params[1].argv);

  klp_matrix        = convert_fftbor2d_output_to_klp_matrix(fftbor2d_data);
  transition_matrix = convert_klp_matrix_to_transition_matrix(&klp_matrix, &mfpt_params);
  mfpt              = compute_mfpt(&klp_matrix, mfpt_params, transition_matrix);

  printf("%f\n", mfpt);

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
