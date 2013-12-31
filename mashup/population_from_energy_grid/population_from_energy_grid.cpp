#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shared/libfftbor2d_header.h"
#include "shared/libspectral_header.h"

int main(int argc, char** argv) {
  FFTBOR2D_PARAMS fftbor2d_params;
  FFTBOR2D_DATA fftbor2d_data;
  FFTBOR2D_THREADED_DATA* threaded_data;
  
  fftbor2d_params = parse_fftbor2d_args(argc, argv);
  fftbor2d_data   = init_fftbor2d_data(fftbor2d_params);
  threaded_data   = init_fftbor2d_threaded_data(fftbor2d_params, fftbor2d_data);
  
  precalculate_energies(fftbor2d_data);
  evaluate_recursions_in_parallel(fftbor2d_params, fftbor2d_data, threaded_data);
  populate_remaining_roots(fftbor2d_data);
  solve_system(fftbor2d_params, fftbor2d_data);
  print_output(fftbor2d_params, fftbor2d_data);
  
  return 0;
}
