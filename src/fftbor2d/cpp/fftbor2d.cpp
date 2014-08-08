#include <stdio.h>
#include "fftbor2d_params.h"
#include "fftbor2d_initializers.h"
#include "fftbor2d_functions.h"
#include "shared/timers.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
  FFTBOR2D_PARAMS parameters;
  parameters = init_fftbor2d_params();
  parse_fftbor2d_args(parameters, argc, argv, &fftbor2d_usage);
  
  START_ALL_TIMERS
  
  FFTBOR2D_DATA data;
  FFTBOR2D_THREADED_DATA* threaded_data;
  data          = init_fftbor2d_data(parameters);
  threaded_data = init_fftbor2d_threaded_data(parameters, data);
  TIMING("initialization")
  
  precalculate_energies(data);
  TIMING("precalculate_energies")
  
  evaluate_recursions_in_parallel(parameters, data, threaded_data);
  TIMING("evaluate_recursions_in_parallel")
  
  populate_remaining_roots(data);
  TIMING("populate_remaining_roots")
  
  solve_system(parameters, data);
  TIMING("solve_system")
  
  if (parameters.verbose) {
    print_fftbor2d_data(data);
  }
  
  print_output(parameters, data);
  TIMING("print_output")
  
  free_fftbor2d_threaded_data(threaded_data, parameters.max_threads);
  free_fftbor2d_data(data);
  free_fftbor2d_params(parameters);
  
  STOP_REMAINING_TIMERS
  
  return 0;
}
