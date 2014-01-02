#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "constants.h"
#include "vienna/functions.h"
#include "functions.h"
#include "params.h"
#include "initializers.h"

#define TIMING(start, stop, task) printf("[benchmarking] %8.2f\ttime in ms for %s\n", (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0), task);

extern double temperature;
SPECTRAL_PARAMS parameters;

int main(int argc, char** argv) {
  struct timeval full_start, full_stop, start, stop;
  paramT*        vienna_params;
  model_detailsT vienna_details;
  set_model_details(&vienna_details);
  parameters          = parse_spectral_args(argc, argv);
  vienna_details.noLP = !parameters.lonely_bp;
  vienna_params       = get_scaled_parameters(temperature, vienna_details);
  
  if (parameters.sequence == NULL) {
    spectral_usage();
  }
  
  if (parameters.benchmark) {
    gettimeofday(&full_start, NULL);
    gettimeofday(&start, NULL);
  }
  
  int i, seq_length, energy_cap, num_structures = 0;
  double mfe_energy;
  double* transition_matrix;
  EIGENSYSTEM eigensystem;
  char* sequence  = parameters.sequence;
  seq_length      = strlen(sequence);
  char* empty_str = malloc((seq_length + 1) * sizeof(char));
  char* mfe_str   = malloc((seq_length + 1) * sizeof(char));
  
  for (i = 0; i < seq_length; ++i) {
    empty_str[i] = '.';
  }
  
  empty_str[i]             = '\0';
  mfe_energy               = (double)fold_par(sequence, mfe_str, vienna_params, 0, 0);
  energy_cap               = parameters.energy_cap ? (int)(2 * abs(mfe_energy * 100)) : 1000000;
  SOLUTION* all_structures = subopt_par(sequence, empty_str, vienna_params, energy_cap, 0, 0, NULL);
  
  while (all_structures[num_structures].structure != NULL) {
    num_structures++;
  }
  
  find_key_structure_indices_in_structure_list(&parameters, all_structures, num_structures, empty_str, mfe_str);
  
  if (parameters.verbose) {
    printf("sequence:\t%s\n", sequence);
    printf("start:\t\t%s\t%+.2f kcal/mol\t(%d)\n", parameters.start_structure, all_structures[parameters.start_index].energy, parameters.start_index);
    printf("stop:\t\t%s\t%+.2f kcal/mol\t(%d)\n", parameters.end_structure, all_structures[parameters.end_index].energy, parameters.end_index);
    printf("energy cap:\t%+.2f kcal/mol above MFE\n", energy_cap / 100.);
    printf("num str:\t%d\n", num_structures);
  }
  
  #ifdef SUPER_HEAVY_DEBUG
  
  for (i = 0; i < num_structures; ++i) {
    printf("%s\t%+.2f\t%d\n", all_structures[i].structure, all_structures[i].energy, i);
  }
  
  printf("\n");
  #endif
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  }
  
  transition_matrix = convert_structures_to_transition_matrix(all_structures, num_structures);
  #ifdef INSANE_DEBUG
  print_matrix("transition_matrix", transition_matrix, num_structures);
  #endif
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "convert_structures_to_transition_matrix")
    gettimeofday(&start, NULL);
  }
  
  eigensystem = convert_transition_matrix_to_eigenvectors(transition_matrix, num_structures);
  
  if (parameters.eigen_only) {
    for (i = 0; i < seq_length; ++i) {
      printf("%+.8f\n", eigensystem.values[i]);
    }
  } else {
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "convert_transition_matrix_to_eigenvectors")
      gettimeofday(&start, NULL);
    }
    
    invert_matrix(&eigensystem);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "invert_matrix")
      gettimeofday(&start, NULL);
    }
    
    #ifdef INSANE_DEBUG
    print_eigensystem(eigensystem);
    #endif
    
    print_population_proportion(parameters, eigensystem);
  }
  
  free_eigensystem(eigensystem);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "generate_population_points")
    gettimeofday(&full_stop, NULL);
    TIMING(full_start, full_stop, "total")
  }
  
  return 0;
}
