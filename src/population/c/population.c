#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "population_constants.h"
#include "population_functions.h"
#include "population_params.h"
#include "population_initializers.h"
#include "shared/libklp_matrix_header.h"

#define TIMING(start, stop, task) printf("[benchmarking] %8.2f\ttime in ms for %s\n", (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0), task);

int main(int argc, char** argv) {
  struct timeval full_start, full_stop, start, stop;
  KLP_PARAMS klp_params;
  POPULATION_PARAMS parameters;
  SOLUTION* all_structures;
  KLP_MATRIX klp_matrix;
  TRANSITION_MATRIX transition_matrix;
  EIGENSYSTEM eigensystem;
  
  int i, num_structures = 0;
  
  klp_params = init_klp_matrix_params();
  parameters = init_population_params();
  parse_klp_matrix_args(&klp_params, argc, argv, &population_usage);
  parse_population_args(&klp_params, &parameters, argc, argv, &population_usage);
  
  if (parameters.benchmark) {
    gettimeofday(&full_start, NULL);
    gettimeofday(&start, NULL);
  }
  
  // If we're not deserializing an input binary eigensystem (which doesn't necessarily mean we're serializing something), we need to setup the eigensystem
  // from the user input.
  if (!DESERIALIZING(parameters)) {
    if (parameters.input_file) {
      klp_matrix        = klp_matrix_from_file(parameters.input_file, parameters.input && !klp_params.energy_based, &population_usage);
      transition_matrix = convert_klp_matrix_to_transition_matrix(&klp_matrix, &klp_params);
      
      if (parameters.benchmark) {
        gettimeofday(&stop, NULL);
        TIMING(start, stop, "csv initialization")
        gettimeofday(&start, NULL);
      }
    } else if (parameters.sequence) {
      ensure_key_structures_and_energies_assigned(&parameters);
    
      all_structures = sample_structures(parameters);
    
      while (all_structures[num_structures].structure != NULL) {
        num_structures++;
      }
    
      find_key_structure_indices_in_structure_list(&klp_params, parameters, all_structures, num_structures);
    
      if (parameters.verbose) {
        printf("sequence:\t%s\n", parameters.sequence);
        printf("start:\t\t%s\t%+.2f kcal/mol\t(%d)\n", all_structures[klp_params.start_state].structure, all_structures[klp_params.start_state].energy, klp_params.start_state);
        printf("stop:\t\t%s\t%+.2f kcal/mol\t(%d)\n", all_structures[klp_params.end_state].structure, all_structures[klp_params.end_state].energy, klp_params.end_state);
        printf("num str:\t%d\n\n", num_structures);
      
        for (i = 0; i < num_structures; ++i) {
          printf("%d\t%s\t%+.2f\n", i, all_structures[i].structure, all_structures[i].energy);
        }
      
        printf("\n");
      }
    
      if (parameters.benchmark) {
        gettimeofday(&stop, NULL);
        TIMING(start, stop, "subopt initialization")
        gettimeofday(&start, NULL);
      }
    }
  }
  
  if (DESERIALIZING(parameters)) {
    eigensystem = deserialize_eigensystem(parameters);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "deserialization")
      gettimeofday(&start, NULL);
    }
  } else {
    transition_matrix = convert_structures_to_transition_matrix(all_structures, num_structures);
#ifdef INSANE_DEBUG
    print_matrix("transition_matrix", transition_matrix.matrix, transition_matrix.row_length);
#endif
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "convert_structures_to_transition_matrix")
      gettimeofday(&start, NULL);
    }
    
    eigensystem = convert_transition_matrix_to_eigenvectors(transition_matrix);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "convert_transition_matrix_to_eigenvectors")
      gettimeofday(&start, NULL);
    }
    
    if (!parameters.eigen_only) {
      invert_matrix(&eigensystem);
      
      if (parameters.benchmark) {
        gettimeofday(&stop, NULL);
        TIMING(start, stop, "invert_matrix")
        gettimeofday(&start, NULL);
      }
    }
  }
  
  if (parameters.eigen_only) {
    print_eigenvalues(eigensystem);
    exit(0);
  } else if (parameters.equilibrium) {
    print_equilibrium(eigensystem, klp_params, parameters);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "estimate_equilibrium")
    }
  } else if (SERIALIZING(parameters)) {
    serialize_eigensystem(eigensystem, parameters);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "serializing")
    }
  } else {
#ifdef INSANE_DEBUG
    print_eigensystem(eigensystem);
#endif
    
    print_population_proportion(eigensystem, klp_params, parameters);
    
    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "generate_population_points")
    }
  }
  
  free_eigensystem(eigensystem);
  
  if (parameters.benchmark) {
    gettimeofday(&full_stop, NULL);
    TIMING(full_start, full_stop, "total")
  }
  
  return 0;
}
