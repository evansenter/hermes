#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "population_constants.h"
#include "population_functions.h"
#include "population_params.h"
#include "population_initializers.h"
#include "shared/libklp_matrix_header.h"
#include "shared/timers.h"

int main(int argc, char** argv) {
  KLP_PARAMS klp_params;
  POPULATION_PARAMS parameters;
  SOLUTION* all_structures;
  KLP_MATRIX klp_matrix;
  TRANSITION_MATRIX transition_matrix;
  EIGENSYSTEM eigensystem;
  
  int i, num_structures  = 0;
  klp_params             = init_klp_matrix_params();
  klp_params.rate_matrix = 1;
  parameters             = init_population_params();
  parse_klp_matrix_args(&klp_params, argc, argv, &population_usage);
  parse_population_args(&klp_params, &parameters, argc, argv, &population_usage);
  
  START_ALL_TIMERS
  
  // If we're not deserializing an input binary eigensystem (which doesn't necessarily mean we're serializing something), we need to setup the eigensystem
  // from the user input.
  if (!DESERIALIZING(parameters)) {
    if (parameters.input_file) {
      klp_matrix        = klp_matrix_from_file(parameters.input_file, klp_params, &population_usage);
      transition_matrix = transpose_matrix(transition_matrix_from_klp_matrix_and_params(&klp_params, &klp_matrix));
      TIMING("csv initialization")
    } else if (parameters.sequence) {
      ensure_key_structures_and_energies_assigned(&parameters);
      all_structures = sample_structures(parameters);
      
      while (all_structures[num_structures].structure != NULL) {
        num_structures++;
      }
      
      find_key_structure_indices_in_structure_list(&klp_params, parameters, all_structures, num_structures);
      
      if (parameters.verbose) {
        printf("sequence:\t%s\n",                      parameters.sequence);
        printf("start:\t\t%s\t%+.2f kcal/mol\t(%d)\n", all_structures[klp_params.start_state].structure, all_structures[klp_params.start_state].energy, klp_params.start_state);
        printf("stop:\t\t%s\t%+.2f kcal/mol\t(%d)\n",  all_structures[klp_params.end_state].structure, all_structures[klp_params.end_state].energy, klp_params.end_state);
        printf("num str:\t%d\n\n",                     num_structures);
        
        for (i = 0; i < num_structures; ++i) {
          printf("%d\t%s\t%+.2f\n", i, all_structures[i].structure, all_structures[i].energy);
        }
        
        printf("\n");
      }
      
      TIMING("subopt initialization")
      
      transition_matrix = convert_structures_to_transition_matrix(all_structures, num_structures);
      TIMING("convert_structures_to_transition_matrix")
    }
  }
  
  if (DESERIALIZING(parameters)) {
    eigensystem = deserialize_eigensystem(parameters);
    TIMING("deserialization")
  } else {
    eigensystem = convert_transition_matrix_to_eigenvectors(transition_matrix);
    TIMING("convert_transition_matrix_to_eigenvectors")
    
    if (!parameters.eigen_only) {
      invert_matrix(&eigensystem);
      TIMING("invert_matrix")
    }
  }
  
  if (parameters.eigen_only) {
    print_eigenvalues(eigensystem);
    TIMING_WITHOUT_RESTART("print_eigenvalues")
  } else if (parameters.equilibrium) {
    print_equilibrium(eigensystem, klp_params, parameters);
    TIMING_WITHOUT_RESTART("print_equilibrium")
  } else if (SERIALIZING(parameters)) {
    serialize_eigensystem(eigensystem, parameters);
    TIMING_WITHOUT_RESTART("serialize_eigensystem")
  } else {
    print_population_proportion(eigensystem, klp_params, parameters);
    TIMING_WITHOUT_RESTART("print_population_proportion")
  }
  
  free_eigensystem(eigensystem);
  
  STOP_REMAINING_TIMERS
  
  return 0;
}
