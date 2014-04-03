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

extern int subopt_sorted;
extern double temperature;

int main(int argc, char** argv) {
  struct timeval full_start, full_stop, start, stop;
  SPECTRAL_PARAMS parameters;
  paramT*         vienna_params;
  model_detailsT  vienna_details;
  set_model_details(&vienna_details);
  vienna_details.noLP = !parameters.lonely_bp;
  vienna_params       = get_scaled_parameters(temperature, vienna_details);
  subopt_sorted       = 1;
  parameters          = init_spectral_params();
  parse_spectral_args(&parameters, argc, argv);

  if (parameters.sequence == NULL) {
    spectral_usage();
  }

  if (parameters.benchmark) {
    gettimeofday(&full_start, NULL);
    gettimeofday(&start, NULL);
  }

  int i, seq_length, energy_cap, num_structures = 0;
  char* sequence;
  char* empty_str;
  char* mfe_str;
  double mfe_energy;
  TRANSITION_MATRIX transition_matrix;
  EIGENSYSTEM eigensystem;
  SOLUTION* all_structures;

  if (!DESERIALIZING(parameters)) {
    sequence   = parameters.sequence;
    seq_length = strlen(sequence);
    empty_str  = malloc((seq_length + 1) * sizeof(char));
    mfe_str    = malloc((seq_length + 1) * sizeof(char));

    for (i = 0; i < seq_length; ++i) {
      empty_str[i] = '.';
    }

    empty_str[i]   = '\0';
    mfe_energy     = (double)fold_par(sequence, mfe_str, vienna_params, 0, 0);
    energy_cap     = parameters.energy_cap ? (int)round(parameters.energy_cap * 100) : 1000000;
    all_structures = subopt_par(sequence, empty_str, vienna_params, energy_cap, 0, 0, NULL);

    while (all_structures[num_structures].structure != NULL) {
      num_structures++;
    }

    find_key_structure_indices_in_structure_list(&parameters, all_structures, num_structures, empty_str, mfe_str);

    if (parameters.verbose) {
      printf("sequence:\t%s\n", sequence);
      printf("start:\t\t%s\t%+.2f kcal/mol\t(%d)\n", all_structures[parameters.start_index].structure, all_structures[parameters.start_index].energy, parameters.start_index);
      printf("stop:\t\t%s\t%+.2f kcal/mol\t(%d)\n", all_structures[parameters.end_index].structure, all_structures[parameters.end_index].energy, parameters.end_index);
      printf("mfe energy:\t%+.2f\n", mfe_energy);
      printf("energy cap:\t%+.2f kcal/mol above MFE\n", energy_cap / 100.);
      printf("num str:\t%d\n\n", num_structures);

      for (i = 0; i < num_structures; ++i) {
        printf("%d\t%s\t%+.2f\n", i, all_structures[i].structure, all_structures[i].energy);
      }

      printf("\n");
    }

    if (!SERIALIZING(parameters) && (parameters.start_index < 0 || parameters.end_index < 0)) {
      if (parameters.start_index < 0) {
        fprintf(stderr, "Error: the starting structure could not be found. Usually this means that you did not specify an explicit starting structure so the empty structure was used, but the energy band was not wide enough for RNAsubopt to sample it.\n");
      }

      if (parameters.end_index < 0) {
        fprintf(stderr, "Error: the ending structure could not be found. Usually this means that you did explicitly provided a suboptimal ending structure, but the energy band was not wide enough for RNAsubopt to sample it.\n");
      }

      abort();
    }

    if (parameters.benchmark) {
      gettimeofday(&stop, NULL);
      TIMING(start, stop, "initialization")
      gettimeofday(&start, NULL);
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

    print_population_proportion(parameters, eigensystem);

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
