#include <stdlib.h>
#include "initializers.h"
#include "functions.h"
#include "vienna/functions.h"
#include "shared/libmfpt_header.h"

EIGENSYSTEM init_eigensystem(int length) {
  EIGENSYSTEM eigensystem = {
    .values          = calloc(length, sizeof(double)),
    .vectors         = calloc(length * length, sizeof(double*)),
    .inverse_vectors = calloc(length * length, sizeof(double*)),
    .length          = length
  };
  return eigensystem;
}

void free_eigensystem(EIGENSYSTEM eigensystem) {
  free(eigensystem.values);
  free(eigensystem.vectors);
  free(eigensystem.inverse_vectors);
}

void print_eigensystem(const EIGENSYSTEM eigensystem) {
  print_array("eigensystem.values", eigensystem.values, eigensystem.length);
  print_matrix("eigensystem.vectors", eigensystem.vectors, eigensystem.length);
  print_matrix("eigensystem.inverse_vectors", eigensystem.inverse_vectors, eigensystem.length);
}

void init_vienna_global_params(const POPULATION_PARAMS parameters) {
  dangles       = 0;
  noLonelyPairs = !parameters.lonely_bp;
  temperature   = parameters.temperature;
  subopt_sorted = 1;
}

paramT* init_vienna_params(const POPULATION_PARAMS parameters) {
  paramT* vienna_params;

  init_vienna_global_params(parameters);
  vienna_params = scale_parameters();

  return vienna_params;
}

pf_paramT* init_vienna_pf_params(const POPULATION_PARAMS parameters) {
  pf_paramT* vienna_pf_params;

  init_vienna_global_params(parameters);
  vienna_pf_params = get_scaled_pf_parameters();

  return vienna_pf_params;
}
