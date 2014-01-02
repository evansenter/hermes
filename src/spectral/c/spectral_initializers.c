#include <stdlib.h>
#include "functions.h"
#include "shared/libmfpt_header.h"

EIGENSYSTEM init_eigensystem(int length) {
  EIGENSYSTEM eigensystem = {
    .values          = calloc(length, sizeof(double)),
    .vectors         = init_transition_matrix(length),
    .inverse_vectors = init_transition_matrix(length),
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
