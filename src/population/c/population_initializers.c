#include <stdlib.h>
#include "functions.h"
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
