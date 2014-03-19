#ifndef SPECTRAL_DATA_STRUCTURES_H
#define SPECTRAL_DATA_STRUCTURES_H

#include "vienna/data_structures.h"

typedef struct {
  short verbose;
  char* sequence;
  char* start_structure;
  char* end_structure;
  int start_index;
  int end_index;
  double temperature;
  double start_time;
  double end_time;
  double step_size;
  short lonely_bp;
  double energy_cap;
  short eigen_only;
  short benchmark;
} SPECTRAL_PARAMS;

typedef struct {
  double* values;
  double* vectors;
  double* inverse_vectors;
  int length;
} EIGENSYSTEM;

#endif
