#ifndef SPECTRAL_DATA_STRUCTURES_H
#define SPECTRAL_DATA_STRUCTURES_H

#include "vienna_data_structures.h"

typedef struct {
  int verbose;
  char* sequence;
  char* energy_grid_file;
  char* start_structure;
  char* end_structure;
  double temperature;
  double start_time;
  double end_time;
  double step_size;
  int lonely_bp;
  int energy_cap;
  int eigen_only;
  int use_min;
} SPECTRAL_PARAMS;

typedef struct {
  double *values;       
  double *vectors;
  double *inverse_vectors;
  int length;
} EIGENSYSTEM;

#endif
