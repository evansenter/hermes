#ifndef POPULATION_DATA_STRUCTURES_H
#define POPULATION_DATA_STRUCTURES_H

#include "shared/libklp_matrix_header.h"
#include "vienna/externs.h"

typedef struct {
  short verbose;
  char* input_file;
  char* sequence;
  char* start_structure;
  char* end_structure;
  char* filename;
  int serialize;
  double temperature;
  double start_time;
  double end_time;
  double step_size;
  short equilibrium;
  short soft_bounds;
  double epsilon;
  double delta;
  int window_size;
  int all_subpop_for_eq;
  short lonely_bp;
  double energy_cap;
  short eigen_only;
  short input;
  short benchmark;
} POPULATION_PARAMS;

typedef struct {
  double* values;
  double* vectors;
  double* inverse_vectors;
  int length;
} EIGENSYSTEM;

#define SERIALIZING(params) (params.serialize == 1)
#define DESERIALIZING(params) (params.serialize == -1)
#define E_ROW_ORDER(m, i, j, n) ((m)[((i) * (n)) + (j)])
#define E_COL_ORDER(m, i, j, n) ((m)[((j) * (n)) + (i)])

#endif
