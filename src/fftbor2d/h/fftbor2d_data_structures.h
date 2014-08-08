#ifndef FFTBOR2D_DATA_STRUCTURES_H
#define FFTBOR2D_DATA_STRUCTURES_H

#include "vienna/externs.h"

#include <complex>
typedef std::complex<double> dcomplex;

typedef struct { // Variables are sorted by the order they get instantiated, do not change (structs suck in C++).
  char* sequence;
  char* structure_1;
  char* structure_2;
  char* energy_file;
  int   seq_length;
  int   precision;
  int   max_threads;
  char  format;
  short benchmark;
  short verbose;
} FFTBOR2D_PARAMS;

#define BASIC_FLAG 'B'
#define SIMPLE_FLAG 'S'
#define MATRIX_FLAG 'M'
#define CSV_FLAG 'C'
#define FORMAT(flag) (parameters.format == flag)
#define HUMANIZED_FORMAT                                      \
  (parameters.format == 'B' ? "basic (with headers)" :        \
    ((parameters.format == 'S' ? "simple (without headers)" : \
      ((parameters.format == 'M' ? "matrix" :                 \
        ((parameters.format == 'C' ? "CSV" : "N/A")))))))

typedef struct { // Variables are sorted by the order they get instantiated, do not change.
  char*     sequence;
  char*     structure_1;
  char*     structure_2;
  int       seq_length;
  double    RT;
  paramT*   vienna_params;
  int**     int_bp;
  char*     precision_format;
  short*    int_sequence;
  int**     can_base_pair;
  int***    num_base_pairs;
  int       bp_dist;
  int       row_length;
  int       run_length;
  int       num_roots;
  dcomplex* solutions;
  double    partition_function;
  dcomplex* roots_of_unity;
  double*   probabilities;
  int*      non_zero_indices;
  int       non_zero_count;
  int**     delta_table;
  int**     j_paired_to_0;
  int**     j_paired_to_1;
  double**  EZ;
  double**  EH;
  double**  EHM;
  double**  EMA;
  double**  EMB;
  double*** EIL;
  double*** EM1;
} FFTBOR2D_DATA;

typedef struct { // Variables are sorted by the order they get instantiated, do not change.
  dcomplex** Z;
  dcomplex** ZB;
  dcomplex** ZM;
  dcomplex** ZM1;
  dcomplex*  root_to_power;
} FFTBOR2D_THREADED_DATA;

#endif
