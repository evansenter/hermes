#ifndef MFPT_DATA_STRUCTURES_H
#define MFPT_DATA_STRUCTURES_H

typedef struct {
  char* input_file;
  short all_mfpt;
  short input;
  short verbose;
} MFPT_PARAMS;

#include "shared/libklp_matrix_header.h"

typedef double(*transition_probability)(const KLP_MATRIX, const double*, int, int, short);

#endif
