#ifndef MFPT_DATA_STRUCTURES_H
#define MFPT_DATA_STRUCTURES_H

typedef struct {
  char* input_file;
  short all_mfpt;
  short input;
  short benchmark;
  short verbose;
} MFPT_PARAMS;

#include "shared/libklp_matrix_header.h"

#endif
