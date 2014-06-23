#ifndef VIENNA_FUNCTIONS_H
#define VIENNA_FUNCTIONS_H

#include <stdio.h>
#include "data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

void read_parameter_file(const char*);
paramT* scale_parameters(void);
pf_paramT* get_scaled_pf_parameters(void);
float pf_fold_par(const char*, char*, pf_paramT*, int, int, int);
float fold_par(char*, char*, paramT*, int, int);
SOLUTION* subopt_par(char*, char*, paramT*, int, int, int, FILE*);
int* get_iindx(unsigned int seq_length);
unsigned int* maximumMatchingConstraint(const char* sequence, short* vienna_bp);

#ifdef __cplusplus
}
#endif

#endif