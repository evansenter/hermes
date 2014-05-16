#ifndef VIENNA_FUNCTIONS_H
#define VIENNA_FUNCTIONS_H

#include <stdio.h>
#include "data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

paramT* scale_parameters(void);
pf_paramT *get_scaled_pf_parameters(void);
float pf_fold_par(const char*, char*, pf_paramT*, int, int, int);
float fold_par(char*, char*, paramT*, int, int);
float energy_of_struct_par(const char*, const char*, paramT*, int);
SOLUTION* subopt_par(char*, char*, paramT*, int, int, int, FILE*);

#ifdef __cplusplus
}
#endif

#endif