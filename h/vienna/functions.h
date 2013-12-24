#ifndef VIENNA_FUNCTIONS_H
#define VIENNA_FUNCTIONS_H

#include "data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

void set_model_details(model_detailsT*);
paramT* get_scaled_parameters(double, model_detailsT);
paramT* scale_parameters(void);
float fold_par(char*, char*, paramT*, int, int);
SOLUTION* subopt_par(char*, char*, paramT*, int, int, int, FILE*);

#ifdef __cplusplus
}
#endif

#endif