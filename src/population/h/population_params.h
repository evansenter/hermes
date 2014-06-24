#ifndef POPULATION_PARAMS_H
#define POPULATION_PARAMS_H

#include "population_data_structures.h"

POPULATION_PARAMS init_population_params();
void parse_population_args(POPULATION_PARAMS*, int, char**, void (*)(int));
int population_error_handling(const POPULATION_PARAMS);
void debug_population_parameters(const POPULATION_PARAMS);
void population_usage(int);

#endif
