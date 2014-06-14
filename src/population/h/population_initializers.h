#ifndef POPULATION_INITIALIZERS_H
#define POPULATION_INITIALIZERS_H

#include "population_data_structures.h"

extern double temperature;
extern int dangles;
extern int noLonelyPairs;
extern int subopt_sorted;

EIGENSYSTEM init_eigensystem(int);
void free_eigensystem(EIGENSYSTEM);
void print_eigensystem(const EIGENSYSTEM);
void init_vienna_global_params(const POPULATION_PARAMS);
paramT* init_vienna_params(const POPULATION_PARAMS);
pf_paramT* init_vienna_pf_params(const POPULATION_PARAMS);

#endif
