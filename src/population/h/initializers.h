#ifndef POPULATION_INITIALIZERS_H
#define POPULATION_INITIALIZERS_H

#include "data_structures.h"

EIGENSYSTEM init_eigensystem(int);
void free_eigensystem(EIGENSYSTEM);
void print_eigensystem(const EIGENSYSTEM);

#endif
