#ifndef SPECTRAL_FUNCTIONS_H
#define SPECTRAL_FUNCTIONS_H

#include "data_structures.h"

void population_proportion_from_row_ordered_transition_matrix(const SPECTRAL_PARAMS, double*, int);
double* convert_structures_to_transition_matrix(const SOLUTION*, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
void invert_matrix(EIGENSYSTEM*);
double probability_at_time(const EIGENSYSTEM, double, int, int);
void find_key_structure_indices_in_structure_list(SPECTRAL_PARAMS*, const SOLUTION*, int, char*, char*);
void print_population_proportion(const SPECTRAL_PARAMS, const EIGENSYSTEM);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);

#endif
